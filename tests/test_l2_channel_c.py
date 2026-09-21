"""Tests for Channel C's guardrails.

Every gene, drug and cell line here is invented, and none of these tests opens the 5 GB
LINCS release: they are about the promises the channel makes, not about its arithmetic,
which is tested on constructed cases in ``test_l2_connectivity.py``.

The promises are: the proxy caveat reaches every row, the conventional quality threshold
is reported rather than silently applied, no curated gene set is substituted when the
knockdown is missing, the null is the one that can fail, and the channel refuses rather
than nominating what it cannot defend.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.l2_channels import channel_c_signature as channel, lincs

REPO_SOURCE = None


def _settings(**overrides) -> dict:
    return channel.settings_for({"l2": {"channel_c": overrides}}) if overrides \
        else channel.settings_for({})


def _meta(gene="INVENTEDGENE", cell_lines=("CELLA", "CELLB"), quality=0.5,
          pert_type=lincs.KNOCKDOWN_TYPE):
    return [{"sig_id": f"{cell}_{gene}", "pert_iname": gene, "pert_type": pert_type,
             "cell_id": cell, "distil_cc_q75": str(quality), "distil_ss": "2.0"}
            for cell in cell_lines]


class _FakeMatrix:
    """Enough of an h5py file for :func:`lincs.read_block` on a handful of rows."""

    def __init__(self, values):
        self._values = np.asarray(values, dtype=np.float64)

    def __getitem__(self, key):
        if key == lincs.MATRIX_PATH:
            return self._values
        raise KeyError(key)


# ------------------------------------------------------------ settings


def test_defaults_are_complete_and_typed() -> None:
    settings = _settings()
    assert set(settings) == set(channel.DEFAULTS)
    assert isinstance(settings["query_size"], int)
    assert isinstance(settings["approved_only"], bool)


def test_an_unknown_setting_is_refused() -> None:
    """A typo must not silently leave a default in place on a channel this delicate."""
    with pytest.raises(ValueError, match="unknown setting"):
        _settings(querysize=50)


def test_a_null_too_small_to_calibrate_is_refused() -> None:
    with pytest.raises(ValueError, match="without a null is not reported"):
        _settings(n_control_genes=1)


def test_deflation_is_on_by_default() -> None:
    """Measured: uncorrected, the ranking does not beat random gene sets."""
    assert channel.DEFAULTS["deflate_generic_axis"] is True


# ------------------------------------------------------------ the proxy


def test_the_proxy_pools_every_cell_line_and_reports_quality() -> None:
    values = np.array([[3.0, 1.0, -1.0, -3.0], [1.0, 3.0, -3.0, -1.0]])
    consensus, report = channel.build_proxy(
        _FakeMatrix(values), ["CELLA_INVENTEDGENE", "CELLB_INVENTEDGENE"],
        _meta(), "INVENTEDGENE", [0, 1, 2, 3], _settings())
    assert np.allclose(consensus, values.mean(axis=0))
    assert report["cell_lines"] == ["CELLA", "CELLB"]
    assert report["cross_cell_line_concordance"]["n_pairs"] == 1
    assert report["proxy_type"] == channel.PROXY_TYPE


def test_the_conventional_threshold_is_reported_not_applied() -> None:
    """Both halves matter: the weak signatures are kept, and the count is stated."""
    weak = 0.05
    assert weak < channel.CONVENTIONAL_CC_Q75
    values = np.array([[3.0, 1.0, -1.0, -3.0], [1.0, 3.0, -3.0, -1.0]])
    _, report = channel.build_proxy(
        _FakeMatrix(values), ["CELLA_INVENTEDGENE", "CELLB_INVENTEDGENE"],
        _meta(quality=weak), "INVENTEDGENE", [0, 1, 2, 3], _settings())
    assert len(report["signatures"]) == 2          # kept
    assert report["n_below_conventional_cc_q75"] == 2   # and counted
    assert report["dropped"] == []


def test_raising_the_threshold_does_filter() -> None:
    """The setting is real, not decorative."""
    values = np.array([[3.0, 1.0, -1.0, -3.0], [1.0, 3.0, -3.0, -1.0]])
    with pytest.raises(ValueError, match="were dropped"):
        channel.build_proxy(
            _FakeMatrix(values), ["CELLA_INVENTEDGENE", "CELLB_INVENTEDGENE"],
            _meta(quality=0.05), "INVENTEDGENE", [0, 1, 2, 3],
            _settings(min_replicate_correlation=0.2))


def test_a_missing_knockdown_refuses_rather_than_substituting_a_gene_set() -> None:
    """The scaffold offered a curated aneuploidy set as a fallback. Taking it silently
    would leave the channel reporting a score whose provenance no longer matches its name.
    """
    with pytest.raises(ValueError, match="does not fall back"):
        channel.build_proxy(_FakeMatrix(np.zeros((1, 4))), ["CELLA_OTHER"],
                            _meta(gene="OTHERGENE"), "INVENTEDGENE", [0, 1, 2, 3],
                            _settings())


def test_a_compound_signature_is_not_mistaken_for_a_knockdown() -> None:
    with pytest.raises(ValueError, match="no trt_sh.cgs consensus signature"):
        channel.build_proxy(_FakeMatrix(np.zeros((1, 4))), ["CELLA_INVENTEDGENE"],
                            _meta(pert_type=lincs.COMPOUND_TYPE), "INVENTEDGENE",
                            [0, 1, 2, 3], _settings())


# ------------------------------------------------------------ deflation


def test_deflation_removes_the_generic_axis() -> None:
    rng = np.random.default_rng(5)
    z_scores = rng.normal(size=(50, 30)) + 2.0      # a strong shared component
    consensus = z_scores.mean(axis=0) + rng.normal(scale=0.1, size=30)
    deflated, report = channel.deflate_generic_axis(consensus, z_scores)
    assert abs(report["r_after"]) < abs(report["r_before"])
    # Exactly zero, not merely small: the axis and the reported statistic are the same
    # geometry, which is why generic_axis() centres.
    assert report["r_after"] == pytest.approx(0.0, abs=1e-6)
    assert report["applied"] is True


def test_deflation_leaves_the_signatures_alone() -> None:
    """The correction is applied to the query. Each compound's measured response is a
    fact about the experiment and is not edited to suit the question."""
    rng = np.random.default_rng(6)
    z_scores = rng.normal(size=(20, 12)) + 1.0
    before = z_scores.copy()
    channel.deflate_generic_axis(z_scores.mean(axis=0), z_scores)
    assert np.array_equal(z_scores, before)


def test_projecting_out_twice_changes_nothing_further() -> None:
    rng = np.random.default_rng(7)
    z_scores = rng.normal(size=(20, 12)) + 1.0
    axis = channel.generic_axis(z_scores)
    once = channel.project_out(rng.normal(size=12), axis)
    assert np.allclose(once, channel.project_out(once, axis), atol=1e-12)


# ------------------------------------------------------------ the null


def test_the_null_is_other_genes_not_random_gene_sets() -> None:
    """Recorded because two random-gene-set nulls were tried and both passed everything.

    The failures are kept in the output so the null is not quietly swapped for an easier
    one by a later session that finds this one inconvenient.
    """
    _, report = channel.calibrate_depth(
        np.array([-1.0, -0.5, 0.0]), {"GENEA": np.array([-0.1, 0.0, 0.1])}, _settings())
    assert "another gene" in report["method"].lower()
    assert len(report["superseded_nulls"]) == 2


def test_a_head_no_better_than_the_controls_yields_depth_zero() -> None:
    """The null has to be able to fail, or it is decoration."""
    controls = {f"GENE{i}": np.array([-2.0, -1.5, -1.0]) for i in range(20)}
    depth, report = channel.calibrate_depth(np.array([-0.5, -0.4, -0.3]), controls,
                                            _settings(top_n=1))
    assert depth == 0
    assert all(not rung["passes"] for rung in report["ladder"])


def test_a_head_better_than_every_control_passes() -> None:
    controls = {f"GENE{i}": np.array([-0.1, 0.0, 0.1]) for i in range(20)}
    depth, report = channel.calibrate_depth(np.array([-9.0, -8.0, -7.0]), controls,
                                            _settings(top_n=1))
    assert depth == 1
    # The ladder rounds p to four places for the record it writes.
    assert report["ladder"][0]["p"] == pytest.approx(1 / 21, abs=1e-4)


def test_the_control_genes_exclude_the_sac_panel_and_the_target() -> None:
    """A spindle-checkpoint gene is not an unrelated control for another one."""
    meta = []
    for gene in ("BUB1B", "BUB1", "CEP57", "TRIP13", "BUB3", "CEP192",
                 *(f"INVENTED{i}" for i in range(10))):
        meta.extend(_meta(gene=gene, cell_lines=("CELLA", "CELLB")))
    controls = channel.control_genes(meta, "BUB1B", _settings(n_control_genes=5), 42)
    assert len(controls) == 5
    assert not {c for c in controls} & {"BUB1B", "BUB1", "CEP57", "TRIP13", "BUB3",
                                        "CEP192"}


def test_control_genes_are_reproducible_from_the_seed() -> None:
    meta = [row for i in range(30)
            for row in _meta(gene=f"INVENTED{i}", cell_lines=("CELLA", "CELLB"))]
    meta.extend(_meta(gene="TARGET", cell_lines=("CELLA", "CELLB")))
    settings = _settings(n_control_genes=6)
    assert channel.control_genes(meta, "TARGET", settings, 42) == \
        channel.control_genes(meta, "TARGET", settings, 42)


def test_too_few_eligible_control_genes_is_refused() -> None:
    meta = _meta(gene="TARGET") + _meta(gene="INVENTED1")
    with pytest.raises(ValueError, match="fewer than the"):
        channel.control_genes(meta, "TARGET", _settings(n_control_genes=40), 42)


# ------------------------------------------------------------ the caveat


def test_the_caveat_says_what_the_score_is_not() -> None:
    caveat = channel.PROXY_CAVEAT
    assert "NOT THE PATIENT" in caveat
    assert "no RNA from the subject" in caveat
    for word in ("knockdown", "cell lines", "hypothesis"):
        assert word in caveat


def test_the_caveat_is_a_column_the_writer_cannot_omit() -> None:
    """It travels per row, not in a footnote a downstream join could drop."""
    source = (__import__("pathlib").Path(channel.__file__)).read_text(encoding="utf-8")
    header_line = [line for line in source.splitlines() if '"proxy_caveat"' in line]
    assert header_line, "proxy_caveat must appear in the candidates.tsv header"
    assert "PROXY_CAVEAT," in source, "and be appended to every row"


def test_every_scaffold_field_still_reaches_the_table() -> None:
    """The stub's guardrail named four fields; none may be dropped by a later edit."""
    source = (__import__("pathlib").Path(channel.__file__)).read_text(encoding="utf-8")
    for field in ("proxy_signature_id", "proxy_type", "cell_line", "proxy_caveat"):
        assert f'"{field}"' in source


def test_the_channel_declares_no_patient_data_path() -> None:
    source = (__import__("pathlib").Path(channel.__file__)).read_text(encoding="utf-8")
    assert "data_dir" not in source.replace("never opens ``data_dir``", "")
