"""Tests for L5's publication boundary and rendering.

Every identifier, drug and phenotype term here is **invented**. The point of this file is
the boundary: L5 writes the files a judge opens, so the checks that must hold are about
what may not appear in them, and about the rendering rules that keep a shortlist from
reading as advice.
"""

from __future__ import annotations

import pytest

from src.l5_report import publish_guard, render
from src.l5_report.run import _contradiction_text, _verdicts_by_candidate


# ------------------------------------------------------------ the publication boundary


@pytest.mark.parametrize("text, description", [
    ("A feature list included HP:1234567 in the draft.", "an HPO term id"),
    ("Variant at 7:117559590 in the panel.", "a genomic coordinate"),
    ("The change c.1234A>G was called.", "an HGVS coding change"),
    ("Cross-referenced to MONDO:0008234 for the disease.", "a disease identifier"),
    ("Reported karyotype: 47 in two lines.", "an explicit karyotype"),
    ("Cells were 47,XY in the sample.", "an aneuploid karyotype string"),
    ("Tagged rs123456 in the report.", "a dbSNP identifier"),
])
def test_identifying_content_is_refused(text, description) -> None:
    with pytest.raises(ValueError, match="refusing to publish"):
        publish_guard.check_text(text, "a test page")


def test_ordinary_report_prose_passes() -> None:
    publish_guard.check_text(
        "INVENTEDAZOLE reaches the shortlist on network proximity alone. The endpoint is "
        "secondary prevention. Its proximity z-score is -2.07.", "a test page")


def test_the_endpoint_must_be_qualified_as_secondary() -> None:
    """D4: the bare word overstates the claim, and this is a report about a child."""
    with pytest.raises(ValueError, match="SECONDARY prevention"):
        publish_guard.check_text("The endpoint for this shortlist is chemoprevention.",
                                 "a test page")


def test_the_qualified_endpoint_passes() -> None:
    publish_guard.check_text(
        "The endpoint is secondary chemoprevention, meaning recurrence risk.", "a page")


def test_the_endpoint_is_checked_per_sentence_not_per_document() -> None:
    """A qualifier in the introduction does not travel with a screenshotted caption."""
    with pytest.raises(ValueError, match="SECONDARY prevention"):
        publish_guard.check_text(
            "This report concerns secondary prevention throughout. "
            "Table 3 ranks candidates for chemoprevention.", "a test page")


def test_burden_fields_are_refused_structurally() -> None:
    """A per-chromosome result always arrives as a key with a value under it."""
    payload = {"l0": {"per_contig": {"7": {"verdict": "flagged"}}}}
    with pytest.raises(ValueError, match="per-chromosome"):
        publish_guard.assert_no_burden_detail(payload, "report.json")


def test_burden_fields_are_refused_however_deeply_nested() -> None:
    payload = {"a": [{"b": [{"mosaic_fraction_if_gain": 0.21}]}]}
    with pytest.raises(ValueError, match="per-chromosome"):
        publish_guard.assert_no_burden_detail(payload, "report.json")


def test_naming_a_withheld_field_is_allowed() -> None:
    """The report says which fields it withholds; withholding silently is worse.

    This is a regression test: an earlier guard matched the field NAMES as text and so
    rejected strip_withheld()'s own output.
    """
    published = publish_guard.strip_withheld(
        {"method": {"statistic": "median |BAF - 0.5|"},
         "sensitivity": {"min_detectable_mosaic_fraction": 0.11, "median_depth": 31.4},
         "per_contig": {"1": {"verdict": "flagged", "mosaic_fraction_if_gain": 0.3}}})
    publish_guard.assert_publishable(published, "report.json")
    assert "per_contig" in published["withheld"]


def test_strip_withheld_keeps_the_method_and_drops_the_result() -> None:
    published = publish_guard.strip_withheld({
        "method": {"statistic": "median |BAF - 0.5|", "window_bp": 1000000},
        "sensitivity": {"min_detectable_mosaic_fraction": 0.11, "median_depth": 31.4},
        "per_contig": {"1": {"verdict": "flagged"}},
        "inferred_sex": "XY",
    })
    assert published["method"]["statistic"] == "median |BAF - 0.5|"
    assert published["sensitivity_envelope"]["min_detectable_mosaic_fraction"] == 0.11
    assert "per_contig" not in published
    assert "inferred_sex" not in published
    # median_depth describes the sample, not the assay's resolving power; not published.
    assert "median_depth" not in published["sensitivity_envelope"]


def test_assert_publishable_accepts_a_string_or_a_structure() -> None:
    publish_guard.assert_publishable("plain prose about secondary prevention", "a page")
    publish_guard.assert_publishable({"note": ["nested", {"deeper": "fine"}]}, "a page")


# ------------------------------------------------------------ rendering


def test_every_value_is_escaped() -> None:
    """Rendered text includes openFDA prose quoted straight from the source."""
    html = render.table(["Label text"], [("<script>alert(1)</script> & co",)])
    assert "<script>" not in html
    assert "&lt;script&gt;" in html and "&amp;" in html


def test_a_page_cannot_be_written_without_the_disclaimer() -> None:
    """Structural, not editorial: page() emits it, callers cannot omit it."""
    html = render.page("A title", "a subtitle", "<p>body</p>")
    assert "Hypothesis generation only" in html
    assert "not medical advice" in html


def test_a_page_carries_no_javascript() -> None:
    """D2: an interactive filter presents as a clinical decision aid."""
    html = render.page("A title", "a subtitle", render.table(["A"], [("b",)]))
    lowered = html.lower()
    assert "<script" not in lowered and "onclick" not in lowered


def test_the_rendered_page_passes_its_own_guard() -> None:
    html = render.page("A title", "a subtitle",
                       render.block("<p>Endpoint is secondary prevention.</p>"))
    publish_guard.check_text(html, "a test page")


def test_the_endpoint_statement_satisfies_the_guard() -> None:
    """The constant every page prints must itself be publishable."""
    publish_guard.check_text(render.ENDPOINT_STATEMENT, "the endpoint statement")
    publish_guard.check_text(render.DISCLAIMER, "the disclaimer")


# ------------------------------------------------------------ the rendering rules


def test_absent_contradicting_evidence_never_renders_as_a_blank() -> None:
    """D2 makes the field required because a blank reads as 'nobody looked'."""
    consulted = _contradiction_text({"contradicting_evidence": "", "model_consulted": True})
    assert "Searched, none found" in consulted

    not_consulted = _contradiction_text({"contradicting_evidence": "",
                                         "model_consulted": False})
    assert "Not searched" in not_consulted
    assert not_consulted.strip()


def test_present_contradicting_evidence_is_rendered_verbatim() -> None:
    text = "The only in vitro record used a transformed line."
    assert _contradiction_text({"contradicting_evidence": text}) == text


def test_verdicts_are_read_by_their_real_shape() -> None:
    grouped = _verdicts_by_candidate({"seed": 42, "candidates": [
        {"chembl_id": "CHEMBL1", "verdicts": [{"rule": "pediatric_use",
                                               "verdict": "pass"}]}]})
    assert grouped["CHEMBL1"][0]["rule"] == "pediatric_use"


def test_a_verdicts_artifact_of_the_wrong_shape_is_refused() -> None:
    """An empty safety section reads as 'no findings', not as 'not loaded'."""
    with pytest.raises(ValueError, match="not the expected"):
        _verdicts_by_candidate([{"chembl_id": "CHEMBL1"}])
    with pytest.raises(ValueError, match="not the expected"):
        _verdicts_by_candidate({"rows": []})
