"""Decision D5: the transcript policy, its two passes, and the variant-call record.

Invented values only -- no patient variant, coordinate, gene call, or phenotype term
appears here.
"""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pytest
import yaml

from src.l0_genomics import transcripts as tx

REPO_ROOT = Path(__file__).resolve().parent.parent


def _cfg(policy) -> dict:
    return {"annotator": {"transcript_policy": policy}}


def _call(**overrides) -> tx.VariantCall:
    fields = dict(
        gene="GENE1",
        transcript="ENST00000000001.1",
        annotation_pass="primary",
        lof_tier="mane_select",
        mane_release="MANE.GRCh38.v0.0",
        transcript_policy="mane_select_tiered",
        consequence="stop_gained",
        impact="HIGH",
        hgvs_c="c.1A>T",
    )
    fields.update(overrides)
    return tx.VariantCall(**fields)


def test_committed_config_uses_the_decided_policy() -> None:
    """The config carries D5's policy and a MANE release, and the decision log records it.

    Same standard as the therapeutic endpoint: a policy nobody wrote down is
    indistinguishable from one somebody guessed.
    """
    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    annotator = config["annotator"]
    assert annotator["transcript_policy"] == "mane_select_tiered"
    assert annotator["mane_release"]
    assert tx.snpeff_passes(config)

    decisions = (REPO_ROOT / "mngmt" / "decisions.md").read_text(encoding="utf-8")
    assert "## D5" in decisions
    assert annotator["transcript_policy"] in decisions


def test_two_passes_mane_select_first_then_all_transcripts() -> None:
    passes = tx.snpeff_passes(_cfg("mane_select_tiered"))
    assert [name for name, _ in passes] == ["primary", "secondary"]
    assert passes[0][1] == ("-tag", "MANE_Select")
    assert passes[1][1] == ()


@pytest.mark.parametrize("policy", ["canon", "all_transcripts"])
def test_refused_policies_raise_with_their_reason(policy: str) -> None:
    with pytest.raises(ValueError, match="refused"):
        tx.snpeff_passes(_cfg(policy))


def test_unknown_policy_raises() -> None:
    with pytest.raises(ValueError, match="unknown"):
        tx.snpeff_passes(_cfg("longest_transcript"))


@pytest.mark.parametrize("config", [{}, {"annotator": {}}, {"annotator": None}])
def test_missing_policy_has_no_default(config: dict) -> None:
    with pytest.raises(ValueError, match="no default"):
        tx.snpeff_passes(config)


@pytest.mark.parametrize(
    ("primary", "secondary", "tier"),
    [(True, True, "mane_select"), (False, True, "non_mane"), (False, False, "not_lof")],
)
def test_lof_tier(primary: bool, secondary: bool, tier: str) -> None:
    assert tx.lof_tier(primary_lof=primary, secondary_lof=secondary) == tier


def test_primary_lof_without_secondary_lof_is_inconsistent() -> None:
    """The all-transcript pass includes MANE Select, so this can only be a bug."""
    with pytest.raises(ValueError, match="inconsistent"):
        tx.lof_tier(primary_lof=True, secondary_lof=False)


def test_variant_call_carries_the_d5_fields() -> None:
    record = asdict(_call())
    for field in ("transcript", "mane_release", "lof_tier", "annotation_pass", "transcript_policy"):
        assert record[field]


def test_non_mane_tier_only_from_the_secondary_pass() -> None:
    """Passes stay separate: a non-MANE LoF can never be recorded as the primary call."""
    assert _call(lof_tier="non_mane", annotation_pass="secondary").lof_tier == "non_mane"
    with pytest.raises(ValueError, match="secondary pass"):
        _call(lof_tier="non_mane", annotation_pass="primary")


def test_mane_select_tier_only_from_the_primary_pass() -> None:
    with pytest.raises(ValueError, match="primary pass"):
        _call(lof_tier="mane_select", annotation_pass="secondary")


def test_not_lof_may_come_from_either_pass() -> None:
    for annotation_pass in tx.PASSES:
        _call(lof_tier="not_lof", annotation_pass=annotation_pass,
              consequence="missense_variant", impact="MODERATE")


@pytest.mark.parametrize(
    "overrides",
    [
        {"transcript": ""},
        {"mane_release": ""},
        {"annotation_pass": "tertiary"},
        {"lof_tier": "probable"},
        {"impact": "SEVERE"},
        {"transcript_policy": "canon"},
    ],
)
def test_variant_call_rejects_invalid_fields(overrides: dict) -> None:
    with pytest.raises(ValueError):
        _call(**overrides)
