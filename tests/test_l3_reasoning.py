"""Tests for L3's local-model reasoning step.

No model is called: every test here covers the parts that must be right *before* a
request is made, or that interpret what comes back. The network path itself is exercised
by hand against a running server (docs/runtime.md), because a mocked HTTP round trip
would assert that our own stub behaves like our own stub.

Phenotype terms used here are invented (CLAUDE.md constraint 7).
"""

from __future__ import annotations

import importlib
import json

import pytest

reasoning = importlib.import_module("src.l3_integrate.claude_reasoning")


# --------------------------------------------------------------------------------------
# The guardrail: nothing patient-derived may enter a packet.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("value", "description"),
    [
        ("HP:9000001", "an HPO term id"),
        ("MONDO:0009999", "a disease identifier"),
        ("OMIM:257300", "a disease identifier"),
        ("7:140753336", "a genomic coordinate"),
        ("chr7:140753336", "a genomic coordinate"),
        ("c.1234A>G", "an HGVS coding change"),
    ],
)
def test_packet_with_patient_derived_field_is_refused(value: str, description: str) -> None:
    """A packet carrying subject-derived content must raise, not be sent.

    Fail closed. The alternative is that a value nested somewhere nobody looked reaches a
    process that writes logs, and the first anyone knows is when the logs are read.
    """
    packet = {"candidate": "INVENTIB", "note": f"observed {value} in the cohort"}
    with pytest.raises(ValueError, match=description):
        reasoning.assert_publishable(packet)


def test_packet_is_checked_after_serialization_not_field_by_field() -> None:
    """A forbidden value nested deep inside the packet is still caught."""
    packet = {"candidate": "INVENTIB", "literature_prior": {"caution": {"detail": ["HP:9000002"]}}}
    with pytest.raises(ValueError, match="HPO term id"):
        reasoning.assert_publishable(packet)


def test_ordinary_packet_passes() -> None:
    """Drug identifiers, axes and DOIs are public and must not trip the guard."""
    packet = {
        "candidate": "INVENTIB",
        "literature_prior": {"axis": "RAF/MEK/ERK dependence", "target": "MEK1/2"},
        "citation_keys": ["DOI:10.1038/s41467-024-00000-0", "PMID:12345678"],
    }
    reasoning.assert_publishable(packet)  # must not raise


def test_phenotype_channel_contributes_only_its_rank() -> None:
    """Channel D's evidence names diseases resembling the subject; only rank may travel.

    L3 itself reads only ``candidates.tsv`` from that channel. This asserts the reasoning
    step keeps the same boundary rather than reaching for the richer file next door.
    """
    row = {
        "drug_name": "INVENTIB", "rank": "12", "supporting_channels": "phenotype|prior",
        "rank_phenotype": "0.31", "n_channels_supporting": "2", "convergence": "mixed",
    }
    packet = reasoning.build_packet(row, channel_rows={}, citations=[])

    assert packet["phenotype_channel"] == {"rank": "0.31"}
    assert "disease" not in json.dumps(packet).lower()
    reasoning.assert_publishable(packet)


# --------------------------------------------------------------------------------------
# Evidence grade: supplied, never generated (D17 correction 3).
# --------------------------------------------------------------------------------------

def test_evidence_grade_is_absent_from_what_the_model_may_return() -> None:
    """The response schema must not let the model grade evidence.

    During the 2026-09-20 spike a 7B graded cell-line work as "clinical". Moving the
    confidence *number* into Python does not help if a wrong categorical field feeds the
    computation, so the field is supplied and the model is never asked.
    """
    assert "evidence_grade" not in reasoning.RUBRIC_SCHEMA["properties"]
    assert "evidence_grade" not in reasoning.RUBRIC_SCHEMA["required"]


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        ({"n_clinical": "2", "n_in_vivo": "1", "n_in_vitro": "5"}, "clinical"),
        ({"n_clinical": "0", "n_in_vivo": "3", "n_in_vitro": "5"}, "in_vivo"),
        ({"n_clinical": "0", "n_in_vivo": "0", "n_in_vitro": "1"}, "in_vitro"),
        ({"n_clinical": "0", "n_in_vivo": "0", "n_in_vitro": "0"}, "ungraded"),
        ({}, "ungraded"),
    ],
)
def test_strongest_grade(row: dict, expected: str) -> None:
    """The strongest grade present is what travels with the claim."""
    assert reasoning.strongest_grade(row) == expected


def test_only_verified_citations_reach_a_packet() -> None:
    """A citation channel E could not verify is not evidence."""
    references = {"compounds": [{
        "chembl_id": "CHEMBL1",
        "curated_citations": [
            {"citation": "DOI:10.1/verified", "status": "verified"},
            {"citation": "DOI:10.1/unresolved", "status": "unresolved"},
        ],
    }]}
    assert reasoning.citation_index(references) == {"CHEMBL1": ["DOI:10.1/verified"]}


# --------------------------------------------------------------------------------------
# Citation validation: the model may not introduce a source.
# --------------------------------------------------------------------------------------

def test_invented_citation_keys_are_dropped_and_recorded() -> None:
    """An invented key is dropped, and the drop is reported rather than hidden."""
    kept, dropped = reasoning.validate_citations(
        ["DOI:10.1/real", "DOI:10.1/hallucinated"], ["DOI:10.1/real"])
    assert kept == ["DOI:10.1/real"]
    assert dropped == ["DOI:10.1/hallucinated"]


def test_bracketed_keys_are_normalised_before_matching() -> None:
    """Models bracket their keys; that must not look like an invented citation."""
    kept, dropped = reasoning.validate_citations(["[DOI:10.1/real]", " DOI:10.1/real "],
                                                 ["DOI:10.1/real"])
    assert kept == ["DOI:10.1/real"]
    assert dropped == []


def test_no_citations_returned_is_not_an_error() -> None:
    """A rationale may legitimately rest on no citation; that is not a failure."""
    assert reasoning.validate_citations(None, ["DOI:10.1/real"]) == ([], [])


# --------------------------------------------------------------------------------------
# Confidence: computed here, bounded, and never rewarded for asking fewer questions.
# --------------------------------------------------------------------------------------

def _rubric(**overrides) -> dict:
    base = {"mechanism_supported": True, "contradiction_found": False}
    base.update(overrides)
    return base


def test_confidence_is_bounded() -> None:
    """Every combination stays inside [0, 1]."""
    for grade in ("clinical", "in_vivo", "in_vitro", "ungraded", "nonsense"):
        for supported in (True, False):
            for contradicted in (True, False):
                score = reasoning.confidence(
                    _rubric(mechanism_supported=supported, contradiction_found=contradicted),
                    grade=grade, n_channels=5, n_adverse=3, discriminating=True)
                assert 0.0 <= score <= 1.0


def test_contradiction_lowers_confidence() -> None:
    """Finding evidence against a candidate must reduce its score, never raise it."""
    kwargs = {"grade": "in_vivo", "n_channels": 2, "n_adverse": 0, "discriminating": True}
    clean = reasoning.confidence(_rubric(), **kwargs)
    contradicted = reasoning.confidence(_rubric(contradiction_found=True), **kwargs)
    assert contradicted < clean


def test_unsupported_mechanism_lowers_confidence_more_than_a_weak_grade() -> None:
    """If the evidence does not support the mechanism, grade cannot rescue it."""
    unsupported_clinical = reasoning.confidence(
        _rubric(mechanism_supported=False), grade="clinical", n_channels=1,
        n_adverse=0, discriminating=False)
    supported_invitro = reasoning.confidence(
        _rubric(), grade="in_vitro", n_channels=1, n_adverse=0, discriminating=False)
    assert unsupported_clinical < supported_invitro


def test_discriminating_convergence_outweighs_broad_convergence() -> None:
    """Agreement among broad channels is agreement that a drug exists."""
    kwargs = {"grade": "in_vitro", "n_channels": 2, "n_adverse": 0}
    assert (reasoning.confidence(_rubric(), discriminating=True, **kwargs)
            > reasoning.confidence(_rubric(), discriminating=False, **kwargs))


def test_adverse_records_lower_confidence() -> None:
    """Channel E's own adverse-direction count must count against the candidate."""
    kwargs = {"grade": "in_vivo", "n_channels": 2, "discriminating": True}
    assert (reasoning.confidence(_rubric(), n_adverse=2, **kwargs)
            < reasoning.confidence(_rubric(), n_adverse=0, **kwargs))


# --------------------------------------------------------------------------------------
# The endpoint guarantee.
# --------------------------------------------------------------------------------------

def test_reasoning_refuses_a_non_loopback_endpoint() -> None:
    """D17: a local model behind a routable socket is a hosted API."""
    with pytest.raises(ValueError, match="loopback"):
        reasoning.resolve_endpoint({"reasoning_endpoint": "http://10.0.0.5:8080"})


def test_missing_l3_table_is_reported_not_guessed(tmp_path) -> None:
    """Without L3's table there is nothing to reason about, and that must say so."""
    with pytest.raises(FileNotFoundError, match="no L3 candidate table"):
        reasoning.reason_over_candidates({"results_dir": str(tmp_path), "seed": 42})


# --------------------------------------------------------------------------------------
# Evidence tiers (D18): the kind of claim, not its strength.
# --------------------------------------------------------------------------------------

def test_a_graded_citation_puts_a_candidate_in_the_literature_tier() -> None:
    packet = {"citation_keys": ["DOI:10.1/real"], "evidence_grade_supplied": "in_vitro"}
    assert reasoning.evidence_tier(packet) == "literature"


@pytest.mark.parametrize(
    "packet",
    [
        {"citation_keys": [], "evidence_grade_supplied": "ungraded"},
        {"citation_keys": [], "evidence_grade_supplied": "in_vitro"},
        # A citation with no grade is not a graded finding: channel E verified the DOI
        # but recorded no record behind it at any grade.
        {"citation_keys": ["DOI:10.1/real"], "evidence_grade_supplied": "ungraded"},
        {},
    ],
)
def test_without_graded_literature_a_candidate_is_network_only(packet: dict) -> None:
    assert reasoning.evidence_tier(packet) == "network_only"


def test_network_only_candidates_state_the_absence_rather_than_asserting_a_judgement() -> None:
    """D18: absence is the finding, and `mechanism_supported` is not a false answer.

    On the first full run 82 of 83 survivors were network-only. Recording `false` for
    them would report "the evidence does not support this" when the truth is "there is no
    evidence" -- a different claim, and the one that matters to a reader.
    """
    assert "No published evidence" in reasoning.NO_LITERATURE
    assert "network proximity" in reasoning.NO_LITERATURE


def test_tiering_is_reflected_in_the_caveats() -> None:
    """The tier distinction is most of the result, so it travels with the output."""
    assert any("two tiers" in c or "KIND of evidence" in c for c in reasoning.CAVEATS)
    assert any("not sent to the model" in c for c in reasoning.CAVEATS)
