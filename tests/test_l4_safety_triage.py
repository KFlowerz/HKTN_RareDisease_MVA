"""Tests for the L4 safety rules.

Every drug, label sentence and class here is invented. Where a real label's phrasing is
reproduced it is a short factual quotation used to pin a parser behaviour -- the
selumetinib and trametinib sentences are the two that decided a real verdict, and a test
that did not contain them would not protect the behaviour that matters most in this layer.
"""

from __future__ import annotations

import pytest

from src.l4_validate import labels as labels_mod
from src.l4_validate import safety_triage as st


def _label(**sections) -> labels_mod.Label:
    return labels_mod.Label(sections={k: v for k, v in sections.items() if v})


CARC = "carcinogenesis_and_mutagenesis_and_impairment_of_fertility"


# --------------------------------------------------------------------- negation


@pytest.mark.parametrize("text, expected", [
    ("Inventib was not mutagenic or clastogenic in vitro.", None),
    ("Inventib was not mutagenic in a bacterial reverse mutation assay.", None),
    ("No evidence of carcinogenicity was observed.", None),
    ("Inventib was not carcinogenic, mutagenic, or clastogenic.", None),
    ("Inventib was carcinogenic in rats.", "carcinogenic_finding"),
    ("Treatment increased the incidence of hepatic tumors in mice.", "increased_tumours"),
    ("Inventib was negative in the Ames assay but was clastogenic in human lymphocytes.",
     "clastogenic_finding"),
    ("Studies were negative. Inventib was clastogenic at high doses.",
     "clastogenic_finding"),
    ("There is an increased risk of secondary malignancies.", "malignancy_risk"),
])
def test_negation_is_scoped_to_the_clause(text, expected):
    """A safety rule that reads 'was not clastogenic' as a finding excludes good drugs."""
    name, _ = st._matched(text, st.GENOTOXIC_PATTERNS)
    assert name == expected


def test_a_comma_list_stays_under_one_negation():
    """'not carcinogenic, mutagenic, or clastogenic' is one negation over three terms."""
    assert st._matched("Drug was not carcinogenic, mutagenic, or clastogenic.",
                       st.GENOTOXIC_PATTERNS) == (None, None)


def test_a_contrast_conjunction_ends_the_negation_scope():
    """Otherwise an earlier 'negative' cancels the finding that follows 'but'."""
    name, _ = st._matched("Negative in the Ames assay, however it was clastogenic.",
                          st.GENOTOXIC_PATTERNS)
    assert name == "clastogenic_finding"


def test_every_occurrence_is_examined_not_just_the_first():
    text = ("Inventib was not mutagenic in bacteria. Inventib was clastogenic in "
            "human lymphocytes.")
    name, found = st._matched(text, st.GENOTOXIC_PATTERNS)
    assert name == "clastogenic_finding"
    assert "clastogenic" in found.group(0)


# ------------------------------------------------------------- the genotoxic rule


def test_a_cytotoxic_pharmacologic_class_excludes_without_a_label():
    """The FDA's own class vocabulary is a structured signal, not a text match."""
    verdict = st.genotoxic_rule(None, ("Proteasome Inhibitor [EPC]",))
    assert verdict.excluded and verdict.reason == "cytotoxic_class"
    assert "Proteasome Inhibitor" in verdict.snippet


def test_no_label_fails_closed():
    """Absence of a label is not evidence of safety."""
    verdict = st.genotoxic_rule(None, ())
    assert verdict.excluded and verdict.reason == "insufficient_evidence"


def test_an_aneugenic_label_excludes():
    """A real verdict this pipeline made, and the one that matters most for MVA.

    Selumetinib's own label reports an aneugenic mode of action. For a child whose
    syndrome *is* chromosome missegregation, that is the worst possible nomination, and
    the curated prior had it as the most paediatric-ready compound on its axis.
    """
    label = _label(**{CARC: (
        "Selumetinib did result in an increase in micronucleated immature erythrocytes "
        "(chromosome aberrations) in mouse micronucleus studies, predominantly via an "
        "aneugenic mode of action, but at doses > 160 mg/kg.")})
    verdict = st.genotoxic_rule(label, ())
    assert verdict.excluded and verdict.reason == "aneugenic_finding"
    assert "aneugenic" in verdict.snippet


def test_a_negative_genotoxicity_section_passes():
    label = _label(**{CARC: (
        "Carcinogenicity studies with trametinib have not been conducted. Trametinib was "
        "not genotoxic in studies evaluating reverse mutations in bacteria, chromosomal "
        "aberrations in mammalian cells, or micronuclei in rats.")})
    verdict = st.genotoxic_rule(label, ())
    assert verdict.verdict == "pass" and verdict.reason == "negative_findings_reported"


def test_a_section_that_says_neither_fails_closed():
    """'No carcinogenicity studies have been conducted' is not a clearance."""
    label = _label(**{CARC: "No carcinogenicity studies have been conducted."})
    verdict = st.genotoxic_rule(label, ())
    assert verdict.excluded and verdict.reason == "insufficient_evidence"


def test_the_parent_section_is_read_when_the_subsection_is_absent():
    label = _label(nonclinical_toxicology="Inventib was clastogenic in human lymphocytes.")
    verdict = st.genotoxic_rule(label, ())
    assert verdict.excluded and verdict.source_field == "nonclinical_toxicology"


def test_a_malignancy_warning_excludes_even_without_a_carcinogenesis_section():
    label = _label(boxed_warning="There is an increased risk of secondary malignancies.")
    verdict = st.genotoxic_rule(label, ())
    assert verdict.excluded and verdict.source_field == "boxed_warning"


def test_the_snippet_quotes_the_sentence_the_verdict_rested_on():
    """A reader must be able to check the verdict, including against a dose caveat."""
    label = _label(**{CARC: (
        "Inventib was negative in the Ames assay. Inventib was clastogenic at 40 mg/kg.")})
    verdict = st.genotoxic_rule(label, ())
    assert "clastogenic at 40 mg/kg" in verdict.snippet
    assert "Ames" not in verdict.snippet


# ------------------------------------------------------------- the paediatric rule


def test_paediatric_use_not_established_excludes():
    label = _label(pediatric_use="Safety and effectiveness in pediatric patients have "
                                 "not been established.")
    verdict = st.pediatric_rule(label)
    assert verdict.excluded and verdict.reason == "safety_not_established"


def test_an_established_age_band_passes():
    label = _label(pediatric_use="The safety and effectiveness of INVENTIB have been "
                                 "established in pediatric patients 1 year of age and "
                                 "older.")
    assert st.pediatric_rule(label).verdict == "pass"


def test_a_denial_beats_an_establishment_in_the_same_section():
    """A label may establish one age band and deny another; the denial binds."""
    label = _label(pediatric_use=("Safety and effectiveness have been established in "
                                  "pediatric patients 12 years and older. Safety and "
                                  "effectiveness in pediatric patients below 12 years "
                                  "have not been established."))
    assert st.pediatric_rule(label).excluded


def test_a_missing_paediatric_section_fails_closed():
    assert st.pediatric_rule(_label()).excluded
    assert st.pediatric_rule(None).reason == "insufficient_evidence"


# -------------------------------------------------------------- the marketing rule


def test_a_molecule_with_no_marketed_product_is_excluded():
    """Open Targets' APPROVAL means 'was approved'; lorcaserin was withdrawn in 2020."""
    marketing = labels_mod.Marketing(statuses=frozenset({"None (Tentative Approval)"}))
    verdict = st.marketing_rule(marketing)
    assert verdict.excluded and verdict.reason == "not_currently_marketed"


def test_any_marketed_product_passes():
    """Older formulations lapsing is normal for a drug still on the market."""
    marketing = labels_mod.Marketing(statuses=frozenset({"Discontinued", "Prescription"}))
    assert st.marketing_rule(marketing).verdict == "pass"


def test_no_marketing_record_fails_closed():
    assert st.marketing_rule(None).reason == "insufficient_evidence"
    assert st.marketing_rule(labels_mod.Marketing()).excluded


# ------------------------------------------------------------------ the endpoint


def test_bbb_is_recorded_as_not_applicable_under_this_endpoint():
    """D4: the BBB penetration is not a hard criterion under this endpoint."""
    verdict = st.bbb_rule("chemoprevention")
    assert verdict.verdict == "not_applicable"
    assert not verdict.excluded


# ------------------------------------------------------------------ composition


def test_a_hard_rule_excludes_and_a_recorded_rule_does_not():
    verdicts = st.triage_candidate(label=None, marketing=None, pharm_classes=(),
                                   endpoint="chemoprevention")
    assert not st.survives(verdicts)
    assert {v.rule for v in verdicts} == set(st.TRIAGE_RULES)
    assert st.first_exclusion(verdicts).rule == "genotoxic_or_cancer_risk"


def test_the_first_exclusion_follows_rule_order_not_evaluation_order():
    """The reported reason must be stable, so the exclusions table can be counted."""
    label = _label(**{CARC: "Inventib was clastogenic."},
                   pediatric_use="Safety and effectiveness have not been established.")
    verdicts = st.triage_candidate(label=label, marketing=None, pharm_classes=(),
                                   endpoint="chemoprevention")
    assert st.first_exclusion(verdicts).rule == "genotoxic_or_cancer_risk"


def test_a_clean_candidate_survives():
    label = _label(**{CARC: "Inventib was not genotoxic in any assay."},
                   pediatric_use="Safety and effectiveness have been established in "
                                 "pediatric patients 2 years of age and older.")
    marketing = labels_mod.Marketing(statuses=frozenset({"Prescription"}))
    verdicts = st.triage_candidate(label=label, marketing=marketing, pharm_classes=(),
                                   endpoint="chemoprevention")
    assert st.survives(verdicts)
    assert st.first_exclusion(verdicts) is None


def test_every_verdict_names_its_source():
    """D2's contract: rule, verdict, source field and a quotable snippet."""
    label = _label(**{CARC: "Inventib was clastogenic."})
    for verdict in st.triage_candidate(label=label, marketing=None, pharm_classes=(),
                                       endpoint="chemoprevention"):
        row = verdict.as_row()
        assert set(row) == {"rule", "verdict", "reason", "source_field",
                            "source_section", "snippet"}
        assert row["verdict"] in {"pass", "excluded", "not_applicable"}
        assert row["snippet"], f"{verdict.rule} produced no quotable basis"
