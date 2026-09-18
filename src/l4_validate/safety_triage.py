"""L4 -- Pediatric safety triage with hard genotoxic exclusion.

Purpose
    Decide which candidates are defensible to propose for a **child with a
    cancer-predisposition syndrome**. Applies, in order:

      1. **Genotoxic / cancer-risk exclusion (hard).** Anything mutagenic, clastogenic,
         aneugenic, or otherwise associated with increased malignancy risk is dropped.
         MVA already produces chromosome missegregation and cancer predisposition;
         proposing a genotoxin is the single worst failure mode this pipeline has.
      2. **Pediatric use (hard).** Documented paediatric use. No adult-only agent passes
         silently, and a label that says safety was *not* established is an exclusion.
      3. **Currently marketed (hard).** A molecule whose every product is discontinued or
         never got past tentative approval is not something to propose for a child.
      4. **Endpoint-specific criteria (recorded).** CNS/BBB penetration where the endpoint
         calls for it -- it does not here (D4).
      5. **Clinical status (recorded).** Boxed warnings and contraindications, quoted.

Inputs
    L3's harmonised candidates with their RxCUIs and UNIIs; openFDA label prose and
    Drugs@FDA product status (:mod:`.labels`); the openFDA pharmacologic classes L3
    attached.

Outputs
    A :class:`Verdict` per candidate per rule -- ``pass`` / ``excluded`` /
    ``not_applicable`` -- each with the rule that fired, the source field it was drawn
    from, and a quotable snippet of the label text.

Guardrail
    Rules 1 to 3 are **exclusions, not penalties**. They are never converted into score
    adjustments, and no channel-consensus score overrides them. A ranked list is read as a
    recommendation however it is captioned.

    **Fail closed.** A candidate whose genotoxicity or paediatric status cannot be
    determined is excluded with reason ``insufficient_evidence``, not passed through.
    Absence of a warning in a label is not evidence of safety; it is usually absence of a
    label. On the real data this excludes the majority of candidates, which is the honest
    outcome of asking a safety question about drugs openFDA does not describe.

    **The decision is lexical, and that is a real limit.** These rules match fixed patterns
    against label prose. They cannot read a sentence the way a pharmacologist does, and
    Section 13.1 in particular is written in careful negatives -- "was not carcinogenic in
    rats" contains the word carcinogenic. Negation is handled explicitly and conservatively:
    positive evidence of harm beats a negative finding in the same section, because the
    cost of wrongly excluding a candidate here is a missed hypothesis, and the cost of
    wrongly including one is proposing a genotoxin for a child who already has a
    chromosome-instability syndrome.

    Currency gotchas encoded here:
      - FDA pregnancy categories A/B/C/D/X are **retired** (PLLR). Do not look for a
        letter grade; parse the narrative Section 8 (Use in Specific Populations).
      - The NLM RxNav drug-interaction API was **discontinued 2024-01-02**. There is no
        structured, CC-BY-clean DDI-with-severity source; either license one or run NLP
        over SPL text, and label DDI output with its provenance either way.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

from . import labels as labels_mod

LOGGER = logging.getLogger(__name__)

#: Ordered so that hard exclusions are evaluated before anything scored.
TRIAGE_RULES = (
    "genotoxic_or_cancer_risk",  # hard exclusion -- non-negotiable for MVA
    "pediatric_use",             # hard exclusion
    "currently_marketed",        # hard exclusion
    "bbb_penetration",           # endpoint-dependent -- see below; not applicable here
    "clinical_status",           # recorded
)
HARD_RULES = TRIAGE_RULES[:3]

#: `bbb_penetration` is scored only when the endpoint calls for CNS exposure. Under
#: `chemoprevention` (gate G2) it does **not**: The documented phenotype does not establish a CNS requirement so the scaffold's
#: assumption of CNS involvement does not hold for this subject. Record the rule as
#: `not_applicable` with that reason rather than scoring it, and never let a
#: BBB-penetrant agent outrank a non-penetrant one on a criterion that does not apply.
ENDPOINTS_REQUIRING_CNS_EXPOSURE = frozenset()

#: Pharmacologic classes that are disqualifying on their face. openFDA's Established
#: Pharmacologic Class and mechanism-of-action fields are curated FDA vocabulary, so this
#: is a structured signal rather than a text match -- the strongest of the three.
CYTOTOXIC_CLASS_PATTERNS = tuple(re.compile(p, re.I) for p in (
    r"\balkylating\b", r"\bantimetabolite\b", r"\btopoisomerase inhibitor\b",
    r"\bnucleoside metabolic inhibitor\b", r"\banthracycline\b",
    r"\bplatinum-based drug\b", r"\bmicrotubule inhibitor\b",
    r"\bdna (?:methyltransferase|synthesis) inhibitor\b", r"\bradiopharmaceutical\b",
    r"\bproteasome inhibitor\b",
))

#: Positive evidence of genotoxicity or carcinogenicity in Section 13.1. Written to match
#: the affirmative constructions labels use; the negation patterns below take precedence
#: only where no affirmative statement is present.
GENOTOXIC_PATTERNS = tuple((name, re.compile(p, re.I)) for name, p in (
    ("carcinogenic_finding", r"\b(?:was|were|is|are)\s+carcinogenic\b"),
    ("increased_tumours", r"\bincreas\w+\s+(?:the\s+)?(?:incidence|number)\s+of\s+"
                          r"(?:\w+\s+){0,3}(?:tumou?rs?|neoplasms?|carcinomas?|adenomas?|"
                          r"lymphomas?|leukemias?)"),
    ("mutagenic_finding", r"\b(?:was|were|is|are)\s+mutagenic\b|\bmutagenic\s+in\b"),
    ("clastogenic_finding", r"\bclastogenic\b"),
    ("aneugenic_finding", r"\baneugen\w+"),
    ("positive_genotoxicity_assay", r"\bpositive\b(?:\s+\w+){0,4}\s+(?:in\s+the\s+)?"
                                    r"(?:ames|bacterial reverse mutation|micronucleus|"
                                    r"chromosom\w+ aberration|mouse lymphoma)"),
    ("genotoxic_statement", r"\b(?:is|was|are|were)\s+genotoxic\b"),
    ("malignancy_risk", r"\bincreased risk of (?:\w+\s+){0,3}(?:malignanc|lymphoma|"
                        r"cancer|secondary (?:malignanc|cancer))"),
    ("secondary_malignancy", r"\bsecondary\s+(?:malignanc\w+|leukemia)\b"),
))

#: Statements that a genotoxicity study was negative. Only consulted when no affirmative
#: pattern matched -- a label reporting both a negative Ames test and a positive
#: micronucleus assay is a positive finding.
NEGATIVE_GENOTOXICITY_PATTERNS = tuple(re.compile(p, re.I) for p in (
    r"\b(?:was|were|is|are)\s+not\s+(?:carcinogenic|mutagenic|genotoxic|clastogenic)\b",
    r"\bno\s+evidence\s+of\s+(?:carcinogenic|mutagenic|genotox|clastogenic)\w*",
    r"\bnegative\b(?:\s+\w+){0,4}\s+(?:in\s+the\s+)?(?:ames|bacterial reverse mutation|"
    r"micronucleus|chromosom\w+ aberration|mouse lymphoma)",
    r"\bnot\s+(?:been\s+)?shown\s+to\s+be\s+(?:carcinogenic|mutagenic|genotoxic)\b",
))

#: Section 8.4 statements that paediatric use is *not* established. Checked before the
#: affirmative patterns: a label commonly establishes use in one age band and denies it in
#: another, and the denial is the part that binds for a child outside the approved band.
PEDIATRIC_NOT_ESTABLISHED = tuple((name, re.compile(p, re.I)) for name, p in (
    # The gap between the subject and the verb is real label text: "Safety and
    # effectiveness in pediatric patients below 12 years of age have not been
    # established". A fixed phrase misses every label that names an age band, which is
    # most of the ones that matter. Bounded and punctuation-free so it cannot reach across
    # a sentence and pair a subject with someone else's verb.
    ("safety_not_established", r"safety\s+and\s+(?:effectiveness|efficacy)[^.;]{0,90}?"
                               r"(?:have|has)\s+not\s+been\s+established"),
    ("not_established", r"\b(?:have|has)\s+not\s+been\s+established[^.;]{0,40}?"
                        r"\bpediatric\b"),
    ("not_indicated", r"\bnot\s+(?:indicated|recommended)\s+(?:for\s+use\s+)?in\s+"
                      r"(?:pediatric|children)"),
    ("contraindicated_in_children", r"\bcontraindicated\b(?:\s+\w+){0,6}\s+"
                                    r"(?:pediatric|children)"),
))

#: Section 8.4 statements that paediatric use *is* established.
PEDIATRIC_ESTABLISHED = tuple((name, re.compile(p, re.I)) for name, p in (
    ("established", r"safety\s+and\s+(?:effectiveness|efficacy)\s+"
                    r"(?:of\s+\w+\s+)?(?:have|has)\s+been\s+established\s+in\s+"
                    r"pediatric"),
    ("indicated_in_children", r"\bindicated\b(?:\s+\w+){0,8}\s+pediatric\s+patients\b"),
    ("age_band", r"\bpediatric\s+patients\s+(?:aged\s+)?\d+\s+(?:years?|months?)\s+"
                 r"(?:of\s+age\s+)?(?:and\s+older|or\s+older|to\s+\d+)"),
))


@dataclass(frozen=True)
class Verdict:
    """One rule's decision about one candidate."""

    rule: str
    verdict: str                 # pass | excluded | not_applicable
    reason: str = ""
    source_field: str = ""       # the openFDA field the decision was drawn from
    source_section: str = ""     # the SPL section a reader would cite
    snippet: str = ""            # quotable label text supporting the verdict

    @property
    def excluded(self) -> bool:
        return self.verdict == "excluded"

    def as_row(self) -> dict:
        return {"rule": self.rule, "verdict": self.verdict, "reason": self.reason,
                "source_field": self.source_field, "source_section": self.source_section,
                "snippet": self.snippet}


#: Words that negate a finding when they appear before it in the same sentence. Section
#: 13.1 is written almost entirely in these constructions -- "was not mutagenic", "no
#: evidence of carcinogenicity", "negative in the Ames assay" -- so a pattern that matches
#: the bare term reads a clean safety record as a positive finding.
NEGATION_CUES = re.compile(
    r"\b(?:not|no|non|never|neither|nor|without|negative|absence|lack(?:s|ed|ing)?|"
    r"free)\b", re.I)
#: Where a negation's scope ends. Sentence punctuation, plus contrast conjunctions: a label
#: reporting "negative in the Ames assay but was clastogenic in human lymphocytes" carries
#: a real positive finding after the "but", and scanning back past it to the earlier
#: "negative" would cancel exactly the finding that matters.
#:
#: Commas are deliberately NOT boundaries. "was not carcinogenic, mutagenic, or
#: clastogenic" is one negation governing a list, and breaking on the commas would read
#: the last two terms as positive findings.
_SENTENCE_START = re.compile(r"[.;:]\s|\b(?:but|however|whereas|nevertheless|yet)\b",
                             re.I)


def _negated(text: str, start: int) -> bool:
    """Whether a match at ``start`` sits after a negation cue in its own sentence.

    Scoped to the sentence so that a negative result in one sentence cannot cancel a
    positive finding in the next. This is the same problem channel D hit with the
    phenotype document, and it is handled the same way: conservatively, and in the
    direction that keeps a real finding rather than discarding it.

    The asymmetry is deliberate. Reading "was not clastogenic" as a positive finding costs
    a candidate -- it excluded selumetinib, the single most paediatric-ready compound in
    this pipeline, on a sentence that says the opposite. Reading a genuine finding as
    negated costs a wrong inclusion, which is why the affirmative patterns are written
    tightly and every exclusion carries the sentence it was drawn from.
    """
    boundaries = [m.end() for m in _SENTENCE_START.finditer(text, 0, start)]
    sentence_start = boundaries[-1] if boundaries else 0
    return bool(NEGATION_CUES.search(text, sentence_start, start))


def _matched(text: str, patterns) -> tuple:
    """First **un-negated** match from a named pattern set, else ``(None, None)``.

    Every occurrence is examined, not just the first: a section reporting a negative Ames
    test and then a positive micronucleus assay is a positive finding, and stopping at the
    first match would return whichever came first in the text.
    """
    for name, pattern in patterns:
        for found in pattern.finditer(text):
            if not _negated(text, found.start()):
                return name, found
    return None, None


def genotoxic_rule(label, pharm_classes) -> Verdict:
    """Hard exclusion on any positive genotoxicity or cancer-risk signal.

    Three independent signals, any one of which excludes: a disqualifying FDA
    pharmacologic class, an affirmative statement in Section 13.1, and an increased-
    malignancy statement anywhere in the boxed warning or warnings sections. Using more
    than one source matters because Section 13.1 is often absent from a label that a
    pharmacologic class alone would condemn.
    """
    joined_classes = " ".join(pharm_classes or ())
    for pattern in CYTOTOXIC_CLASS_PATTERNS:
        if pattern.search(joined_classes):
            return Verdict("genotoxic_or_cancer_risk", "excluded", "cytotoxic_class",
                           "openfda.pharm_class_epc / pharm_class_moa",
                           "FDA Established Pharmacologic Class",
                           labels_mod.snippet(joined_classes, pattern))

    if label is None:
        return Verdict("genotoxic_or_cancer_risk", "excluded", "insufficient_evidence",
                       "", "", "No openFDA label found for this candidate. Absence of a "
                               "label is not evidence of safety, so this rule fails "
                               "closed.")

    for field_name in ("boxed_warning", "warnings_and_cautions", "warnings"):
        text = label.text(field_name)
        name, found = _matched(text, [p for p in GENOTOXIC_PATTERNS
                                        if p[0] in ("malignancy_risk",
                                                    "secondary_malignancy")])
        if name:
            return Verdict("genotoxic_or_cancer_risk", "excluded", name, field_name,
                           labels_mod.SECTIONS.get(field_name, ""),
                           labels_mod.snippet(text, found))

    # 13.1 first, then its parent section: some labels populate only the outer heading.
    sections = ("carcinogenesis_and_mutagenesis_and_impairment_of_fertility",
                "nonclinical_toxicology")
    section = next((s for s in sections if label.text(s)), sections[0])
    text = label.text(section)
    if not text:
        return Verdict("genotoxic_or_cancer_risk", "excluded", "insufficient_evidence",
                       section, labels_mod.SECTIONS[section],
                       "The label carries no carcinogenesis and mutagenesis section, so "
                       "genotoxicity could not be determined. This rule fails closed.")

    name, found = _matched(text, GENOTOXIC_PATTERNS)
    if name:
        return Verdict("genotoxic_or_cancer_risk", "excluded", name, section,
                       labels_mod.SECTIONS[section], labels_mod.snippet(text, found))

    for pattern in NEGATIVE_GENOTOXICITY_PATTERNS:
        if pattern.search(text):
            return Verdict("genotoxic_or_cancer_risk", "pass", "negative_findings_reported",
                           section, labels_mod.SECTIONS[section],
                           labels_mod.snippet(text, found))

    # A section that says neither is not a clearance.
    return Verdict("genotoxic_or_cancer_risk", "excluded", "insufficient_evidence", section,
                   labels_mod.SECTIONS[section],
                   labels_mod.snippet(text) or "The section reports no finding this rule "
                                               "could read either way.")


def pediatric_rule(label) -> Verdict:
    """Hard exclusion unless paediatric use is documented."""
    if label is None:
        return Verdict("pediatric_use", "excluded", "insufficient_evidence", "", "",
                       "No openFDA label found for this candidate, so paediatric use "
                       "could not be established. This rule fails closed.")

    section = "pediatric_use"
    text = label.text(section)
    if not text:
        return Verdict("pediatric_use", "excluded", "insufficient_evidence", section,
                       labels_mod.SECTIONS[section],
                       "The label carries no paediatric-use section. Absence is not "
                       "evidence of safety in children; this rule fails closed.")

    name, found = _matched(text, PEDIATRIC_NOT_ESTABLISHED)
    if name:
        return Verdict("pediatric_use", "excluded", name, section,
                       labels_mod.SECTIONS[section], labels_mod.snippet(text, found))

    name, found = _matched(text, PEDIATRIC_ESTABLISHED)
    if name:
        return Verdict("pediatric_use", "pass", name, section,
                       labels_mod.SECTIONS[section], labels_mod.snippet(text, found))

    return Verdict("pediatric_use", "excluded", "insufficient_evidence", section,
                   labels_mod.SECTIONS[section], labels_mod.snippet(text))


def marketing_rule(marketing) -> Verdict:
    """Hard exclusion when nothing is currently marketed.

    Open Targets records a maximum clinical stage of ``APPROVAL`` for a drug that *was*
    approved, including one later withdrawn. On the real data this rule is what separates
    lorcaserin -- withdrawn in 2020 over an increased incidence of cancer, and carrying
    only a tentative-approval product record -- from a drug still on the market.
    """
    if marketing is None or not marketing.known:
        return Verdict("currently_marketed", "excluded", "insufficient_evidence",
                       "drugsfda.products.marketing_status", "Drugs@FDA",
                       "No Drugs@FDA product record was found, so current marketing "
                       "status could not be established. This rule fails closed.")
    if marketing.marketed:
        return Verdict("currently_marketed", "pass", "marketed",
                       "drugsfda.products.marketing_status", "Drugs@FDA",
                       f"Product marketing status: {', '.join(sorted(marketing.statuses))}.")
    return Verdict("currently_marketed", "excluded", "not_currently_marketed",
                   "drugsfda.products.marketing_status", "Drugs@FDA",
                   f"No product is currently prescription or over-the-counter. Statuses on "
                   f"record: {', '.join(sorted(marketing.statuses))}.")


def bbb_rule(endpoint: str) -> Verdict:
    """Recorded, and not applicable under this project's endpoint (D4)."""
    if endpoint in ENDPOINTS_REQUIRING_CNS_EXPOSURE:
        return Verdict("bbb_penetration", "excluded", "not_evaluated", "", "",
                       "This endpoint requires CNS exposure, and no CNS-penetration "
                       "evidence source is wired in.")
    return Verdict("bbb_penetration", "not_applicable", "endpoint_has_no_cns_requirement",
                   "", "", "The documented phenotype does not establish a CNS requirement, so "
                           "blood-brain-barrier penetration is not a criterion here "
                           "(decision D4). Recorded rather than scored, so a penetrant "
                           "agent never outranks a non-penetrant one on a criterion that "
                           "does not apply.")


def clinical_status_rule(label) -> Verdict:
    """Recorded, never excluding: what a prescriber would be warned about."""
    if label is None:
        return Verdict("clinical_status", "not_applicable", "no_label", "", "",
                       "No openFDA label found.")
    if label.has("boxed_warning"):
        return Verdict("clinical_status", "pass", "boxed_warning_present", "boxed_warning",
                       labels_mod.SECTIONS["boxed_warning"],
                       labels_mod.snippet(label.text("boxed_warning")))
    for field_name in ("contraindications", "warnings_and_cautions", "warnings"):
        if label.has(field_name):
            return Verdict("clinical_status", "pass", "warnings_recorded", field_name,
                           labels_mod.SECTIONS.get(field_name, ""),
                           labels_mod.snippet(label.text(field_name)))
    return Verdict("clinical_status", "pass", "no_warning_sections", "", "",
                   "The label carries no boxed warning, contraindications or warnings "
                   "section. That is a gap in the record, not a safety finding.")


def triage_candidate(*, label, marketing, pharm_classes, endpoint: str) -> list:
    """Every rule's verdict for one candidate, in :data:`TRIAGE_RULES` order."""
    return [
        genotoxic_rule(label, pharm_classes),
        pediatric_rule(label),
        marketing_rule(marketing),
        bbb_rule(endpoint),
        clinical_status_rule(label),
    ]


def survives(verdicts) -> bool:
    """True when no hard rule excluded the candidate."""
    return not any(v.excluded and v.rule in HARD_RULES for v in verdicts)


def first_exclusion(verdicts):
    """The hard rule that dropped a candidate, in rule order, or ``None``."""
    for rule in HARD_RULES:
        for verdict in verdicts:
            if verdict.rule == rule and verdict.excluded:
                return verdict
    return None
