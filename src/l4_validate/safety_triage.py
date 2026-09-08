"""L4 -- Pediatric safety triage with hard genotoxic exclusion.

Purpose
    Decide which candidates are defensible to propose for a **child with a
    cancer-predisposition syndrome**. Applies, in order:

      1. **Genotoxic / cancer-risk exclusion (hard).** Anything mutagenic, clastogenic,
         aneugenic, or otherwise associated with increased malignancy risk is dropped.
         MVA already produces chromosome missegregation and cancer predisposition;
         proposing a genotoxin is the single worst failure mode this pipeline has.
      2. **Pediatric use (hard).** Approved or documented pediatric use, or an
         explicitly justified pediatric rationale. No adult-only agent passes silently.
      3. **Endpoint-specific criteria (soft, recorded).** CNS/BBB penetration where the
         endpoint calls for it; route and formulation feasibility.
      4. **Clinical status (recorded).** Approval status, boxed warnings,
         contraindications, relevant warnings -- extracted from openFDA label prose.

Inputs
    L3 candidates keyed by RxCUI; openFDA label + FAERS + drugsfda; RxClass / MED-RT
    mechanism and pharmacologic class; UNII/NDC.

Outputs
    A verdict per candidate -- ``pass`` / ``excluded`` -- each with the rule that fired,
    the source field it was drawn from, and a quotable snippet of the label text.

Guardrail
    Rules 1 and 2 are **exclusions, not penalties**. Never convert them into score
    adjustments, and never allow a high channel-consensus score to override them.

    Fail closed: a candidate whose genotoxicity or pediatric status cannot be determined
    is excluded with reason ``insufficient_evidence``, not passed through. Absence of a
    warning in a label is not evidence of safety.

    Currency gotchas encoded here:
      - FDA pregnancy categories A/B/C/D/X are **retired** (PLLR). Do not look for a
        letter grade; parse the narrative Section 8 (Use in Specific Populations).
      - The NLM RxNav drug-interaction API was **discontinued 2024-01-02**. There is no
        structured, CC-BY-clean DDI-with-severity source; either license one or run NLP
        over SPL text, and label DDI output with its provenance either way.
"""

from __future__ import annotations

#: Ordered so that hard exclusions are evaluated before anything scored.
TRIAGE_RULES = (
    "genotoxic_or_cancer_risk",  # hard exclusion -- non-negotiable for MVA
    "pediatric_use",             # hard exclusion
    "bbb_penetration",           # endpoint-dependent -- see below; not applicable here
    "clinical_status",           # recorded
)

#: `bbb_penetration` is scored only when the endpoint calls for CNS exposure. Under
#: `chemoprevention` (gate G2) it does **not**: The documented phenotype does not establish a CNS requirement so the scaffold's
#: assumption of CNS involvement does not hold for this subject. Record the rule as
#: `not_applicable` with that reason rather than scoring it, and never let a
#: BBB-penetrant agent outrank a non-penetrant one on a criterion that does not apply.
ENDPOINTS_REQUIRING_CNS_EXPOSURE = frozenset()


def triage(config: dict) -> None:
    """Apply the pediatric safety triage to the candidate list.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir`` and
            ``therapeutic_endpoint``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: build the genotoxic exclusion set from multiple independent signals rather
    # than one flag -- ATC L01 antineoplastics, MED-RT mechanism classes, IARC/NTP
    # carcinogen listings, and label-derived carcinogenesis/mutagenesis text (Section
    # 13.1) -- and treat any positive as excluding.
    #
    # Extract pediatric use and the remaining criteria from openFDA label prose (the
    # Claude layer does the extraction; this module owns the decision rules and keeps
    # the quotable snippet). Emit one row per candidate per rule so every verdict is
    # traceable to a source field.
    raise NotImplementedError("safety_triage.triage is a scaffold stub")
