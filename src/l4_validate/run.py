"""L4 -- Validation and safety.

Purpose
    Take the reasoned candidate list from L3 and decide what is defensible to publish:
    mechanistic plausibility, druggability, and clinical status, then a hard pediatric
    safety triage, then a blinded internal benchmark of the pipeline itself.

Inputs
    L3 output (ranked, RxCUI-keyed, with rationale and confidence). Drug annotations per
    ``sources.md`` -- license-clean core for anything redistributed, segregated
    enrichment zone for the rest.

Outputs
    Written under ``config["results_dir"]/l4/``:
      - the surviving candidate table with per-candidate safety verdicts and reasons
      - an explicit ``excluded`` table -- what was dropped and on which rule
      - benchmark results (recovery of known aneuploidy/SAC-relevant compounds)

Guardrail
    Exclusions are **hard gates, not weights**. A candidate that fails the genotoxic or
    pediatric rule is dropped, never merely down-ranked -- MVA is cancer-predisposing,
    and a ranked list is read as a recommendation no matter how it is captioned.

    The excluded table ships with the results. A safety filter whose decisions are
    invisible cannot be audited.

    No NC/ShareAlike-derived annotation may reach a redistributed output path from this
    layer; enrichment fields are joined by RxCUI into the non-redistributed zone.
"""

from __future__ import annotations


def run(config: dict) -> None:
    """Execute layer L4.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir``, ``seed``, and
            ``therapeutic_endpoint`` (for BBB / endpoint-specific criteria).

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: annotate each candidate from the CC-BY-clean core (openFDA label, RxClass,
    # MED-RT), then call `safety_triage.triage` for the hard gates, then
    # `benchmark.run_blinded` for the internal benchmark. Write the survivors and the
    # exclusions as two separate tables, each with the rule that produced the verdict.
    raise NotImplementedError("l4_validate.run is a scaffold stub")
