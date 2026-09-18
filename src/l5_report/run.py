"""L5 -- Reporting.

Purpose
    Assemble the submission artifacts: the candidate table, the figures that carry the
    scientific argument, and the assets for the 3-minute video. Organized against the
    judging rubric -- Scientific Rigor 35%, Impact 25%, Innovation 25%, Scalability 15%
    -- so each claim has a figure and each figure has a claim.

Inputs
    L4 output (surviving candidates + exclusions + benchmark), plus intermediate
    artifacts from L0-L3 for the method figures.

Outputs
    Written under ``config["results_dir"]/l5/``: figures, the candidate table in
    report-ready form, the exclusions table, video stills/assets, and the
    **candidate dossier** -- one static page per surviving candidate, plus a companion
    exclusions page.

    The dossier is the human-facing deliverable (decision D2 in ``mngmt/decisions.md``).
    It is static by design: an interactive drug-filtering tool presents as a clinical
    decision aid whatever disclaimer is attached to it, and clinical recommendation is
    out of scope. A dated, citation-bearing document is read as an argument; a filterable
    app is read as advice.

    **This layer renders; it does not compute.** Every field the dossier shows must
    already exist in an L3 or L4 artifact, so the dossier can be regenerated without
    re-running the paid reasoning step. The per-candidate field contract L3 and L4 are
    obliged to emit is in ``mngmt/decisions.md`` -- if a field is missing at render time
    that is an upstream bug, not something to compute here.

Guardrail
    **Nothing published may re-identify the child or family.** In a ~50-patient
    worldwide population, a plot of per-chromosome aneuploidy burden, an exact variant
    coordinate, or a phenotype narrative can be identifying on its own. Publish
    aggregate or categorical forms; never raw per-sample genomic detail. The boundary is
    set by what the family already shares publicly through their own blog posts
    (``COMPLIANCE.md``) -- this work must not widen it.

    **Say "secondary prevention", not "chemoprevention", wherever the endpoint is named.**
    The ranked list addresses recurrence and second-primary risk -- not prevention of a
    first cancer (``mngmt/decisions.md`` D4).
    The unqualified word overstates the claim, and this is a report about a child.

    Every candidate table carries the **hypothesis-generation-only** disclaimer, the
    Channel C proxy-signature caveat, and the ``n_channels_supporting`` count. A figure
    that shows a ranking without its uncertainty is a recommendation in disguise.

    Show the exclusions. The safety triage is a headline result of this work, not an
    appendix -- what the pipeline refused to propose is as informative as what it did.
"""

from __future__ import annotations


def run(config: dict) -> None:
    """Execute layer L5.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir`` and ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: emit, at minimum:
    #  - Rigor: the multi-channel convergence figure (per-channel ranks vs. consensus)
    #    and the blinded benchmark curves with confidence intervals.
    #  - Impact: the candidate dossier -- one static page per surviving candidate
    #    rendering the D2 field contract (rationale, contradicting evidence, calibrated
    #    confidence, per-rule safety verdict with its source field and quotable snippet,
    #    per-claim references, provenance, caveats), plus the companion exclusions page
    #    carrying the same provenance standard. Render only; never compute a missing
    #    field. `contradicting_evidence` that is empty renders as "searched, none found",
    #    never as a blank -- a blank reads as "not looked for".
    #  - Innovation: the aneuploidy-burden feature (aggregate form only) and the
    #    Claude-in-the-loop contradiction-search trace.
    #  - Scalability: the same pipeline run end-to-end on a second monogenic disease,
    #    with only the config changed -- that is the whole Scalability argument.
    #
    # Read only from `results_dir`; never reach back into `data_dir` for a figure.
    raise NotImplementedError("l5_report.run is a scaffold stub")
