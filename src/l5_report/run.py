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
    report-ready form, the exclusions table, and video stills/assets.

Guardrail
    **Nothing published may re-identify the child or family.** In a ~50-patient
    worldwide population, a plot of per-chromosome aneuploidy burden, an exact variant
    coordinate, or a phenotype narrative can be identifying on its own. Publish
    aggregate or categorical forms; never raw per-sample genomic detail.

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
    #  - Impact: the safety-triaged candidate table with rationale, provenance, and the
    #    companion exclusions table.
    #  - Innovation: the aneuploidy-burden feature (aggregate form only) and the
    #    Claude-in-the-loop contradiction-search trace.
    #  - Scalability: the same pipeline run end-to-end on a second monogenic disease,
    #    with only the config changed -- that is the whole Scalability argument.
    #
    # Read only from `results_dir`; never reach back into `data_dir` for a figure.
    raise NotImplementedError("l5_report.run is a scaffold stub")
