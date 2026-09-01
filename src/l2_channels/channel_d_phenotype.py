"""L2 Channel D -- Phenotype / HPO-driven.

Purpose
    Nominate candidates from the *phenotype* rather than the mechanism: map the subject's
    HPO terms, read from the phenotype document at run time, to drugs indicated for those manifestations or for
    phenotypically-overlapping diseases.

    This is the channel that reaches **symptomatic** candidates -- the ones a purely
    mechanism-first pipeline structurally cannot see.

Inputs
    Phenotype terms from the dataset (HPO-coded); Monarch, Orphanet, Open Targets.

Outputs
    Ranked candidates annotated with the HPO term(s) that produced them and the
    phenotype-similarity measure used.

Guardrail
    Symptomatic relief is **not** disease modification, and the distinction must survive
    into the final report -- a drug nominated here addresses a manifestation, not the
    underlying chromosome missegregation. Tag every candidate with
    ``endpoint="symptomatic"`` so L3 and L5 cannot present it as mechanistic.

    Phenotype terms are patient-derived: use HPO codes and aggregate counts only. Never
    write free-text clinical narrative, ages, dates, or anything else that could
    contribute to re-identification of a child in a ~50-patient population.
"""

from __future__ import annotations


def generate(config: dict) -> None:
    """Generate symptomatic candidates from phenotype terms.

    Args:
        config: Parsed pipeline configuration; uses ``data_dir`` (phenotype file),
            ``seed``, and ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: load and normalize the phenotype to HPO term IDs; compute phenotypic
    # similarity (Resnik/Phenodigm) against diseases in Monarch/Orphanet to find
    # phenotypically adjacent disorders that DO have treatments; pull drugs indicated for
    # those disorders and for the individual HPO manifestations via Open Targets; rank by
    # similarity weighted by the clinical severity of the manifestation addressed.
    #
    # Tag every row with the source HPO term(s) and `endpoint="symptomatic"`.
    raise NotImplementedError("channel_d_phenotype.generate is a scaffold stub")
