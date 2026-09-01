"""L2 Channel E -- Aneuploidy-stress literature prior + agentic mining.

Purpose
    Capture what the published literature already knows about tolerating, exploiting, or
    buffering aneuploidy -- biology that no structured resource encodes. Aneuploid cells
    carry characteristic proteotoxic, metabolic, and replication stress; compounds that
    selectively pressure or relieve those states are strong prior candidates.

    Combines a curated prior with agentic literature mining (RareAgent-style): an LLM
    agent plans queries, reads abstracts/full text, and extracts drug-mechanism-evidence
    triples with citations.

Inputs
    PubMed / Europe PMC; a curated seed set of aneuploidy-stress literature; the L0
    causal gene and aneuploidy burden.

Outputs
    Ranked candidates, each with **PMID-level citations**, an extracted evidence
    statement, and an evidence-strength grade (in-vitro vs. in-vivo vs. clinical).

Guardrail
    Every claim carries a citation, and the citation must be **verified to exist** -- an
    LLM-extracted PMID is a hypothesis about the literature until resolved against a real
    record. Discard unresolvable citations rather than reporting them; a fabricated
    reference in a rare-disease report is worse than a missing candidate.

    Send **no patient data** to any external API in this channel: queries are built from
    gene symbols and public disease/HPO terms only, never from the subject's variants,
    phenotype narrative, or any identifier.

    Grade evidence honestly. A single-cell-line in-vitro result and a randomized trial do
    not enter L3 at the same weight.
"""

from __future__ import annotations


def generate(config: dict) -> None:
    """Generate candidates from the aneuploidy-stress literature prior.

    Args:
        config: Parsed pipeline configuration; uses the L0 causal gene, ``seed``, and
            ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: assemble the curated prior (aneuploidy tolerance, CIN, proteotoxic stress,
    # HSP90/proteasome/autophagy modulators, senolytics) as a seeded starting set; then
    # run the agentic pass -- plan queries over Europe PMC, retrieve, extract
    # (drug, mechanism, evidence, PMID) tuples, and iterate on gaps.
    #
    # Resolve every extracted PMID against the Europe PMC API and DROP any that does not
    # resolve. Grade each surviving claim (in-vitro / in-vivo / clinical) and rank by
    # graded evidence strength, not by citation count.
    raise NotImplementedError("channel_e_prior.generate is a scaffold stub")
