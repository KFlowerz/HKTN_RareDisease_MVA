"""L4 -- Blinded internal benchmark.

Purpose
    Answer the question a judge will ask: *does this pipeline recover compounds we
    already believe are relevant?* There is **no external ground truth for Track 2** --
    MVA has no approved therapy and no repurposing gold standard exists -- so the only
    honest validation is an internally-held-out set of aneuploidy/SAC-relevant compounds
    that the pipeline was not tuned against.

Inputs
    A held-out positive set (compounds with published aneuploidy-stress, CIN-selective,
    or SAC-modulating activity), plus a degree-matched negative/decoy set. The full
    ranked candidate list from L3.

Outputs
    Recovery metrics under ``config["results_dir"]/l4/benchmark/``: rank of each held-out
    positive, enrichment at k, AUROC/AUPRC against the decoy set, and a per-channel
    breakdown showing which channels carried the signal.

Guardrail
    **Blinded means blinded.** The positive set is constructed and frozen before the
    pipeline is tuned, and no channel, weight, or threshold may be adjusted after seeing
    benchmark results without re-freezing and disclosing that the set was reused. A
    benchmark optimized against is a description of the pipeline, not a test of it.

    Report the negative result if that is what happens. An honest miss is worth more to
    Scientific Rigor (35% of the rubric) than a tuned hit, and a benchmark that always
    passes is evidence of leakage, not quality.

    This measures **pipeline recovery**, not clinical efficacy. Nothing here licenses a
    claim that any candidate works.
"""

from __future__ import annotations


def run_blinded(config: dict) -> None:
    """Run the blinded internal benchmark against the held-out positive set.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir`` and ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: load the frozen positive set (with its freeze date and provenance) and a
    # degree/popularity-matched decoy set -- decoys matched on drug-target degree and
    # literature volume, so the benchmark measures mechanism rather than fame.
    #
    # Compute rank-based recovery (enrichment at k, AUROC, AUPRC) for the combined
    # ranking and for each channel independently, so a single dominant channel cannot
    # hide behind the consensus score. Bootstrap CIs with `config["seed"]`.
    #
    # Record the freeze date and a hash of the positive set alongside the metrics; a
    # benchmark whose contents can drift is not blinded.
    raise NotImplementedError("benchmark.run_blinded is a scaffold stub")
