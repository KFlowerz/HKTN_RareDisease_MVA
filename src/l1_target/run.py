"""L1 -- Target definition.

Purpose
    Expand the single causal gene from L0 into something druggable: the protein and the
    complex it sits in, then the surrounding interactome neighbourhood (the **disease
    module**). Split the result into two target sets, because the two therapeutic
    strategies point at different proteins:

      - **upstream**  -- restore mitotic fidelity (repair the SAC itself). Scientifically
        the most direct, and largely undruggable; carried so it can be discussed
        honestly rather than quietly dropped.
      - **downstream** -- buffer the consequences of aneuploidy (proteotoxic/replication
        stress, senescence, pre-malignant clone survival). This is where
        chemoprevention-style candidates come from.

Inputs
    L0 output (causal gene, variant effect, affected pathway module).
    Interactome and pathway resources: Reactome, STRING, Open Targets.

Outputs
    Written under ``config["results_dir"]``:
      - the disease module as a node/edge set with per-node provenance and confidence
      - ``targets_upstream`` and ``targets_downstream`` gene/protein sets
      - the seed-expansion parameters used, so the module is reproducible

Guardrail
    Module boundaries must be *derived and recorded*, not curated by hand to fit a
    desired answer -- record the expansion algorithm, its parameters, and the STRING
    confidence cutoff alongside the output. Nothing MVA-specific is hardcoded: this
    layer must work for any monogenic causal gene, which is what makes the pipeline
    generalize (Scalability, 15% of the rubric).
"""

from __future__ import annotations


def run(config: dict) -> None:
    """Execute layer L1.

    Args:
        config: Parsed pipeline configuration. Uses ``causal_gene`` (from L0 output,
            not assumed), ``therapeutic_endpoint``, ``results_dir``, and ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: implement L1 as four steps.
    #  1. Resolve the causal gene to UniProt, then to the complexes containing it
    #     (Reactome / CORUM) -- for BUB1B that is the mitotic checkpoint complex; keep
    #     the resolution generic.
    #  2. Build the disease module: seed the interactome (STRING, confidence >= a
    #     configurable cutoff) with the complex members and expand by a documented
    #     algorithm (DIAMOnD or random-walk-with-restart). Record every parameter.
    #  3. Partition into upstream (SAC machinery, kinetochore, mitotic checkpoint) and
    #     downstream (aneuploidy-stress response: proteostasis, replication stress,
    #     p53/senescence, autophagy) sets. Partition by pathway membership, not by
    #     hand-picking.
    #  4. Attach druggability evidence per node from Open Targets (tractability buckets)
    #     so L2 channels can weight reachable targets -- but do NOT filter here; L4 owns
    #     exclusion decisions.
    # If `therapeutic_endpoint` is None, emit BOTH target sets and let L3 carry them
    # forward; do not silently pick one.
    raise NotImplementedError("l1_target.run is a scaffold stub")
