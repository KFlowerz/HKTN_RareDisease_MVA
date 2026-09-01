"""L2 Channel A -- Knowledge-graph link prediction.

Purpose
    Predict drug-target/drug-gene links in a biomedical knowledge graph and rank drugs
    by predicted association with the causal gene. Uses explainable meta-paths so every
    nomination can be traced through named intermediate nodes rather than reported as a
    bare embedding score.

Inputs
    L1 target sets; a KG built from Open Targets / PrimeKG / Hetionet.

Outputs
    Ranked candidates, each with the supporting meta-path(s) -- e.g.
    ``drug -[targets]-> protein -[interacts]-> BUB1B`` -- and the model's score.

Guardrail
    **Anchor the query at the GENE node, never at the disease node.** MVA is absent or
    near-empty in every public KG; a disease-anchored query would return nothing and be
    misread as "no candidates exist". Gene-level anchoring is this channel's entire
    cold-start defense and must not be refactored away.

    Report meta-paths for every candidate. A KG score without a traversable path is not
    admissible evidence in this pipeline.
"""

from __future__ import annotations


def generate(config: dict) -> None:
    """Generate candidates from knowledge-graph link prediction.

    Args:
        config: Parsed pipeline configuration; uses the L1 target sets, ``seed``, and
            ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: build or load the KG (Open Targets + Hetionet/PrimeKG merge, deduplicated on
    # node identity); train or load a link-prediction model (RGCN / CompGCN / TransE) with
    # negative sampling that excludes known positives; score `drug -> target` edges for
    # every node in the L1 upstream and downstream sets; and extract the top meta-paths
    # per prediction (GNNExplainer or a path-ranking pass) as the explanation.
    #
    # Split edges by time or by source to avoid leakage, and record the KG snapshot
    # version -- a judge must be able to rebuild the identical graph.
    raise NotImplementedError("channel_a_kg.generate is a scaffold stub")
