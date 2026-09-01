"""L2 Channel B -- Network proximity.

Purpose
    Score drugs by the shortest-path proximity of their target sets to the L1 disease
    module in the human interactome (Guney/Barabasi). Drugs whose targets sit inside or
    adjacent to the module are nominated.

Inputs
    L1 disease module (node set); a drug-target bipartite set (targets resolved to the
    same protein identifier space as the interactome); the interactome itself.

Outputs
    Ranked candidates with the raw proximity ``d_c``, a degree-preserving null
    distribution, and the resulting z-score.

Guardrail
    This channel exists precisely because it **does not require the disease to appear in
    any knowledge graph** -- it needs only a target set and a module. That property is
    what makes it trustworthy for a hyper-rare disease, so it must never be
    reimplemented in terms of a disease node.

    Proximity must always be reported against a **degree-preserving** random null.
    Well-studied hub proteins are close to everything; a raw distance without the
    degree-matched z-score is a popularity measurement, not evidence.
"""

from __future__ import annotations


def generate(config: dict) -> None:
    """Generate candidates by network proximity to the disease module.

    Args:
        config: Parsed pipeline configuration; uses the L1 disease module, ``seed``, and
            ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: load the interactome as a `networkx` graph (largest connected component);
    # for each drug, compute the closest-measure proximity d_c between its target set and
    # the module; build a degree-preserving null by sampling target sets matched on the
    # binned degree distribution (>= 1000 permutations, seeded from `config["seed"]`);
    # convert to a z-score and rank ascending. Report both `d_c` and z so the ranking is
    # inspectable.
    raise NotImplementedError("channel_b_proximity.generate is a scaffold stub")
