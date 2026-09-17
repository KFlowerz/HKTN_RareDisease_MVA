"""L2 Channel B -- the network-proximity measure and its degree-preserving null.

Purpose
    Score how close a drug's target set sits to a disease module on the interactome, and
    say whether that closeness is more than the target set's own popularity would predict.
    Pure graph and set operations: no I/O, no network calls, and nothing disease-specific
    or drug-specific.

Inputs
    An interactome as a :class:`networkx.Graph`, a module node set, and a target set per
    drug -- all supplied by :mod:`src.l2_channels.channel_b_proximity`.

Outputs
    Per target set: the closest measure ``d_c``, the mean and standard deviation of a
    degree-matched null, the resulting z-score, and how many targets were usable.

Method
    The **closest measure** of Guney et al. [guney2016] doi:10.1038/ncomms10331: the mean,
    over a drug's targets, of the shortest-path distance from that target to its nearest
    module protein. Because the module is the same for every drug, one multi-source
    breadth-first search from the module yields that distance for every protein in the
    graph at once, and each drug's score is then an average over its own targets.

    Distances are **unweighted hop counts**, as in the source method. The confidence
    cutoff already decided which interactions exist (L1's ``string_score_min``); treating
    a high-confidence edge as also *shorter* would apply the same evidence twice.

Guardrail
    **A raw distance is a popularity measurement.** Well-studied hub proteins are close to
    everything, and drugs are studied against well-studied proteins, so ``d_c`` alone would
    rank by how thoroughly a drug has been characterised. Every score here is reported
    against a null built by resampling each target from other proteins of comparable
    degree [guney2016], and the ranking the channel emits is by z-score.

    **Unreachable targets are reported, never dropped quietly.** A target outside the
    interactome's largest component has no finite distance to the module. Each result
    carries how many targets were usable out of how many were supplied, so a drug scored
    on one of nine targets cannot be mistaken for a drug scored on nine.

    **Nothing here knows what a drug or a disease is.** It takes node sets. That is what
    lets this channel work for a disease that appears in no knowledge graph, and it is
    also what makes the Scalability claim testable.
"""

from __future__ import annotations

import collections
import math
import random
from dataclasses import dataclass, field

import networkx as nx

#: Minimum nodes per degree bin. Degrees are long-tailed, so exact-degree matching would
#: leave a hub alone in its own bin and resample it as itself -- a null identical to the
#: observation. Guney et al. bin degrees to keep each pool sampleable [guney2016].
DEFAULT_MIN_BIN_SIZE = 100
#: The scaffold's floor, and Guney et al.'s order of magnitude. Fewer permutations make
#: the tails of the null -- exactly where an interesting drug sits -- unstable.
DEFAULT_PERMUTATIONS = 1000


@dataclass(frozen=True)
class ProximityResult:
    """One target set's proximity to the module."""

    d_c: float | None          # the closest measure; None when no target was usable
    z: float | None            # (d_c - null mean) / null sd; None when sd is 0 or d_c is
    null_mean: float | None
    null_sd: float | None
    n_targets: int             # targets supplied
    n_usable: int              # targets present in the graph and able to reach the module
    n_in_module: int           # usable targets that are themselves module members
    permutations: int
    note: str = ""             # why a score is missing, when one is
    distances: tuple = field(default_factory=tuple)   # per usable target, for inspection


def largest_component(graph: nx.Graph) -> nx.Graph:
    """The largest connected component, as a subgraph view.

    Distance is undefined between components, and the interactome's small components are
    mostly annotation artefacts. Restricting once, here, keeps every later distance finite
    or explicitly absent rather than silently infinite.
    """
    if graph.number_of_nodes() == 0:
        return graph
    return graph.subgraph(max(nx.connected_components(graph), key=len))


def distances_from_module(graph: nx.Graph, module_nodes) -> dict:
    """``{node: hops to the nearest module node}`` by one multi-source BFS.

    Module members are at distance 0. Nodes the module cannot reach are absent from the
    result rather than mapped to infinity, so a caller must decide what to do about them.
    """
    seeds = [n for n in module_nodes if n in graph]
    if not seeds:
        return {}
    distance = {node: 0 for node in seeds}
    queue = collections.deque(seeds)
    while queue:
        node = queue.popleft()
        for neighbour in graph.neighbors(node):
            if neighbour not in distance:
                distance[neighbour] = distance[node] + 1
                queue.append(neighbour)
    return distance


def degree_bins(graph: nx.Graph, *, min_bin_size: int = DEFAULT_MIN_BIN_SIZE) -> dict:
    """``{node: (node, ...)}`` -- the pool each node may be resampled from.

    Nodes are ordered by degree and cut into consecutive bins of at least ``min_bin_size``;
    a trailing bin too small to stand alone is merged into the previous one. Every node
    maps to the whole pool it belongs to, including itself: a degree-matched draw that
    excluded the node itself would bias the null away from the observation.
    """
    ordered = sorted(graph.nodes(), key=lambda n: (graph.degree(n), n))
    bins, current = [], []
    for node in ordered:
        current.append(node)
        if len(current) >= min_bin_size:
            bins.append(current)
            current = []
    if current:
        if bins:
            bins[-1].extend(current)
        else:
            bins.append(current)

    pools = {}
    for members in bins:
        pool = tuple(members)
        for node in members:
            pools[node] = pool
    return pools


def closest_measure(distance: dict, targets) -> tuple:
    """``(d_c, usable)`` -- the mean distance to the module over reachable targets."""
    usable = [t for t in targets if t in distance]
    if not usable:
        return None, []
    return sum(distance[t] for t in usable) / len(usable), usable


def proximity(distance: dict, pools: dict, targets, *, rng: random.Random,
              permutations: int = DEFAULT_PERMUTATIONS) -> ProximityResult:
    """Score one target set against a degree-matched null.

    Args:
        distance: :func:`distances_from_module` for the module being scored.
        pools: :func:`degree_bins` for the same graph.
        targets: the drug's targets, in the graph's identifier space.
        rng: seeded RNG; the caller owns it, so a run is reproducible.
        permutations: resamples of the target set.

    Returns:
        A :class:`ProximityResult`. ``z`` is negative when the targets sit closer to the
        module than degree-matched proteins do -- the direction that nominates a drug.
    """
    # Sorted, not merely deduplicated: targets arrive as a frozenset of protein ids, whose
    # iteration order depends on string hashing. Drawing from the seeded RNG in that order
    # makes the null -- and therefore the ranking -- differ between runs at the same seed
    # whenever hash randomisation is active, which src/pipeline.py reports as the norm.
    targets = sorted(dict.fromkeys(targets), key=str)
    d_c, usable = closest_measure(distance, targets)
    if d_c is None:
        return ProximityResult(
            d_c=None, z=None, null_mean=None, null_sd=None, n_targets=len(targets),
            n_usable=0, n_in_module=0, permutations=0,
            note="no target is in the interactome's largest component, or none can reach "
                 "the module",
        )

    # Resample only the target set: the module is the same for every drug, so holding it
    # fixed asks "is this drug closer than a drug of comparable target degree?", which is
    # the question the channel is for. Guney et al. randomise both sides, because they
    # compare across many diseases [guney2016].
    sampleable = [t for t in usable if t in pools]
    draws = []
    for _ in range(permutations):
        total, count = 0, 0
        for target in sampleable:
            pick = rng.choice(pools[target])
            if pick in distance:
                total += distance[pick]
                count += 1
        if count:
            draws.append(total / count)

    note = ""
    if len(usable) < len(targets):
        note = (f"{len(targets) - len(usable)} of {len(targets)} target(s) are outside the "
                "interactome's largest component or cannot reach the module, and are not "
                "scored")

    if len(draws) < 2:
        return ProximityResult(
            d_c=d_c, z=None, null_mean=None, null_sd=None, n_targets=len(targets),
            n_usable=len(usable), n_in_module=sum(distance[t] == 0 for t in usable),
            permutations=len(draws),
            note="; ".join(filter(None, [note, "the null could not be built"])),
            distances=tuple(distance[t] for t in usable),
        )

    mean = sum(draws) / len(draws)
    variance = sum((value - mean) ** 2 for value in draws) / (len(draws) - 1)
    sd = math.sqrt(variance)
    if sd == 0:
        z = None
        note = "; ".join(filter(None, [
            note, "every degree-matched resample gave the same distance, so a z-score "
                  "would divide by zero"]))
    else:
        z = (d_c - mean) / sd

    return ProximityResult(
        d_c=d_c, z=z, null_mean=mean, null_sd=sd, n_targets=len(targets),
        n_usable=len(usable), n_in_module=sum(distance[t] == 0 for t in usable),
        permutations=len(draws), note=note,
        distances=tuple(distance[t] for t in usable),
    )


def rank(results: dict) -> list:
    """``[(name, result), ...]`` best first: most negative z, then smallest ``d_c``.

    A target set with no z-score is ranked last whatever its ``d_c``, because an
    unnormalised distance is not comparable with a normalised one.
    """
    def key(item):
        _, result = item
        return (result.z is None, result.z if result.z is not None else 0.0,
                result.d_c if result.d_c is not None else math.inf, item[0])

    return sorted(results.items(), key=key)
