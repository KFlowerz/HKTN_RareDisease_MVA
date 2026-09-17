"""Tests for the Channel B proximity measure.

Every input is an invented graph with invented node names. No gene symbol, drug name,
pathway or disease term appears here -- the same rule L1's tests follow, and for the same
reason: this channel carries the Scalability claim, which is only true if the method knows
nothing about this case.
"""

from __future__ import annotations

import random
from pathlib import Path

import networkx as nx
import pytest

from src.l2_channels import proximity

REPO_ROOT = Path(__file__).resolve().parent.parent


def _path_graph(n: int) -> nx.Graph:
    """``n0 - n1 - ... - n(n-1)``, so distances are known by construction."""
    graph = nx.Graph()
    graph.add_edges_from((f"n{i}", f"n{i + 1}") for i in range(n - 1))
    return graph


# ------------------------------------------------------------------ distances


def test_module_members_are_at_distance_zero() -> None:
    distance = proximity.distances_from_module(_path_graph(5), ["n0"])
    assert distance["n0"] == 0


def test_distance_is_hops_to_the_nearest_module_node() -> None:
    distance = proximity.distances_from_module(_path_graph(6), ["n0", "n5"])
    assert [distance[f"n{i}"] for i in range(6)] == [0, 1, 2, 2, 1, 0]


def test_unreachable_nodes_are_absent_not_infinite() -> None:
    """Absent forces the caller to decide; infinity silently poisons an average."""
    graph = _path_graph(3)
    graph.add_edge("x0", "x1")
    distance = proximity.distances_from_module(graph, ["n0"])
    assert "x0" not in distance and "x1" not in distance


def test_distances_from_a_module_outside_the_graph_are_empty() -> None:
    assert proximity.distances_from_module(_path_graph(3), ["absent"]) == {}


def test_weighted_edges_do_not_shorten_a_path() -> None:
    """Hop counts, not confidence-weighted distance.

    The confidence cutoff already decided which edges exist; letting a strong edge also be
    a *short* one would count the same evidence twice.
    """
    graph = _path_graph(4)
    for _, _, data in graph.edges(data=True):
        data["weight"] = 0.001
    assert proximity.distances_from_module(graph, ["n0"])["n3"] == 3


# ------------------------------------------------------------- the components


def test_largest_component_is_kept() -> None:
    graph = _path_graph(5)
    graph.add_edge("x0", "x1")
    assert set(proximity.largest_component(graph)) == {f"n{i}" for i in range(5)}


def test_largest_component_of_an_empty_graph_is_empty() -> None:
    assert proximity.largest_component(nx.Graph()).number_of_nodes() == 0


# --------------------------------------------------------------- degree bins


def test_every_node_lands_in_a_pool() -> None:
    graph = nx.barabasi_albert_graph(60, 2, seed=42)
    pools = proximity.degree_bins(graph, min_bin_size=10)
    assert set(pools) == set(graph.nodes())


def test_a_pool_contains_the_node_itself() -> None:
    """Excluding the node would bias the null away from the observation."""
    graph = nx.barabasi_albert_graph(40, 2, seed=42)
    pools = proximity.degree_bins(graph, min_bin_size=10)
    assert all(node in pools[node] for node in graph)


def test_pools_group_comparable_degrees() -> None:
    graph = nx.barabasi_albert_graph(120, 3, seed=42)
    pools = proximity.degree_bins(graph, min_bin_size=20)
    for node, pool in pools.items():
        degrees = [graph.degree(other) for other in pool]
        assert min(degrees) <= graph.degree(node) <= max(degrees)


def test_no_pool_is_a_single_node(  ) -> None:
    """A hub alone in its own bin would resample as itself: a null with no variance."""
    graph = nx.barabasi_albert_graph(200, 2, seed=42)
    pools = proximity.degree_bins(graph, min_bin_size=100)
    assert all(len(pool) >= 100 for pool in pools.values())


def test_a_graph_smaller_than_one_bin_is_one_pool() -> None:
    graph = _path_graph(5)
    pools = proximity.degree_bins(graph, min_bin_size=100)
    assert all(len(pool) == 5 for pool in pools.values())


# ----------------------------------------------------------- closest measure


def test_closest_measure_averages_over_targets() -> None:
    distance = proximity.distances_from_module(_path_graph(6), ["n0"])
    d_c, usable = proximity.closest_measure(distance, ["n1", "n3"])
    assert d_c == 2.0 and usable == ["n1", "n3"]


def test_closest_measure_ignores_unreachable_targets() -> None:
    distance = proximity.distances_from_module(_path_graph(4), ["n0"])
    d_c, usable = proximity.closest_measure(distance, ["n1", "absent"])
    assert d_c == 1.0 and usable == ["n1"]


def test_closest_measure_of_nothing_usable_is_none() -> None:
    distance = proximity.distances_from_module(_path_graph(4), ["n0"])
    assert proximity.closest_measure(distance, ["absent"]) == (None, [])


# ------------------------------------------------------------------ scoring


def _scored(graph, module_nodes, targets, *, permutations=200, seed=42, min_bin_size=10):
    distance = proximity.distances_from_module(graph, module_nodes)
    pools = proximity.degree_bins(graph, min_bin_size=min_bin_size)
    return proximity.proximity(distance, pools, targets, rng=random.Random(seed),
                               permutations=permutations)


def test_targets_inside_the_module_score_zero_distance() -> None:
    graph = nx.barabasi_albert_graph(120, 3, seed=1)
    result = _scored(graph, [0, 1, 2], [0, 1])
    assert result.d_c == 0.0 and result.n_in_module == 2


def test_a_close_target_set_scores_more_negative_than_a_far_one() -> None:
    """The direction that nominates a drug is a negative z."""
    graph = nx.barabasi_albert_graph(300, 3, seed=7)
    module_nodes = list(nx.single_source_shortest_path_length(graph, 0, cutoff=1))
    distance = proximity.distances_from_module(graph, module_nodes)
    near = [n for n, d in distance.items() if d == 1][:5]
    far = sorted((n for n, d in distance.items() if d == max(distance.values())))[:5]

    close_result = _scored(graph, module_nodes, near)
    far_result = _scored(graph, module_nodes, far)
    assert close_result.d_c < far_result.d_c
    assert close_result.z < far_result.z


def test_an_empty_target_set_is_reported_not_scored() -> None:
    result = _scored(nx.barabasi_albert_graph(60, 2, seed=1), [0], [])
    assert result.d_c is None and result.z is None and result.n_usable == 0
    assert "no target" in result.note


def test_unusable_targets_are_counted_in_the_note() -> None:
    """A drug scored on one of nine targets must not read like one scored on nine."""
    graph = nx.barabasi_albert_graph(80, 2, seed=1)
    result = _scored(graph, [0], [1, "absent-a", "absent-b"])
    assert result.n_targets == 3 and result.n_usable == 1
    assert "2 of 3" in result.note


def test_duplicate_targets_are_counted_once() -> None:
    graph = nx.barabasi_albert_graph(80, 2, seed=1)
    assert _scored(graph, [0], [5, 5, 5]).n_targets == 1


def test_a_flat_null_yields_no_z_score_rather_than_a_division_by_zero() -> None:
    graph = _path_graph(4)
    distance = proximity.distances_from_module(graph, ["n0", "n1", "n2", "n3"])
    pools = proximity.degree_bins(graph, min_bin_size=100)
    result = proximity.proximity(distance, pools, ["n1"], rng=random.Random(1),
                                 permutations=50)
    assert result.d_c == 0.0 and result.z is None
    assert "divide by zero" in result.note


def test_scoring_is_reproducible_from_the_seed() -> None:
    graph = nx.barabasi_albert_graph(200, 3, seed=3)
    first = _scored(graph, [0, 1], [10, 20, 30], seed=42)
    second = _scored(graph, [0, 1], [10, 20, 30], seed=42)
    assert (first.z, first.null_mean, first.null_sd) == (second.z, second.null_mean,
                                                         second.null_sd)


def test_a_different_seed_moves_the_null_but_not_the_observation() -> None:
    graph = nx.barabasi_albert_graph(200, 3, seed=3)
    first = _scored(graph, [0, 1], [10, 20, 30], seed=1)
    second = _scored(graph, [0, 1], [10, 20, 30], seed=2)
    assert first.d_c == second.d_c
    assert first.null_mean != second.null_mean


def test_the_null_is_degree_matched_not_uniform() -> None:
    """The whole point of the null: hubs are close to everything.

    Scoring a hub against a uniform null would call it proximal on popularity alone. Drawn
    from its own degree pool, a hub's null mean sits near its own distance.
    """
    graph = nx.barabasi_albert_graph(400, 3, seed=11)
    module_nodes = list(nx.single_source_shortest_path_length(graph, 0, cutoff=1))
    hub = max(graph.degree, key=lambda item: item[1])[0]
    result = _scored(graph, module_nodes, [hub], permutations=500, min_bin_size=100)
    assert result.z is not None
    assert abs(result.z) < 3, "a hub should not look proximal merely for being a hub"


# ------------------------------------------------------------------ ranking


def test_rank_puts_the_most_negative_z_first() -> None:
    results = {
        "far": proximity.ProximityResult(d_c=3.0, z=1.5, null_mean=2.0, null_sd=0.5,
                                         n_targets=1, n_usable=1, n_in_module=0,
                                         permutations=10),
        "near": proximity.ProximityResult(d_c=1.0, z=-2.5, null_mean=2.0, null_sd=0.4,
                                          n_targets=1, n_usable=1, n_in_module=0,
                                          permutations=10),
    }
    assert [name for name, _ in proximity.rank(results)] == ["near", "far"]


def test_rank_puts_unscored_sets_last_however_close_they_look() -> None:
    """An unnormalised distance is not comparable with a normalised one."""
    results = {
        "scored": proximity.ProximityResult(d_c=3.0, z=0.5, null_mean=2.0, null_sd=0.5,
                                            n_targets=1, n_usable=1, n_in_module=0,
                                            permutations=10),
        "unscored": proximity.ProximityResult(d_c=0.0, z=None, null_mean=None, null_sd=None,
                                              n_targets=1, n_usable=1, n_in_module=1,
                                              permutations=0),
    }
    assert [name for name, _ in proximity.rank(results)] == ["scored", "unscored"]


def test_rank_is_deterministic_on_ties() -> None:
    same = dict(d_c=1.0, z=-1.0, null_mean=2.0, null_sd=1.0, n_targets=1, n_usable=1,
                n_in_module=0, permutations=10)
    results = {name: proximity.ProximityResult(**same) for name in ("b", "a", "c")}
    assert [name for name, _ in proximity.rank(results)] == ["a", "b", "c"]


# ---------------------------------------------------------------- guardrails


def test_the_method_names_no_gene_drug_or_disease() -> None:
    """Channel B must work for a disease in no knowledge graph -- so it names none.

    Guards the Scalability claim the same way L1's test does: a panel symbol or disease
    term in this file would mean the method is tuned to this case, whatever the docs say.
    """
    forbidden = ("BUB1B", "CEP57", "TRIP13", "BUB3", "CEP192", "aneuploidy",
                 "mosaic variegated", "MVA", "spindle")
    source = (REPO_ROOT / "src" / "l2_channels" / "proximity.py").read_text(encoding="utf-8")
    for term in forbidden:
        assert term.lower() not in source.lower(), f"{term!r} must not appear in the method"


@pytest.mark.parametrize("permutations", [0, 1])
def test_too_few_permutations_gives_no_z_rather_than_a_fake_one(permutations: int) -> None:
    graph = nx.barabasi_albert_graph(80, 2, seed=1)
    result = _scored(graph, [0], [5], permutations=permutations)
    assert result.z is None and "null could not be built" in result.note


# ------------------------------------------- regressions found in code review


def test_the_null_does_not_depend_on_set_iteration_order() -> None:
    """Targets arrive as a frozenset, whose order depends on string hashing.

    Drawing from the seeded RNG in that order made z and the ranking differ between runs
    at the same seed once hash randomisation was active -- which src/pipeline.py reports
    as the normal case. d_c was stable throughout, so the bug was invisible in the
    observation and lived entirely in the null.
    """
    import networkx as nx

    base = nx.barabasi_albert_graph(300, 3, seed=7)
    graph = nx.relabel_nodes(base, {n: f"9606.ENSP{n:08d}" for n in base})
    module = list(nx.single_source_shortest_path_length(graph, "9606.ENSP00000000", cutoff=1))
    distance = proximity.distances_from_module(graph, module)
    pools = proximity.degree_bins(graph, min_bin_size=20)
    targets = [f"9606.ENSP{n:08d}" for n in (11, 27, 43, 58, 61)]

    runs = [proximity.proximity(distance, pools, order, rng=random.Random(42),
                                permutations=200)
            for order in (targets, list(reversed(targets)), sorted(targets, key=len))]
    assert len({r.z for r in runs}) == 1, "the null must not depend on target order"
    assert len({r.null_mean for r in runs}) == 1
