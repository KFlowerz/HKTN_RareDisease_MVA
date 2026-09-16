"""Tests for the L1 disease module: graph expansion and the upstream/downstream split.

Every fixture is an invented graph with invented protein and pathway ids. No real gene,
protein, pathway or disease appears -- which is also the point of
``test_l1_names_no_gene_or_disease``: L1 carries the Scalability claim, and that claim is
only true if the code contains nothing specific to this disease.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.l1_target import module as mod

REPO_ROOT = Path(__file__).resolve().parent.parent

# A small invented interactome: P0 sits in a tight cluster (P1, P2), which connects
# through P3 to a looser neighbourhood (P4, P5), plus an unrelated pair (P8, P9).
EDGES = [
    ("P0", "P1", 950), ("P0", "P2", 920), ("P1", "P2", 900),
    ("P0", "P3", 800), ("P3", "P4", 750), ("P4", "P5", 720),
    ("P5", "P6", 400),                      # below the default cutoff
    ("P8", "P9", 990),                      # separate component
    ("P0", "P0", 999),                      # self-loop
]


def _graph(score_min: int = mod.DEFAULT_SCORE_MIN):
    return mod.build_graph(EDGES, score_min=score_min)


# ---------------------------------------------------------------------------- graph


def test_build_graph_applies_the_cutoff_and_drops_self_loops() -> None:
    graph = _graph()
    assert ("P5", "P6") not in graph.edges          # 400 < 700
    assert not any(a == b for a, b in graph.edges)
    assert graph["P0"]["P1"]["weight"] == pytest.approx(0.95)
    assert graph["P0"]["P1"]["score"] == 950


def test_a_lower_cutoff_admits_weaker_evidence() -> None:
    assert ("P5", "P6") in mod.build_graph(EDGES, score_min=300).edges


# ----------------------------------------------------------------------------- seeds


def test_seed_set_is_the_protein_plus_its_strongest_partners() -> None:
    assert mod.seed_set(_graph(), "P0", seed_score_min=900) == ["P0", "P1", "P2"]


def test_seed_cutoff_controls_the_complex_proxy() -> None:
    assert mod.seed_set(_graph(), "P0", seed_score_min=999) == ["P0"]


def test_absent_protein_raises_rather_than_seeding_noise() -> None:
    with pytest.raises(KeyError, match="not in the interaction graph"):
        mod.seed_set(_graph(), "P_NOT_PRESENT")


# ------------------------------------------------------------------------------ walk


def test_walk_scores_fall_off_with_distance_from_the_seeds() -> None:
    graph = _graph()
    scores = mod.random_walk(graph, ["P0"])
    assert scores["P1"] > scores["P4"] > scores["P5"]
    assert scores["P8"] == pytest.approx(0.0, abs=1e-9)   # unreachable component


def test_walk_is_deterministic() -> None:
    graph = _graph()
    assert mod.random_walk(graph, ["P0"]) == mod.random_walk(graph, ["P0"])


def test_restart_probability_changes_how_far_the_module_reaches() -> None:
    graph = _graph()
    tight = mod.random_walk(graph, ["P0"], restart=0.9)
    loose = mod.random_walk(graph, ["P0"], restart=0.3)
    assert loose["P5"] > tight["P5"]


def test_walk_without_a_present_seed_raises() -> None:
    with pytest.raises(ValueError, match="no seed"):
        mod.random_walk(_graph(), ["P_NOT_PRESENT"])


# --------------------------------------------------------------------------- module


def test_module_keeps_seeds_and_respects_its_size() -> None:
    graph = _graph()
    scores = mod.random_walk(graph, ["P0"])
    chosen = mod.select_module(scores, ["P0", "P1"], size=4)
    assert chosen[:2] == ["P0", "P1"]
    assert len(chosen) == 4


def test_module_never_drops_a_seed_even_when_size_is_smaller() -> None:
    scores = mod.random_walk(_graph(), ["P0"])
    assert mod.select_module(scores, ["P0", "P1", "P2"], size=1) == ["P0", "P1", "P2"]


# ---------------------------------------------------------------------------- split

PATHWAYS = {
    "P0": {"PW-SPECIFIC", "PW-HUGE"},
    "P1": {"PW-SPECIFIC"},                  # shares a specific pathway with P0
    "P4": {"PW-OTHER", "PW-HUGE"},          # shares only the oversized one
    "P5": set(),                            # no pathway at all
}
# PW-HUGE has more members than the bound allows, so it must not make anything upstream.
PATHWAYS.update({f"FILLER{i}": {"PW-HUGE"} for i in range(10)})


def test_sharing_a_specific_pathway_puts_a_node_upstream_with_its_evidence() -> None:
    upstream, downstream, evidence = mod.split_targets(
        ["P0", "P1", "P4", "P5"], PATHWAYS, "P0", max_pathway_size=5)
    assert upstream == ["P0", "P1"]
    assert downstream == ["P4", "P5"]
    assert evidence["P1"] == ["PW-SPECIFIC"]


def test_an_oversized_pathway_is_too_generic_to_define_upstream() -> None:
    """Without the size bound, PW-HUGE would drag P4 upstream and collapse the split."""
    upstream, _, _ = mod.split_targets(["P0", "P4"], PATHWAYS, "P0", max_pathway_size=5)
    assert "P4" not in upstream

    permissive, _, _ = mod.split_targets(["P0", "P4"], PATHWAYS, "P0", max_pathway_size=100)
    assert "P4" in permissive


def test_every_module_node_lands_in_exactly_one_set() -> None:
    module = ["P0", "P1", "P4", "P5"]
    upstream, downstream, _ = mod.split_targets(module, PATHWAYS, "P0", max_pathway_size=5)
    assert sorted(upstream + downstream) == sorted(module)
    assert not set(upstream) & set(downstream)


def test_a_causal_protein_without_a_specific_pathway_is_not_forced_upstream() -> None:
    """The rule is derived, so it reports an empty upstream set rather than inventing one."""
    upstream, downstream, _ = mod.split_targets(["P0", "P4"], PATHWAYS, "P0", max_pathway_size=1)
    assert upstream == []
    assert downstream == ["P0", "P4"]


# ----------------------------------------------------------------------- scalability


def test_l1_names_no_gene_or_disease() -> None:
    """L1 must run unchanged for another monogenic disease -- so it names none.

    Guards the Scalability claim directly: a panel symbol or disease name appearing in L1
    source would mean the module is tuned to this case, whatever the docs say.
    """
    # Gene symbols and the biology of this particular defect. The project's own name is
    # not on the list: it appears legitimately in an env-var name and a user-agent string,
    # and it says nothing about which disease the method is tuned for.
    forbidden = ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192",
                 "aneuploidy", "mitotic", "spindle", "kinetochore", "checkpoint")
    for path in sorted((REPO_ROOT / "src" / "l1_target").glob("*.py")):
        text = path.read_text(encoding="utf-8").lower()
        for term in forbidden:
            assert term.lower() not in text, (
                f"{path.name} names {term!r}; L1 must stay disease-agnostic"
            )
