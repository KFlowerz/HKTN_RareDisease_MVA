"""Tests for L3 rank aggregation.

Every drug key here is invented (``INV1``, ``INV2``, …). The method is checked against
values computed from its own definition rather than against a remembered number, so a test
cannot quietly encode the same mistake as the implementation.
"""

from __future__ import annotations

import math

import pytest

from src.l3_integrate import aggregate as agg


# ------------------------------------------------------------------- the beta CDF


@pytest.mark.parametrize("x, a, b", [
    (0.25, 1, 1), (0.5, 1, 3), (0.1, 2, 2), (0.75, 3, 1), (0.4, 2, 4), (0.9, 5, 2),
])
def test_beta_cdf_matches_the_binomial_sum_it_claims_to_be(x, a, b):
    """For integer shapes I_x(a,b) is a binomial tail; check against it independently."""
    n = a + b - 1
    expected = sum(math.comb(n, j) * x ** j * (1 - x) ** (n - j) for j in range(a, n + 1))
    assert agg._beta_cdf(x, a, b) == pytest.approx(expected)


def test_beta_cdf_of_a_uniform_is_the_identity():
    """Beta(1,1) is Uniform(0,1), so its CDF must be x."""
    for x in (0.0, 0.1, 0.5, 0.9, 1.0):
        assert agg._beta_cdf(x, 1, 1) == pytest.approx(x)


@pytest.mark.parametrize("x", [-1.0, 0.0, 1.0, 2.0])
def test_beta_cdf_is_bounded(x):
    value = agg._beta_cdf(x, 2, 3)
    assert 0.0 <= value <= 1.0


# ---------------------------------------------------------------------- rho score


def test_a_top_rank_in_every_channel_scores_far_better_than_a_middling_one():
    best, _ = agg.rho_score([0.001, 0.001, 0.001], 3)
    middling, _ = agg.rho_score([0.5, 0.5, 0.5], 3)
    assert best < middling < 1.0


def test_one_excellent_rank_can_beat_several_mediocre_ones():
    """RRA takes the minimum over prefixes, which is the point of the method."""
    single, _ = agg.rho_score([0.0005], 3)
    several, _ = agg.rho_score([0.4, 0.45, 0.5], 3)
    assert single < several


def test_a_candidate_no_channel_ranked_scores_worst():
    assert agg.rho_score([], 3) == (1.0, 1.0)


def test_the_score_is_bonferroni_corrected_by_the_channel_count():
    corrected, raw = agg.rho_score([0.01], 4)
    assert corrected == pytest.approx(min(1.0, raw * 4))


def test_the_score_never_exceeds_one():
    corrected, _ = agg.rho_score([0.9, 0.95], 5)
    assert corrected <= 1.0


# ----------------------------------------------------------------- normalisation


def test_ranks_are_normalised_within_a_channel():
    """A rank means 'in the top fraction of what this channel ranked'."""
    assert agg.normalise_ranks(["INV1", "INV2", "INV3", "INV4"], 4) == {
        "INV1": 0.25, "INV2": 0.5, "INV3": 0.75, "INV4": 1.0}


def test_a_duplicate_key_keeps_its_best_rank():
    """Two formulations collapsing onto one identity must not penalise it."""
    assert agg.normalise_ranks(["INV1", "INV2", "INV1"], 3)["INV1"] == pytest.approx(1 / 3)


def test_an_empty_channel_normalises_to_nothing():
    assert agg.normalise_ranks([], 0) == {}


# ------------------------------------------------------------------- convergence


def _lists(**channels):
    return {name: list(keys) for name, keys in channels.items()}


def test_absence_from_a_channel_is_not_a_penalty():
    """Absence is missing evidence, not evidence against (channel E ranks only 8 drugs)."""
    results, _ = agg.aggregate(_lists(
        narrow=["INV1"],
        wide=["INV9", "INV8", "INV7", "INV6", "INV1"]))
    top = {a.key: a for a in results}["INV1"]
    assert top.ranks == {"narrow": 1.0, "wide": 1.0}
    assert set(top.channels) == {"narrow", "wide"}


def test_a_channel_ranking_most_of_the_pool_is_broad():
    """Channel B scored 2,566 of the drugs it could reach; agreeing with it means little."""
    wide = [f"INV{i}" for i in range(10)]
    results, stats = agg.aggregate(_lists(wide=wide, narrow=["INV0"]))
    assert stats["broad_channels"] == ["wide"]
    assert {a.key: a for a in results}["INV0"].convergence == "mixed"


def test_two_discriminating_channels_are_real_convergence():
    results, stats = agg.aggregate(_lists(
        a=["INV0", "INV1"], b=["INV0", "INV2"],
        wide=[f"INV{i}" for i in range(10)]))
    by_key = {x.key: x for x in results}
    assert stats["broad_channels"] == ["wide"]
    assert by_key["INV0"].convergence == "discriminating"
    assert by_key["INV1"].convergence == "mixed"


def test_support_from_broad_channels_only_is_labelled_as_such():
    wide_a = [f"INV{i}" for i in range(10)]
    wide_b = list(reversed(wide_a))
    results, stats = agg.aggregate(_lists(wa=wide_a, wb=wide_b))
    assert set(stats["broad_channels"]) == {"wa", "wb"}
    assert all(a.convergence == "broad_only" for a in results)


def test_a_single_channel_is_not_convergence():
    results, _ = agg.aggregate(_lists(only=["INV1", "INV2"]))
    assert all(a.convergence == "none" for a in results)


def test_convergence_breaks_a_tie_in_favour_of_the_discriminating_pair():
    """Two narrow channels agreeing must outrank the same score reached by one channel."""
    results, _ = agg.aggregate(_lists(a=["INV0"], b=["INV0"], c=["INV1"]))
    assert results[0].key == "INV0"
    assert results[0].convergence == "discriminating"


def test_the_order_is_deterministic_regardless_of_input_order():
    first, _ = agg.aggregate(_lists(a=["INV1", "INV2"], b=["INV2", "INV1"]))
    second, _ = agg.aggregate(_lists(b=["INV2", "INV1"], a=["INV1", "INV2"]))
    assert [x.key for x in first] == [x.key for x in second]


def test_stats_report_what_each_channel_contributed():
    _, stats = agg.aggregate(_lists(a=["INV0", "INV1"], b=["INV0"]))
    assert stats["per_channel_ranked"] == {"a": 2, "b": 1}
    assert stats["candidates"] == 2
    assert sum(stats["by_convergence"].values()) == 2


def test_every_candidate_records_how_many_channels_could_have_ranked_it():
    """n_channels_possible is what the Bonferroni correction uses; it must be the total."""
    results, _ = agg.aggregate(_lists(a=["INV0"], b=["INV1"], c=["INV2"]))
    assert all(x.n_channels_possible == 3 for x in results)
