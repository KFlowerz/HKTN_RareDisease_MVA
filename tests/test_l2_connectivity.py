"""Tests for the connectivity numerics behind Channel C.

Every number here is invented. The point of :mod:`src.l2_channels.connectivity` being
data-free is that its claims can be checked against constructed cases rather than
asserted in a docstring, so this file constructs the cases: a signature that perfectly
reverses a query, one that perfectly matches it, one that does neither, and the
degenerate inputs that would otherwise produce a confident number from nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.l2_channels import connectivity

N_GENES = 20


def _masks(up, down, n_genes=N_GENES):
    """Boolean masks from index lists."""
    up_mask = np.zeros(n_genes, dtype=bool)
    down_mask = np.zeros(n_genes, dtype=bool)
    up_mask[list(up)] = True
    down_mask[list(down)] = True
    return up_mask, down_mask


def _signature(high, low, n_genes=N_GENES):
    """A signature with ``high`` genes strongly up and ``low`` strongly down."""
    values = np.zeros(n_genes)
    values[list(high)] = 3.0
    values[list(low)] = -3.0
    return values.reshape(1, -1)


# ------------------------------------------------------------ the enrichment score


def test_a_set_at_the_top_scores_positive() -> None:
    values = _signature(high=[0, 1, 2], low=[17, 18, 19])
    mask, _ = _masks([0, 1, 2], [])
    es, usable = connectivity.enrichment_score(values, mask)
    assert usable[0]
    assert es[0] > 0.9


def test_a_set_at_the_bottom_scores_negative() -> None:
    values = _signature(high=[0, 1, 2], low=[17, 18, 19])
    mask, _ = _masks([17, 18, 19], [])
    es, _ = connectivity.enrichment_score(values, mask)
    assert es[0] < -0.9


def test_a_spread_set_scores_far_below_a_concentrated_one() -> None:
    """Stated as a comparison, because the statistic is not centred on zero.

    A maximum-deviation statistic over a 20-gene space is positively biased even for a
    set drawn at random -- the running sum gets many chances to wander. The claim that
    matters is the contrast: concentration at the top is what the score responds to.
    """
    values = np.arange(N_GENES, dtype=float)[::-1].reshape(1, -1)
    spread, _ = _masks(range(0, N_GENES, 2), [])
    concentrated, _ = _masks(range(0, 10), [])
    spread_es, _ = connectivity.enrichment_score(values, spread)
    concentrated_es, _ = connectivity.enrichment_score(values, concentrated)
    assert spread_es[0] < 0.5 * concentrated_es[0]


def test_an_all_zero_query_is_unusable_and_scores_zero() -> None:
    """The statistic is undefined with no weight; it must not become a NaN in a ranking."""
    values = np.zeros((1, N_GENES))
    values[0, 5] = 1.0
    mask, _ = _masks([0, 1], [])          # both weights are zero
    es, usable = connectivity.enrichment_score(values, mask)
    assert not usable[0]
    assert es[0] == 0.0
    assert not np.isnan(es[0])


@pytest.mark.parametrize("members", [[], list(range(N_GENES))])
def test_an_empty_or_full_query_is_refused(members) -> None:
    values = _signature(high=[0], low=[19])
    mask, _ = _masks(members, [])
    with pytest.raises(ValueError, match="background"):
        connectivity.enrichment_score(values, mask)


def test_the_score_is_weighted_not_just_positional() -> None:
    """Two sets at the same ranks score differently when their magnitudes differ.

    This is the difference between weighted and classical GSEA, and the whole reason the
    z-scores are carried through rather than reduced to ranks.

    Both signatures put the same two query genes at the same two positions -- ranks 0 and
    8 -- so a purely positional statistic would score them identically. They differ only
    in how the weight is divided between those two members.
    """
    steep = np.array([[10.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.2, -1.0, -2.0]])
    flat = np.array([[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, -1.0, -2.0]])
    mask = np.zeros(10, dtype=bool)
    mask[[0, 8]] = True
    steep_es, _ = connectivity.enrichment_score(steep, mask)
    flat_es, _ = connectivity.enrichment_score(flat, mask)
    assert steep_es[0] > flat_es[0] + 0.3


def test_precomputed_order_gives_the_same_score() -> None:
    """The null's optimisation must not change the statistic."""
    rng = np.random.default_rng(0)
    values = rng.normal(size=(8, N_GENES))
    up_mask, down_mask = _masks([0, 1, 2], [10, 11, 12])
    plain = connectivity.weighted_connectivity(values, up_mask, down_mask)[0]
    reused = connectivity.weighted_connectivity(
        values, up_mask, down_mask, order=connectivity.descending_order(values))[0]
    assert np.allclose(plain, reused)


# ------------------------------------------------------------ WTCS


def test_a_reversing_signature_scores_negative() -> None:
    """The query's up-genes are down in the signature and vice versa."""
    values = _signature(high=[10, 11, 12], low=[0, 1, 2])
    up_mask, down_mask = _masks([0, 1, 2], [10, 11, 12])
    wtcs, es_up, es_down, usable = connectivity.weighted_connectivity(
        values, up_mask, down_mask)
    assert usable[0]
    assert es_up[0] < 0 < es_down[0]
    assert wtcs[0] < -0.9


def test_a_matching_signature_scores_positive() -> None:
    values = _signature(high=[0, 1, 2], low=[10, 11, 12])
    up_mask, down_mask = _masks([0, 1, 2], [10, 11, 12])
    wtcs, _, _, _ = connectivity.weighted_connectivity(values, up_mask, down_mask)
    assert wtcs[0] > 0.9


def test_both_tails_moving_together_scores_zero() -> None:
    """Not a reversal and not a match -- the contrast the query encodes did not move."""
    values = _signature(high=[0, 1, 2, 10, 11, 12], low=[17, 18, 19])
    up_mask, down_mask = _masks([0, 1, 2], [10, 11, 12])
    wtcs, es_up, es_down, _ = connectivity.weighted_connectivity(values, up_mask, down_mask)
    assert np.sign(es_up[0]) == np.sign(es_down[0])
    assert wtcs[0] == 0.0


def test_overlapping_query_tails_are_refused() -> None:
    up_mask, down_mask = _masks([0, 1, 2], [2, 3, 4])
    with pytest.raises(ValueError, match="both query tails"):
        connectivity.weighted_connectivity(_signature([0], [19]), up_mask, down_mask)


def test_wtcs_is_bounded_by_one() -> None:
    rng = np.random.default_rng(7)
    values = rng.normal(size=(200, N_GENES))
    up_mask, down_mask = _masks([0, 1, 2], [10, 11, 12])
    wtcs, _, _, _ = connectivity.weighted_connectivity(values, up_mask, down_mask)
    assert np.all(np.abs(wtcs) <= 1.0)


# ------------------------------------------------------------ normalisation


def test_normalisation_puts_cell_lines_on_one_scale() -> None:
    """A line that scores large on everything must not outrank one that scores small."""
    wtcs = np.array([-0.8, -0.4, -0.08, -0.04])
    groups = np.array(["LOUD", "LOUD", "QUIET", "QUIET"])
    ncs = connectivity.normalise(wtcs, groups)
    assert ncs[0] == pytest.approx(ncs[2])
    assert ncs[1] == pytest.approx(ncs[3])


def test_normalisation_keeps_the_sign() -> None:
    wtcs = np.array([-0.5, 0.5, -0.1, 0.1])
    ncs = connectivity.normalise(wtcs, np.array(["A", "A", "A", "A"]))
    assert np.all(np.sign(ncs) == np.sign(wtcs))


def test_normalisation_treats_the_two_signs_separately() -> None:
    """One strong reverser must not be rescaled by a crowd of weak mimics."""
    wtcs = np.array([-0.6, 0.02, 0.02, 0.02])
    ncs = connectivity.normalise(wtcs, np.array(["A"] * 4))
    assert ncs[0] == pytest.approx(-1.0)


def test_a_zero_score_stays_zero() -> None:
    """Nothing to normalise against must not become a large number."""
    ncs = connectivity.normalise(np.array([0.0, 0.0]), np.array(["A", "A"]))
    assert np.all(ncs == 0.0)


# ------------------------------------------------------------ aggregation


def test_max_quantile_keeps_a_consistent_effect() -> None:
    assert connectivity.max_quantile([-0.9, -0.8, -0.85, -0.95]) < -0.8


def test_max_quantile_discards_a_lone_outlier() -> None:
    """One lucky signature among many nulls must not carry the molecule."""
    values = [-0.99] + [0.0] * 9
    assert connectivity.max_quantile(values) == pytest.approx(0.0, abs=0.05)


def test_max_quantile_of_nothing_is_zero() -> None:
    assert connectivity.max_quantile([]) == 0.0


def test_max_quantile_prefers_the_larger_magnitude() -> None:
    """A mean would cancel these; the aggregation reports the side that moved."""
    values = [-0.9, -0.9, -0.9, -0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    assert connectivity.max_quantile(values) < 0


# ------------------------------------------------------------ the query and its null


def test_query_masks_take_the_extremes_and_stay_disjoint() -> None:
    consensus = np.arange(N_GENES, dtype=float)
    up_mask, down_mask = connectivity.query_masks(consensus, 3)
    assert set(np.flatnonzero(up_mask)) == {17, 18, 19}
    assert set(np.flatnonzero(down_mask)) == {0, 1, 2}
    assert not np.logical_and(up_mask, down_mask).any()


@pytest.mark.parametrize("size", [0, N_GENES // 2, N_GENES])
def test_a_query_leaving_no_background_is_refused(size) -> None:
    with pytest.raises(ValueError, match="background"):
        connectivity.query_masks(np.arange(N_GENES, dtype=float), size)


def test_random_query_masks_are_disjoint_and_the_right_size() -> None:
    rng = np.random.default_rng(3)
    up_mask, down_mask = connectivity.random_query_masks(N_GENES, 4, rng)
    assert up_mask.sum() == down_mask.sum() == 4
    assert not np.logical_and(up_mask, down_mask).any()


def test_random_query_masks_are_reproducible_from_the_seed() -> None:
    first = connectivity.random_query_masks(N_GENES, 4, np.random.default_rng(42))
    second = connectivity.random_query_masks(N_GENES, 4, np.random.default_rng(42))
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])


def test_empirical_p_counts_the_observation_in_its_own_null() -> None:
    """A p of exactly zero claims more resolution than the draws can carry."""
    p = connectivity.empirical_p(-1.0, np.zeros(99), tail="lower")
    assert p == pytest.approx(1 / 100)
    assert p > 0


def test_empirical_p_of_an_unremarkable_score_is_large() -> None:
    rng = np.random.default_rng(1)
    null = rng.normal(size=999)
    assert connectivity.empirical_p(0.0, null, tail="lower") > 0.4


def test_empirical_p_needs_a_null() -> None:
    with pytest.raises(ValueError, match="needs a null"):
        connectivity.empirical_p(-1.0, [])


def test_empirical_p_refuses_an_unknown_tail() -> None:
    with pytest.raises(ValueError, match="tail must be"):
        connectivity.empirical_p(-1.0, np.zeros(10), tail="sideways")


# ------------------------------------------------------------ end to end, on invented data


def test_a_planted_reverser_beats_the_null_and_noise_does_not() -> None:
    """The whole chain: score, normalise, aggregate, test.

    Two molecules over a 200-gene space -- one built to reverse the query, one drawn from
    noise. The planted one must come out significant and the noise one must not, or the
    pipeline would be reporting rank order as evidence.
    """
    rng = np.random.default_rng(11)
    n_genes, size = 200, 20
    consensus = rng.normal(size=n_genes)
    up_mask, down_mask = connectivity.query_masks(consensus, size)

    reverser = np.vstack([-consensus + rng.normal(scale=0.3, size=n_genes)
                          for _ in range(4)])
    noise = rng.normal(size=(4, n_genes))
    values = np.vstack([reverser, noise])
    order = connectivity.descending_order(values)
    wtcs, _, _, _ = connectivity.weighted_connectivity(values, up_mask, down_mask,
                                                       order=order)
    ncs = connectivity.normalise(wtcs, np.array(["LINE"] * 8))
    assert connectivity.max_quantile(ncs[:4]) < connectivity.max_quantile(ncs[4:])

    null_rng = np.random.default_rng(99)
    draws = []
    for _ in range(199):
        random_up, random_down = connectivity.random_query_masks(n_genes, size, null_rng)
        drawn, _, _, _ = connectivity.weighted_connectivity(values, random_up,
                                                            random_down, order=order)
        draws.append((connectivity.max_quantile(drawn[:4]),
                      connectivity.max_quantile(drawn[4:])))
    planted_null = [d[0] for d in draws]
    noise_null = [d[1] for d in draws]

    assert connectivity.empirical_p(connectivity.max_quantile(wtcs[:4]),
                                    planted_null, tail="lower") <= 0.05
    assert connectivity.empirical_p(connectivity.max_quantile(wtcs[4:]),
                                    noise_null, tail="lower") > 0.05
