"""Weighted connectivity scoring -- the numerics behind Channel C.

Purpose
    Score how strongly a perturbation's transcriptional response *reverses* a query
    signature, by the weighted connectivity score (WTCS) of the L1000 Connectivity Map
    [subramanian2017] doi:10.1016/j.cell.2017.10.049 -- a two-tailed, weighted
    Kolmogorov-Smirnov statistic over the query's up- and down-regulated gene sets, built
    on the running-sum enrichment score of GSEA [subramanian2005]
    doi:10.1073/pnas.0506580102.

Inputs
    A matrix of differential-expression z-scores (signatures x genes) and two disjoint
    query gene sets, given as boolean masks over the gene axis.

Outputs
    Per signature: an enrichment score for each tail, the combined WTCS, its
    within-group normalised form (NCS), an aggregation across the signatures of one
    perturbagen, and an empirical p against resampled queries.

Guardrail
    **No data reaches this module at all** -- patient-derived or otherwise. It takes
    arrays and returns arrays; it opens no file, reads no config, and knows nothing about
    drugs, cell lines or genes. That is deliberate: it is the part of Channel C that can
    be tested exhaustively on invented numbers, so every claim the channel makes about
    its own method is checked here rather than asserted in a docstring.

    **A connectivity score is not a p-value.** WTCS and NCS are *relative* measures with
    no null of their own -- NCS divides by a mean magnitude, which rescales but never
    calibrates. :func:`empirical_p` supplies the null by resampling the query, and the
    channel must not report a ranking without one. This is the rule Channel B already
    follows for proximity: a raw distance without a degree-matched null is a popularity
    measurement, not evidence.
"""

from __future__ import annotations

import numpy as np

#: CMap's aggregation quantiles across the signatures of one perturbagen. The score of
#: larger magnitude wins, so a compound that reverses strongly in a few cell lines is not
#: averaged away by the lines where it does nothing.
#: [subramanian2017] doi:10.1016/j.cell.2017.10.049
MAX_QUANTILE = (0.33, 0.67)


def enrichment_score(values: np.ndarray, member_mask: np.ndarray) -> tuple:
    """Weighted running-sum enrichment score of one gene set in each signature.

    The statistic walks each signature's genes from most up-regulated to most
    down-regulated, adding ``|z|`` at every query gene and subtracting a constant step at
    every other, then takes the running sum's greatest excursion from zero. Weighting by
    ``|z|`` is what makes it *weighted*: a query gene at the very top of the ranking
    counts for more than one that merely appears in the top half.

    Args:
        values: ``(n_signatures, n_genes)`` z-scores. Peak memory is a small multiple of
            the input, so the caller feeds this in chunks.
        member_mask: ``(n_genes,)`` boolean -- the query set.

    Returns:
        ``(es, usable)``. ``es`` is ``(n_signatures,)`` and signed: positive when the set
        concentrates among the up-regulated genes, negative among the down-regulated.
        ``usable`` is False for signatures whose query weights summed to zero, where the
        statistic is undefined; their score is set to 0 rather than to a NaN that would
        propagate silently into a ranking.

    Raises:
        ValueError: If the mask is empty or covers every gene. Both leave the miss step
            undefined, and both mean the caller built the query wrongly.
    """
    values = np.asarray(values, dtype=np.float64)
    # Descending by z. A stable sort makes ties resolve by gene index, so the score does
    # not depend on the platform's sort implementation -- this pipeline is seeded and its
    # results have to reproduce.
    return _score_in_order(values, np.argsort(-values, axis=1, kind="stable"), member_mask)


def _score_in_order(values: np.ndarray, order: np.ndarray,
                    member_mask: np.ndarray) -> tuple:
    """:func:`enrichment_score` with the descending order already computed.

    The order depends only on the signature, never on the query, so the permutation null
    -- which rescores one fixed set of signatures against thousands of random queries --
    sorts once instead of once per draw.
    """
    member_mask = np.asarray(member_mask, dtype=bool)
    n_genes = values.shape[1]
    n_members = int(member_mask.sum())
    if n_members == 0 or n_members >= n_genes:
        raise ValueError(
            f"query set covers {n_members} of {n_genes} gene(s); an empty set has no "
            "enrichment and a full set leaves no background to be enriched against"
        )

    hits = member_mask[order]
    weights = np.abs(np.take_along_axis(values, order, axis=1)) * hits

    total = weights.sum(axis=1, keepdims=True)
    usable = total[:, 0] > 0
    # Divide by 1 where the total is zero; those rows are zeroed out below anyway.
    hit_cum = np.cumsum(weights, axis=1) / np.where(total > 0, total, 1.0)
    miss_cum = np.cumsum(~hits, axis=1) / float(n_genes - n_members)

    deviation = hit_cum - miss_cum
    extreme = np.argmax(np.abs(deviation), axis=1)
    es = deviation[np.arange(values.shape[0]), extreme]
    return np.where(usable, es, 0.0), usable


def descending_order(values: np.ndarray) -> np.ndarray:
    """Each signature's genes ordered most up-regulated first.

    Exposed so a caller scoring one fixed set of signatures against many queries -- the
    permutation null -- can sort once and pass the result back in.
    """
    return np.argsort(-np.asarray(values, dtype=np.float64), axis=1, kind="stable")


def weighted_connectivity(values: np.ndarray, up_mask: np.ndarray,
                          down_mask: np.ndarray, *, order: np.ndarray = None) -> tuple:
    """WTCS of each signature against an up/down query.

    Positive means the signature *resembles* the query; negative means it **reverses**
    it, which is what Channel C is looking for.

    The two tails must disagree for the score to be non-zero. A perturbation that pushes
    the query's up-genes and its down-genes in the *same* direction has not reversed
    anything -- it has moved the whole query set together, which is a property of the
    signature's overall shape rather than of the contrast the query encodes. CMap scores
    that case as 0, and so does this.

    Args:
        values: ``(n_signatures, n_genes)`` z-scores.
        up_mask: ``(n_genes,)`` boolean -- query genes up in the proxy state.
        down_mask: ``(n_genes,)`` boolean -- query genes down in it.
        order: an optional descending order from :func:`descending_order`, reused across
            queries over the same signatures.

    Returns:
        ``(wtcs, es_up, es_down, usable)``.

    Raises:
        ValueError: If the two query sets overlap. A gene cannot be both tails of one
            contrast, and an overlap would double-count it in both running sums.
    """
    values = np.asarray(values, dtype=np.float64)
    up_mask = np.asarray(up_mask, dtype=bool)
    down_mask = np.asarray(down_mask, dtype=bool)
    overlap = int(np.logical_and(up_mask, down_mask).sum())
    if overlap:
        raise ValueError(f"{overlap} gene(s) appear in both query tails; the up and down "
                         "sets of one contrast are disjoint by construction")

    if order is None:
        order = descending_order(values)
    es_up, usable_up = _score_in_order(values, order, up_mask)
    es_down, usable_down = _score_in_order(values, order, down_mask)
    opposed = np.sign(es_up) != np.sign(es_down)
    wtcs = np.where(opposed, (es_up - es_down) / 2.0, 0.0)
    return wtcs, es_up, es_down, np.logical_and(usable_up, usable_down)


def normalise(wtcs: np.ndarray, groups) -> np.ndarray:
    """Normalised connectivity: WTCS divided by its mean magnitude within group and sign.

    Signatures from different cell lines are not on one scale -- some lines respond to
    almost everything and produce large scores throughout, which would let the cell line,
    rather than the drug, decide the ranking. CMap divides each score by the mean
    magnitude of the scores sharing its cell line *and* its sign, separately because the
    positive and negative tails are not symmetric.
    [subramanian2017] doi:10.1016/j.cell.2017.10.049

    Args:
        wtcs: ``(n_signatures,)`` connectivity scores.
        groups: ``(n_signatures,)`` group label per signature -- the cell line.

    Returns:
        ``(n_signatures,)`` normalised scores. A signature whose group-and-sign stratum
        has zero mean magnitude keeps a score of 0: there is nothing to normalise
        against, and inventing a divisor would manufacture a large score out of noise.
    """
    wtcs = np.asarray(wtcs, dtype=np.float64)
    groups = np.asarray(groups)
    ncs = np.zeros_like(wtcs)
    for group in np.unique(groups):
        in_group = groups == group
        for positive in (True, False):
            stratum = in_group & ((wtcs > 0) if positive else (wtcs < 0))
            if not stratum.any():
                continue
            scale = float(np.abs(wtcs[stratum]).mean())
            if scale > 0:
                ncs[stratum] = wtcs[stratum] / scale
    return ncs


def max_quantile(values) -> float:
    """Aggregate one perturbagen's signatures: the more extreme of two quantiles.

    A mean would let the cell lines where a drug does nothing cancel the ones where it
    reverses the query strongly, and a maximum would report the single luckiest
    signature. CMap's compromise takes the 33rd and 67th percentiles and keeps whichever
    is larger in magnitude, so a consistent effect in either direction survives while a
    lone outlier does not.
    [subramanian2017] doi:10.1016/j.cell.2017.10.049

    Returns:
        The aggregated score, or ``0.0`` for an empty input.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    low, high = np.percentile(values, [MAX_QUANTILE[0] * 100, MAX_QUANTILE[1] * 100])
    return float(high if abs(high) >= abs(low) else low)


def query_masks(consensus: np.ndarray, size: int) -> tuple:
    """The query's two tails: the ``size`` most up- and most down-regulated genes.

    Args:
        consensus: ``(n_genes,)`` z-scores of the proxy signature.
        size: how many genes each tail takes.

    Returns:
        ``(up_mask, down_mask)``, disjoint boolean arrays.

    Raises:
        ValueError: If the two tails would meet or overlap. Taking more than half the
            genes into each tail leaves no background, and the score would then be
            measuring the partition rather than the signature.
    """
    consensus = np.asarray(consensus, dtype=np.float64)
    n_genes = consensus.size
    if size < 1 or 2 * size >= n_genes:
        raise ValueError(
            f"query size {size} against {n_genes} gene(s): each tail needs at least one "
            "gene, and the two together must leave a background to be enriched against"
        )
    order = np.argsort(-consensus, kind="stable")
    up_mask = np.zeros(n_genes, dtype=bool)
    down_mask = np.zeros(n_genes, dtype=bool)
    up_mask[order[:size]] = True
    down_mask[order[-size:]] = True
    return up_mask, down_mask


def random_query_masks(n_genes: int, size: int, rng) -> tuple:
    """A random query of the same shape -- the null :func:`query_masks` is tested against.

    The null asks whether *this* gene set reverses the signature more than an arbitrary
    set of the same size would. It therefore resamples the query, not the signature: the
    signature's own structure (its overall spread, its tail shape) is a property of the
    experiment and has to be held fixed, or the comparison stops being about the query.

    Args:
        n_genes: the gene axis both tails are drawn from.
        size: genes per tail, matching the observed query.
        rng: a seeded :class:`numpy.random.Generator`.
    """
    if 2 * size >= n_genes:
        raise ValueError(f"cannot draw two disjoint tails of {size} from {n_genes} gene(s)")
    drawn = rng.choice(n_genes, size=2 * size, replace=False)
    up_mask = np.zeros(n_genes, dtype=bool)
    down_mask = np.zeros(n_genes, dtype=bool)
    up_mask[drawn[:size]] = True
    down_mask[drawn[size:]] = True
    return up_mask, down_mask


def empirical_p(observed: float, null, *, tail: str = "lower") -> float:
    """One-sided empirical p, with the observation counted in its own null.

    The ``+1`` in both terms is not a rounding choice: leaving it out lets a p of exactly
    0 be reported, which claims more resolution than ``len(null)`` draws can carry. The
    smallest attainable value is ``1 / (len(null) + 1)``, and the channel records the
    permutation count beside every p so that floor stays visible.

    Args:
        observed: the score to test.
        null: scores from resampled queries.
        tail: ``"lower"`` for reversal, where a strongly negative score is the extreme
            one; ``"upper"`` for resemblance.

    Raises:
        ValueError: If the null is empty, or ``tail`` is neither value.
    """
    null = np.asarray(null, dtype=np.float64)
    if null.size == 0:
        raise ValueError("an empirical p needs a null; none was supplied")
    if tail not in ("lower", "upper"):
        raise ValueError(f"tail must be 'lower' or 'upper', not {tail!r}")
    extreme = (null <= observed) if tail == "lower" else (null >= observed)
    return float((int(extreme.sum()) + 1) / (null.size + 1))
