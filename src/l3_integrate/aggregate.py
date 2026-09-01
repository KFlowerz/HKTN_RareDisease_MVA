"""L3 -- Rank aggregation across channels.

Purpose
    Combine the per-channel rankings into one consensus ranking, preferring candidates
    that **multiple independent channels** surface. Robust Rank Aggregation gives a
    p-value against the null that a candidate's ranks are randomly distributed; Borda
    is carried as a simpler, more transparent comparator.

Inputs
    RxCUI-harmonized per-channel ranked tables from ``harmonize``.

Outputs
    One ranked table with per-channel ranks retained, the aggregated score, an
    ``n_channels_supporting`` count, and the aggregation method used.

Guardrail
    Cross-channel convergence is the point of the whole architecture: a candidate ranked
    modestly by four independent channels is stronger evidence than one ranked first by a
    single channel. The aggregation must reward that, and ``n_channels_supporting`` must
    be visible on every downstream row.

    Handle **partial lists honestly**. Channels cover different drug universes, so a
    candidate's absence from a channel usually means "not evaluated", not "ranked last".
    Imputing a worst rank for absence would systematically punish drugs that only one
    resource covers.

    Never weight channels to make a preferred compound win. Weights, if any, must be set
    before seeing results and recorded in the output.
"""

from __future__ import annotations


def rank_aggregate(config: dict) -> None:
    """Aggregate per-channel rankings into a consensus ranking.

    Args:
        config: Parsed pipeline configuration; uses ``channels``, ``results_dir``,
            ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: normalize each channel's ranks to [0, 1]; run Robust Rank Aggregation
    # (beta-distribution order statistics) to get a per-candidate p-value, with an
    # explicit rule for partial lists -- restrict each channel's normalization to the
    # candidates it actually evaluated rather than imputing a worst rank; compute Borda
    # as a comparator and report BOTH so the ranking's sensitivity to method is visible;
    # correct for multiple testing (BH) across candidates.
    #
    # Emit per-channel ranks, `n_channels_supporting`, RRA p/q, Borda score, and the
    # aggregation parameters on every row.
    raise NotImplementedError("aggregate.rank_aggregate is a scaffold stub")
