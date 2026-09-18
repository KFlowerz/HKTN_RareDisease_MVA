"""L3 -- Rank aggregation across channels, with convergence reported honestly.

Purpose
    Combine the channels' independent rankings into one, preferring candidates that more
    than one channel found. The method is Robust Rank Aggregation (Kolde et al., 2012)
    [kolde2012] doi:10.1093/bioinformatics/btr709: each channel contributes a *normalised*
    rank in (0, 1], and a candidate's score asks how unlikely its rank vector would be if
    the channels had ranked at random.

Inputs
    One ranked list per channel, already harmonised onto a shared identity key
    (:mod:`.harmonize`), plus the length of each channel's list.

Outputs
    One :class:`Aggregate` per candidate: its RRA p-value, the corrected score used for
    ranking, its per-channel normalised ranks, which channels supported it, and what
    *kind* of convergence produced it.

Method
    For a candidate appearing in some channels, its normalised ranks are sorted ascending
    to r(1) ≤ … ≤ r(k). Under the null that a channel ranks at random, r(i) is distributed
    Beta(i, n - i + 1) where n is the number of channels the candidate *could* have
    appeared in. The RRA score is the smallest of those Beta CDFs, Bonferroni-corrected by
    n. A small score means: no random ranking would plausibly have placed this candidate
    this well, this often.

    A channel that did not rank a candidate contributes nothing rather than a zero. An
    absent rank is missing evidence, not evidence of absence -- Channel E ranks 8
    compounds from a curated file, and the 2,560 drugs it never considered are not
    evidence against them.

Guardrail
    **Convergence is not a count.** Channel B scored 2,566 of the approved drugs it could
    reach, so its agreeing with anything is close to uninformative, while Channel D and
    Channel E each nominate a few dozen from different evidence entirely. Reporting
    "supported by 2 channels" without saying *which two* would let the weakest possible
    agreement read as the strongest. Every row therefore carries
    :attr:`Aggregate.convergence`, and the discriminating channels are named in config
    rather than hard-coded, so the judgement is visible and arguable.

    **Channels are never re-weighted to produce a nicer ranking.** The selectivity weights
    below describe how much of the candidate pool each channel ranks, which is a property
    of the channel, not a preference for its results.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

LOGGER = logging.getLogger(__name__)

#: A channel that ranks most of the candidate pool cannot discriminate, whatever its
#: internal score says. Above this fraction a channel is "broad" and its agreement is not
#: counted as convergence on its own.
BROAD_CHANNEL_FRACTION = 0.5


@dataclass(frozen=True)
class Aggregate:
    """One candidate's aggregated position."""

    key: str
    score: float                 # corrected RRA score; smaller is better
    p_value: float               # uncorrected minimum Beta CDF
    ranks: dict                  # {channel: normalised rank in (0, 1]}
    channels: tuple              # channels that ranked it, sorted
    discriminating: tuple        # of those, the ones that are not "broad"
    n_channels_possible: int

    @property
    def n_channels_supporting(self) -> int:
        return len(self.channels)

    @property
    def convergence(self) -> str:
        """What kind of agreement produced this row.

        ``none`` -- one channel only.
        ``broad_only`` -- more than one channel, but only broad ones; weak evidence.
        ``discriminating`` -- at least two channels that each rank a small slice.
        ``mixed`` -- one discriminating channel plus a broad one.
        """
        if self.n_channels_supporting < 2:
            return "none"
        if len(self.discriminating) >= 2:
            return "discriminating"
        if len(self.discriminating) == 1:
            return "mixed"
        return "broad_only"


def _beta_cdf(x: float, a: int, b: int) -> float:
    """Regularised incomplete beta I_x(a, b) for positive integer shapes.

    For integer ``a`` and ``b`` this is exactly the binomial sum
    ``sum_{j=a}^{a+b-1} C(n, j) x^j (1-x)^(n-j)`` with ``n = a + b - 1``, which is
    computed here directly. Written out rather than taken from SciPy so the whole method
    is inspectable in this file and the pipeline keeps one fewer dependency in the path
    that produces its headline ranking.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    n = a + b - 1
    return sum(math.comb(n, j) * x ** j * (1.0 - x) ** (n - j) for j in range(a, n + 1))


def rho_score(normalised_ranks, n_possible: int) -> tuple:
    """Robust Rank Aggregation score. Returns ``(corrected_score, p_value)``.

    Args:
        normalised_ranks: this candidate's ranks in (0, 1], one per channel that ranked it.
        n_possible: how many channels could have ranked it.

    The minimum is taken over prefixes: a candidate ranked superbly by one channel and
    ignored by the rest scores on that one, and a candidate ranked decently by several
    scores on the group. Kolde's Bonferroni correction by ``n_possible`` keeps the two
    comparable.
    """
    ordered = sorted(normalised_ranks)
    if not ordered:
        return 1.0, 1.0
    n = max(n_possible, len(ordered))
    p = min(_beta_cdf(r, i, n - i + 1) for i, r in enumerate(ordered, start=1))
    return min(1.0, p * n), p


def normalise_ranks(table, total: int) -> dict:
    """``{key: normalised rank}`` from an ordered list of keys.

    ``total`` is the length of the channel's own list, so a rank is read as "this
    candidate was in the top r fraction of what this channel ranked". Ties are not
    expected -- each channel emits a strict order -- and a duplicate key keeps its best
    (smallest) rank, which is what happens when two formulations collapse onto one
    identity.
    """
    if total <= 0:
        return {}
    out: dict = {}
    for position, key in enumerate(table, start=1):
        value = position / total
        if key not in out or value < out[key]:
            out[key] = value
    return out


def aggregate(per_channel: dict, *, broad_fraction: float = BROAD_CHANNEL_FRACTION,
              pool_size: int = 0) -> tuple:
    """Aggregate ranked lists. Returns ``(sorted Aggregates, stats)``.

    Args:
        per_channel: ``{channel: [key, ...]}`` in each channel's own rank order, already
            harmonised so a key means the same molecule everywhere.
        broad_fraction: a channel ranking more than this fraction of the pool is "broad"
            and its support alone is not treated as convergence.
        pool_size: the number of distinct candidates across all channels; defaults to the
            observed union.
    """
    channels = sorted(per_channel)
    normalised = {c: normalise_ranks(per_channel[c], len(per_channel[c])) for c in channels}
    union = sorted({k for ranks in normalised.values() for k in ranks})
    pool = pool_size or len(union)

    broad = {c for c in channels if pool and len(normalised[c]) / pool > broad_fraction}
    if broad:
        LOGGER.info("broad channel(s) ranking more than %.0f%% of the pool: %s",
                    broad_fraction * 100, ", ".join(sorted(broad)))

    results = []
    for key in union:
        supporting = tuple(c for c in channels if key in normalised[c])
        score, p = rho_score([normalised[c][key] for c in supporting], len(channels))
        results.append(Aggregate(
            key=key, score=score, p_value=p,
            ranks={c: round(normalised[c][key], 6) for c in supporting},
            channels=supporting,
            discriminating=tuple(c for c in supporting if c not in broad),
            n_channels_possible=len(channels)))

    # Ties on score break on convergence quality, then support count, then key -- so a
    # candidate two discriminating channels found outranks one a broad channel happened to
    # place identically, and the order never depends on dict iteration.
    order = {"discriminating": 0, "mixed": 1, "broad_only": 2, "none": 3}
    results.sort(key=lambda a: (a.score, order[a.convergence], -a.n_channels_supporting,
                                a.key))

    stats = {
        "channels": channels,
        "broad_channels": sorted(broad),
        "broad_fraction": broad_fraction,
        "candidates": len(results),
        "by_convergence": {kind: sum(1 for a in results if a.convergence == kind)
                           for kind in ("discriminating", "mixed", "broad_only", "none")},
        "per_channel_ranked": {c: len(normalised[c]) for c in channels},
    }
    return results, stats
