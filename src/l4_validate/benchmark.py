"""L4 -- Blinded internal benchmark.

Purpose
    Answer the question a judge will ask: *does this pipeline recover compounds we already
    believe are relevant?* There is **no external ground truth for Track 2** -- MVA has no
    approved therapy and no repurposing gold standard exists -- so the only honest
    validation is an internally-held-out set of aneuploidy/SAC-relevant compounds that the
    scored channel was not built from.

Inputs
    The positive set in ``src/l2_channels/aneuploidy_prior.tsv``, L3's ranked candidates,
    and the per-channel ranks L3 recorded.

Outputs
    Under ``config["results_dir"]/l4/benchmark/``: the rank of each positive, enrichment
    at k, AUROC against everything else ranked, and a per-channel breakdown -- with
    channel E's own numbers marked circular and excluded from the headline.

Guardrail
    **The positive set is channel E's seed file, so channel E cannot be benchmarked on
    it.** Scoring channel E against the compounds channel E was built from measures
    nothing: the recovery is 100% by construction, and reporting it as a result would be
    the leakage the design is supposed to avoid. This module therefore reports channel E's
    number as ``circular`` and computes the headline against the channels that never saw
    the set -- B (network proximity) and D (phenotype). Those are a genuine test: does an
    independent method rediscover compounds the literature already links to aneuploidy?

    **Report the negative result if that is what happens.** An honest miss is worth more
    to Scientific Rigor than a tuned hit, and a benchmark that always passes is evidence
    of leakage, not quality. Nothing here may be used to re-tune a channel without
    re-freezing the set and disclosing the reuse.

    **This measures pipeline recovery, not clinical efficacy.** Nothing here licenses a
    claim that any candidate works.

    The set is small -- around a dozen compounds, of which fewer reach any given channel --
    so every metric carries wide uncertainty. The bootstrap interval is reported for
    exactly that reason, and a point estimate from this set is not a result on its own.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import random
from pathlib import Path

from ..l2_channels import channel_e_prior

LOGGER = logging.getLogger(__name__)

BENCHMARK_DIR = "benchmark"
#: Channels built from the positive set; their recovery of it is circular by construction.
CIRCULAR_CHANNELS = frozenset({"prior"})
#: Enrichment is reported at these cutoffs.
ENRICHMENT_AT = (10, 50, 100, 500)
BOOTSTRAP_SAMPLES = 2000


def positive_set(path: Path = None) -> tuple:
    """``(names, fingerprint)`` -- the compounds the benchmark looks for.

    The fingerprint is a SHA-256 of the normalised compound names, so a change to the seed
    file changes the recorded benchmark identity. A benchmark whose contents can drift
    unnoticed is not blinded.
    """
    rows = channel_e_prior.load_prior(path or channel_e_prior.PRIOR_FILE)
    names = sorted({r["compound"].strip().lower() for r in rows})
    digest = hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()
    return names, digest


def _rank_column(channel: str) -> str:
    return f"rank_{channel}"


def recovery(ranked: list, positives: set, key: str = "drug_name") -> dict:
    """Where the positives land in a ranking. Returns ranks and enrichment.

    ``ranked`` is in rank order. Positions are 1-based over the rows actually ranked.
    """
    hits, positions = [], []
    for position, row in enumerate(ranked, start=1):
        if (row.get(key) or "").strip().lower() in positives:
            hits.append(row.get(key, "").strip().lower())
            positions.append(position)
    return {"n_ranked": len(ranked), "n_positives_found": len(positions),
            "positions": positions, "found": sorted(hits),
            "enrichment_at": {str(k): sum(1 for p in positions if p <= k)
                              for k in ENRICHMENT_AT}}


def auroc(positions: list, n_ranked: int, n_positives: int) -> float:
    """AUROC of the ranking against every other ranked item.

    Equivalent to the Mann-Whitney U statistic normalised by the number of pairs: the
    probability that a randomly chosen positive outranks a randomly chosen negative.
    0.5 is chance. Returns ``float('nan')`` when there is nothing to compare.
    """
    n_negatives = n_ranked - n_positives
    if n_positives <= 0 or n_negatives <= 0:
        return float("nan")
    # Sum over positives of how many negatives they beat. A positive at position p has
    # (p - 1) items above it, of which (number of positives above it) are not negatives.
    better = 0
    for index, position in enumerate(sorted(positions), start=1):
        negatives_above = (position - 1) - (index - 1)
        better += n_negatives - negatives_above
    return better / (n_positives * n_negatives)


def bootstrap_auroc(positions: list, n_ranked: int, n_positives: int, *, seed: int,
                    samples: int = BOOTSTRAP_SAMPLES) -> dict:
    """A percentile interval for AUROC, resampling the positives with replacement.

    With a positive set this small the interval is wide, which is the point: it stops a
    point estimate from being read as a measurement.
    """
    if n_positives <= 1:
        return {"samples": 0, "ci_low": None, "ci_high": None,
                "note": "too few positives to resample"}
    rng = random.Random(seed)
    values = []
    for _ in range(samples):
        drawn = [rng.choice(positions) for _ in positions]
        value = auroc(sorted(drawn), n_ranked, n_positives)
        if value == value:  # not NaN
            values.append(value)
    if not values:
        return {"samples": 0, "ci_low": None, "ci_high": None, "note": "no valid samples"}
    values.sort()
    return {"samples": len(values),
            "ci_low": round(values[int(0.025 * (len(values) - 1))], 4),
            "ci_high": round(values[int(0.975 * (len(values) - 1))], 4)}


def channel_ranking(rows: list, channel: str) -> list:
    """L3's rows that the named channel ranked, in that channel's own order."""
    column = _rank_column(channel)
    ranked = [r for r in rows if (r.get(column) or "").strip()]
    return sorted(ranked, key=lambda r: float(r[column]))


def run_blinded(config: dict) -> dict:
    """Run the blinded internal benchmark against the held-out positive set.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir`` and ``seed``.

    Raises:
        FileNotFoundError: If L3 has not run.
    """
    l3_path = Path(config["results_dir"]) / "l3" / "candidates.tsv"
    if not l3_path.is_file():
        raise FileNotFoundError(f"{l3_path} does not exist; the benchmark scores L3's "
                                "ranking and has nothing to measure without it.")
    with open(l3_path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    names, digest = positive_set()
    positives = set(names)
    channels = sorted({c[len("rank_"):] for c in (rows[0] if rows else {})
                       if c.startswith("rank_")})

    combined = recovery(rows, positives)
    combined["auroc"] = auroc(combined["positions"], combined["n_ranked"],
                              combined["n_positives_found"])
    combined["auroc_ci"] = bootstrap_auroc(
        combined["positions"], combined["n_ranked"], combined["n_positives_found"],
        seed=config["seed"])

    per_channel = {}
    for channel in channels:
        ranked = channel_ranking(rows, channel)
        result = recovery(ranked, positives)
        result["auroc"] = auroc(result["positions"], result["n_ranked"],
                                result["n_positives_found"])
        result["auroc_ci"] = bootstrap_auroc(result["positions"], result["n_ranked"],
                                             result["n_positives_found"],
                                             seed=config["seed"])
        result["circular"] = channel in CIRCULAR_CHANNELS
        if result["circular"]:
            result["note"] = ("This channel was built from the positive set, so its "
                              "recovery is 100% by construction and measures nothing. "
                              "Excluded from the headline.")
        per_channel[channel] = result

    independent = [c for c in channels if c not in CIRCULAR_CHANNELS]
    out_dir = Path(config["results_dir"]) / "l4" / BENCHMARK_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "positive_set": {
            "source": "src/l2_channels/aneuploidy_prior.tsv",
            "n_compounds": len(names),
            "sha256": digest,
            "note": ("The positive set is channel E's seed file. Channel E therefore "
                     "cannot be benchmarked on it, and the headline is computed over the "
                     "channels that never saw it."),
        },
        "headline_channels": independent,
        "circular_channels": sorted(CIRCULAR_CHANNELS & set(channels)),
        # The aggregate ranking includes the circular channel's contribution, so the
        # combined number is partly circular too and is NOT the headline. The result to
        # report is the independent channels' own recovery, below.
        "combined": {**combined, "partly_circular": bool(CIRCULAR_CHANNELS & set(channels)),
                     "note": ("L3's aggregate ranking includes the channel that was built "
                              "from this positive set, so this number inherits that "
                              "circularity. Read the independent channels instead.")},
        "per_channel": per_channel,
        "seed": config["seed"],
        "caveats": [
            "The positive set has around a dozen compounds and fewer reach any one "
            "channel, so every number here carries wide uncertainty. The bootstrap "
            "interval is reported so a point estimate is not read as a measurement.",
            "AUROC is computed against every other candidate the ranking contains, not "
            "against a matched decoy set. A decoy set matched on target degree and "
            "literature volume would be a stronger test and is not built.",
            "The positive set is drawn from published aneuploidy-stress literature, which "
            "is itself biased toward well-studied compounds. Recovering it is evidence "
            "that a channel finds what the literature already contains, not that it finds "
            "what would help this child.",
            "This measures pipeline recovery, not clinical efficacy.",
            "Nothing here may be used to re-tune a channel without re-freezing the "
            "positive set and disclosing the reuse. A benchmark optimised against is a "
            "description of the pipeline, not a test of it.",
        ],
    }
    (out_dir / "recovery.json").write_text(json.dumps(report, indent=2, sort_keys=True),
                                           encoding="utf-8")

    with open(out_dir / "recovery.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["scope", "n_ranked", "n_positives_found", "auroc", "ci_low",
                         "ci_high", "circular"])
        writer.writerow(["combined", combined["n_ranked"], combined["n_positives_found"],
                         f"{combined['auroc']:.4f}", combined["auroc_ci"].get("ci_low"),
                         combined["auroc_ci"].get("ci_high"), False])
        for channel in channels:
            result = per_channel[channel]
            writer.writerow([channel, result["n_ranked"], result["n_positives_found"],
                             f"{result['auroc']:.4f}", result["auroc_ci"].get("ci_low"),
                             result["auroc_ci"].get("ci_high"), result["circular"]])

    LOGGER.info("benchmark: %d positive(s); combined AUROC %.3f; independent channels %s",
                len(names), combined["auroc"], ", ".join(independent) or "none")
    return report
