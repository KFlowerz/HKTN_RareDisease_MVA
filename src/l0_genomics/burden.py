"""L0 -- aneuploidy burden from B-allele frequency.

Purpose
    Quantify per-chromosome and per-window allelic imbalance from a single-sample VCF,
    and convert it into a mosaic-fraction estimate. MVA's phenotype is mosaic aneuploidy,
    so this is a feature of the disease, not background context.

Inputs
    A bgzipped, tabix-indexed single-sample VCF carrying ``FORMAT/AD`` and ``FORMAT/DP``.

Outputs
    Per-window and per-chromosome records: informative-site counts, mean |BAF - 0.5|,
    a z-score against the sample's own autosomal baseline, and mosaic-fraction estimates
    under both the gain and loss models.

Method
    In a disomic region, heterozygous sites cluster near an allele ratio of 0.5. A mosaic
    gain or loss splits that cluster, and the magnitude of the split scales with the
    fraction of cells carrying the alteration, so one statistic both detects the event
    and estimates its cell fraction [conlin2010] doi:10.1093/hmg/ddq003. The phase-aware
    form of the same signal detects mosaic chromosomal alterations at low cell fractions
    [loh2018] doi:10.1038/s41586-018-0321-x.

    An allele ratio is measured between two alleles at one locus, so read depth, GC
    content, and mappability affect numerator and denominator alike and cancel. None of
    the corrections a depth-based caller requires arise here.

Guardrail
    Reads only from the VCF path given; writes nothing. Emits no variant coordinates --
    windows are genomic intervals, and the values returned are aggregate statistics over
    many sites. This module must stay free of anything that could identify the subject.

    **A deviation is not a diagnosis.** Low-level sample contamination and low-fraction
    mosaicism produce the same BAF signature, so :func:`contamination_indicator` must be
    consulted before any mosaic claim, and its verdict travels with the result.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field

#: Primary assembly contigs, no `chr` prefix -- matching this dataset's VCF, which
#: declares 2,580 contigs including hs38d1 decoys and alts. Restricting here keeps
#: decoy noise out of the burden vector entirely.
AUTOSOMES = tuple(str(i) for i in range(1, 23))
SEX_CONTIGS = ("X", "Y")
PRIMARY_CONTIGS = AUTOSOMES + SEX_CONTIGS

#: Below this many informative heterozygous sites a window says `insufficient_sites`.
#: Never zero: zero reads as "no aneuploidy detected" when the truth is "not measurable".
MIN_INFORMATIVE_SITES = 200

#: Minimum FORMAT/DP for a site to be informative. Low-depth sites have BAF variance that
#: swamps the mosaic signal being measured.
DEFAULT_MIN_DP = 10

#: Heterozygous fraction below which a contig is treated as hemizygous. A diploid autosome
#: in this dataset runs ~60%; a single X runs ~6%.
HEMIZYGOUS_HET_FRACTION = 0.20

#: Largest |BAF - 0.5| a fully clonal event can produce under each model, at f = 1.
#: Gain: dev = f / (2(2+f)) -> 1/6. Loss: dev = f / (2(2-f)) -> 1/2.
#: A deviation above the ceiling is not a more extreme mosaic -- it is evidence the model
#: does not fit, and the caller is told so rather than handed a saturated number.
MAX_DEV_GAIN = 1.0 / 6.0
MAX_DEV_LOSS = 0.5


@dataclass(frozen=True)
class WindowStat:
    """Aggregate allelic statistics for one genomic window."""

    contig: str
    start: int
    end: int
    n_sites: int
    n_het: int
    mean_abs_dev: float | None
    median_dp: float | None
    z: float | None = None
    verdict: str = "pending"

    @property
    def informative(self) -> bool:
        return self.n_het >= MIN_INFORMATIVE_SITES


@dataclass
class BurdenResult:
    """The full per-sample burden output."""

    baseline_mean_dev: float
    baseline_sd: float
    inferred_sex: str
    windows: list = field(default_factory=list)
    contigs: dict = field(default_factory=dict)
    contamination: str = "unknown"
    notes: list = field(default_factory=list)


def baf(ref_depth: int, alt_depth: int) -> float | None:
    """B-allele frequency at one site, or None when there is no coverage."""
    total = ref_depth + alt_depth
    return None if total <= 0 else alt_depth / total


def mosaic_fraction(abs_dev: float, model: str) -> float | None:
    """Convert a mean |BAF - 0.5| into an estimated fraction of altered cells.

    Args:
        abs_dev: Mean absolute deviation of heterozygous BAF from 0.5.
        model: ``"gain"`` (one extra copy in the altered clone) or ``"loss"``
            (one copy lost).

    Returns:
        Estimated cell fraction in [0, 1], or None if the deviation is out of range
        for the model.

    Note:
        For a clone at fraction ``f`` carrying an extra copy, the gained allele's
        expected ratio is ``(1 + f) / (2 + f)``, giving ``dev = f / (2(2 + f))`` and so
        ``f = 4·dev / (1 - 2·dev)``. For a clonal loss the ratio is ``(1 - f) / (2 - f)``,
        giving ``dev = f / (2(2 - f))`` and ``f = 4·dev / (1 + 2·dev)``.

        **The two models are not distinguishable from BAF alone** -- the same deviation
        implies a different cell fraction depending on whether copies were gained or
        lost, which is why both are reported and depth is consulted as a tie-breaker.
        Reporting one number here would hide a real ambiguity behind a false precision.
        At ``f = 1`` the gain model gives ``dev = 1/6 ≈ 0.167``; deviations near that
        ceiling indicate a near-clonal event rather than a mosaic one.
    """
    if model not in ("gain", "loss"):
        raise ValueError(f"unknown model {model!r}; expected 'gain' or 'loss'")
    if abs_dev <= 0:
        return 0.0

    # Each model has a hard ceiling at f = 1, and a deviation above it cannot have been
    # produced by that model at any cell fraction. Clamping to 1.0 instead of returning
    # None would report a fully clonal event where the real answer is "these data are
    # not explained by a mosaic gain" -- most often reference bias, mismapping, or
    # contamination. Saying so is the useful output.
    ceiling = MAX_DEV_GAIN if model == "gain" else MAX_DEV_LOSS
    if abs_dev > ceiling:
        return None

    if model == "gain":
        return min(1.0, 4 * abs_dev / (1 - 2 * abs_dev))
    return min(1.0, 4 * abs_dev / (1 + 2 * abs_dev))


def infer_sex(het_fraction_x: float, n_sites_x: int) -> str:
    """Infer sex from X heterozygosity, or decline to.

    Args:
        het_fraction_x: Heterozygous fraction among informative X sites.
        n_sites_x: Number of informative X sites the fraction is based on.

    Returns:
        ``"XY"``, ``"XX"``, or ``"undetermined"``.

    Note:
        A hemizygous X yields almost no heterozygous calls; a diploid X behaves like an
        autosome. Below :data:`MIN_INFORMATIVE_SITES` the call is declined rather than
        guessed -- an inferred sex drives which model is applied to X, so a wrong call
        propagates silently.
    """
    if n_sites_x < MIN_INFORMATIVE_SITES:
        return "undetermined"
    return "XY" if het_fraction_x < HEMIZYGOUS_HET_FRACTION else "XX"


def contamination_indicator(baseline_mean_dev: float, baseline_sd: float) -> str:
    """Judge whether the genome-wide baseline looks contaminated.

    Returns:
        ``"clean"``, ``"suspect"``, or ``"unknown"``.

    Note:
        Contamination raises |BAF - 0.5| across **every** chromosome at once, because
        foreign reads perturb allele ratios genome-wide. Mosaic aneuploidy raises it on
        **some** and leaves the rest at baseline. So a globally elevated baseline is the
        signature to worry about, and it must be ruled out before anything is called
        mosaic: calling contamination "mosaicism" would misstate the central feature of
        this child's disease, and calling mosaicism "contamination" would discard the
        finding the pipeline exists to make.

        This is a coarse screen, not a contamination estimator. A dedicated tool
        (VerifyBamID-style) needs alignments, and this dataset ships none.
    """
    if baseline_mean_dev <= 0:
        return "unknown"
    if baseline_mean_dev > 0.12:
        return "suspect"
    return "clean"


def summarize(windows: list) -> tuple:
    """Compute the sample's own diploid baseline from its autosomal windows.

    Args:
        windows: All :class:`WindowStat` records, any contig.

    Returns:
        ``(mean, sd)`` of ``mean_abs_dev`` over informative autosomal windows.

    Note:
        The baseline is derived from the sample rather than a literature constant, so it
        absorbs this sample's coverage, contamination, and reference-bias behaviour. It
        is computed from the **median** of autosomal windows rather than the mean, so a
        genuinely aneuploid chromosome does not inflate the null it is scored against.
    """
    devs = [
        w.mean_abs_dev
        for w in windows
        if w.contig in AUTOSOMES and w.informative and w.mean_abs_dev is not None
    ]
    if len(devs) < 2:
        return (0.0, 0.0)
    centre = statistics.median(devs)
    # Median absolute deviation, scaled to a normal-consistent sd. Robust to the
    # aneuploid windows we are trying to detect.
    mad = statistics.median([abs(d - centre) for d in devs])
    return (centre, mad * 1.4826)


#: Smallest mosaic fraction this method claims to resolve. Below it, BAF deviation is not
#: separable from reference bias and mapping noise in a VCF without alignments.
MIN_REPORTABLE_FRACTION = 0.10

#: The absolute BAF deviation MIN_REPORTABLE_FRACTION implies under the gain model:
#: dev = f / (2(2+f)). An effect-size floor is required alongside any z-score, because
#: each window averages thousands of sites, so the window-mean SD is ~0.002 and a
#: biologically trivial difference clears any z threshold. Statistical significance and
#: biological meaning come apart badly here, and only the effect size tracks the latter.
MIN_EFFECT = MIN_REPORTABLE_FRACTION / (2 * (2 + MIN_REPORTABLE_FRACTION))

#: Contiguous informative windows a chromosome needs before it is called imbalanced.
#: A real segmental event spans neighbouring windows; a lone spike is far more often a
#: mapping artifact, a region of homozygosity, or a centromere-adjacent oddity.
MIN_SUPPORTING_WINDOWS = 2


def score_windows(
    windows: list,
    baseline_mean: float,
    baseline_sd: float,
    *,
    min_effect: float = MIN_EFFECT,
) -> list:
    """Attach a z-score and verdict to each window, against the sample's own baseline.

    A window is ``imbalanced`` only if it clears **both** the statistical test and the
    effect-size floor. Requiring both is the point: with thousands of sites per window
    the z-score alone flags differences far too small to be a mosaic event, and an
    effect size alone ignores how well-determined the window is.
    """
    scored = []
    for w in windows:
        if not w.informative or w.mean_abs_dev is None:
            scored.append(
                WindowStat(
                    w.contig, w.start, w.end, w.n_sites, w.n_het,
                    w.mean_abs_dev, w.median_dp, None, "insufficient_sites",
                )
            )
            continue

        effect = w.mean_abs_dev - baseline_mean
        z = (effect / baseline_sd) if baseline_sd > 0 else None

        if z is None:
            verdict = "unscored"
        elif effect >= min_effect and z >= 5:
            verdict = "imbalanced"
        elif effect >= min_effect / 2 and z >= 5:
            verdict = "equivocal"
        else:
            verdict = "baseline"

        scored.append(
            WindowStat(
                w.contig, w.start, w.end, w.n_sites, w.n_het,
                w.mean_abs_dev, w.median_dp, z, verdict,
            )
        )
    return scored


def _longest_run(flags: list) -> int:
    """Length of the longest run of True in a list."""
    best = run = 0
    for f in flags:
        run = run + 1 if f else 0
        best = max(best, run)
    return best


def rollup_contig(windows: list, baseline_mean: float) -> dict:
    """Summarize one contig from its scored windows.

    Args:
        windows: Scored :class:`WindowStat` records for a single contig, in position order.
        baseline_mean: The sample's diploid baseline.

    Returns:
        A verdict dict for the contig.

    Note:
        The representative statistic is the **median** of informative windows, not the
        maximum. Taking the max over ~14 windows per chromosome, across 22 chromosomes,
        is an uncorrected multiple-comparisons procedure: it will flag nearly every
        chromosome regardless of truth. The median asks the different and correct
        question -- is this chromosome as a whole shifted?

        A chromosome is called imbalanced only when at least
        :data:`MIN_SUPPORTING_WINDOWS` *contiguous* windows are individually imbalanced,
        because a real segmental event is contiguous and an isolated spike usually is not.
    """
    informative = [w for w in windows if w.informative and w.mean_abs_dev is not None]
    if not informative:
        return {"verdict": "insufficient_sites", "n_windows": len(windows)}

    devs = sorted(w.mean_abs_dev for w in informative)
    mid = devs[len(devs) // 2] if len(devs) % 2 else (devs[len(devs) // 2 - 1] + devs[len(devs) // 2]) / 2
    excess = max(0.0, mid - baseline_mean)

    flags = [w.verdict == "imbalanced" for w in informative]
    run = _longest_run(flags)
    n_flagged = sum(flags)

    # Contiguity gates BOTH calls. An isolated flagged window among twenty is noise --
    # mapping artifacts, regions of homozygosity, centromere-adjacent oddities -- and
    # letting one flag make a chromosome "equivocal" restates the multiple-comparisons
    # error one notch quieter. The flagged count is still reported, so an isolated
    # window is visible without being promoted to a finding.
    if run < MIN_SUPPORTING_WINDOWS:
        verdict = "baseline"
    elif excess >= MIN_EFFECT:
        verdict = "imbalanced"
    else:
        verdict = "segmental_only"

    return {
        "verdict": verdict,
        "n_windows": len(windows),
        "n_informative_windows": len(informative),
        "n_flagged_windows": n_flagged,
        "longest_contiguous_flagged": run,
        "median_abs_dev": mid,
        "excess_over_baseline": excess,
        "mosaic_fraction_if_gain": mosaic_fraction(excess, "gain"),
        "mosaic_fraction_if_loss": mosaic_fraction(excess, "loss"),
    }
