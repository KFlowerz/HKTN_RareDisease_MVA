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

import math
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
    median_depth: float | None = None
    sigma: float | None = None
    #: Smallest mosaic fraction this sample can resolve, measured from its own depth and
    #: scatter. Travels with the result so the report states what the method could have
    #: seen, rather than leaving a null finding to be read as "nothing is there".
    detection_limit: float | None = None


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


#: z above the sample's own window-to-window scatter required to call a window.
Z_THRESHOLD = 5.0

#: Contiguous informative windows a chromosome needs before it is called imbalanced.
#: A real segmental event spans neighbouring windows; a lone spike is far more often a
#: mapping artifact, a region of homozygosity, or a centromere-adjacent oddity.
MIN_SUPPORTING_WINDOWS = 2


def sigma_at_depth(depth: float) -> float:
    """Per-site standard deviation of BAF at a given read depth.

    For ``X ~ Binomial(d, 1/2)`` and ``BAF = X/d``, ``sd(BAF) = 0.5/sqrt(d)``.
    """
    return 0.5 / math.sqrt(depth) if depth and depth > 0 else float("nan")


def folded_mean(delta: float, sigma: float) -> float:
    """Expected ``|BAF - 0.5|`` when the true shift is ``delta`` and noise is ``sigma``.

    For ``X ~ Normal(0.5 ± delta, sigma)``::

        E|X - 0.5| = sigma·sqrt(2/pi)·exp(-delta²/2sigma²) + delta·erf(delta/(sigma·sqrt2))

    Note:
        This is the relationship the naive inversion ignored. The measured statistic is
        **not** the true shift: at ``delta = 0`` it already equals ``sigma·sqrt(2/pi)``
        -- the noise floor -- and it responds *quadratically* to small shifts, because
        sampling noise dominates them. Treating the excess over baseline as the shift
        therefore under-reports the mosaic fraction badly: at a true fraction of 0.40 it
        gives 0.12, and at 0.10 it gives 0.
    """
    if sigma <= 0:
        return abs(delta)
    z = delta / (sigma * math.sqrt(2.0))
    return sigma * math.sqrt(2.0 / math.pi) * math.exp(-(z * z)) + delta * math.erf(z)


def deconvolve_shift(statistic: float, sigma: float, *, tol: float = 1e-9) -> float | None:
    """Recover the true BAF shift from an observed folded mean.

    Args:
        statistic: Observed mean ``|BAF - 0.5|`` over a window.
        sigma: Per-site BAF standard deviation, from :func:`sigma_at_depth`.

    Returns:
        The shift ``delta`` in [0, 0.5), or None when the statistic lies below the noise
        floor (no shift can explain a value *smaller* than pure noise) or above 0.5.

    Note:
        Solved by bisection rather than SciPy's root finders, to keep this module free of
        heavy dependencies -- it is the numeric core and is exercised in every test run.
    """
    if sigma <= 0 or not math.isfinite(statistic):
        return None
    floor = folded_mean(0.0, sigma)
    if statistic <= floor:
        return 0.0
    hi = 0.4999
    if statistic >= folded_mean(hi, sigma):
        return None
    lo = 0.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if folded_mean(mid, sigma) < statistic:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2


def excess_shift(statistic: float, baseline: float, sigma: float) -> float | None:
    """Shift attributable to a region, over and above the sample's own floor.

    Args:
        statistic: The region's observed mean ``|BAF - 0.5|``.
        baseline: The sample's diploid baseline statistic.
        sigma: Per-site BAF standard deviation at this sample's depth.

    Returns:
        ``delta`` in excess of baseline, floored at 0, or None if undeterminable.

    Note:
        Deconvolving the raw statistic from zero is wrong, and wrong in a way that looks
        plausible: the baseline itself sits *above* the pure-binomial floor -- 0.0671
        observed against 0.0602 predicted in this dataset -- because of reference bias
        and mapping behaviour. Inverting the raw value therefore assigns that systematic
        excess to mosaicism and reports a fraction near 0.17 for **every** chromosome,
        including ones that are plainly at baseline.

        Both quantities are converted to shift space first, then subtracted. Subtracting
        in statistic space instead would be the original error over again, since the
        statistic is not linear in the shift.
    """
    if sigma is None or sigma <= 0:
        return None
    d_obs = deconvolve_shift(statistic, sigma)
    d_base = deconvolve_shift(baseline, sigma)
    if d_obs is None or d_base is None:
        return None
    return max(0.0, d_obs - d_base)


def detection_limit(baseline: float, sd: float, sigma: float, *, z: float = Z_THRESHOLD) -> float | None:
    """Smallest mosaic fraction (gain model) this sample can actually resolve.

    Args:
        baseline: The sample's diploid baseline statistic.
        sd: Window-to-window scatter of that statistic.
        sigma: Per-site BAF standard deviation at this sample's depth.
        z: Detection threshold in units of ``sd``.

    Returns:
        The fraction at which a window would just reach ``z``, or None if undeterminable.

    Note:
        Derived from **this sample's** depth and scatter rather than configured as a
        constant. A hardcoded sensitivity claim is a claim about a dataset nobody
        measured; this one is measured, and it belongs in the output so the report can
        state what the method could and could not have seen.

        **Model-dependent, and known to be so.** This assumes the sample's systematic
        baseline excess and a mosaic shift combine additively *in shift space*. That is
        an approximation: reference bias perturbs individual sites, while a mosaic gain
        shifts every heterozygous site in a region, and the two are not the same physical
        process. Computed without the additivity assumption the limit comes out roughly
        2.5x higher, and a spike-in simulation at this sample's depth put the z = 5
        crossing near 0.25. Treat the returned value as the optimistic end of a range
        whose pessimistic end is ~2.5x larger, and calibrate against simulated spike-ins
        before quoting a single number in the report.
    """
    if sd <= 0 or sigma <= 0:
        return None
    delta = excess_shift(baseline + z * sd, baseline, sigma)
    if delta is None:
        return None
    return mosaic_fraction(delta, "gain")


def score_windows(
    windows: list,
    baseline_mean: float,
    baseline_sd: float,
    *,
    z_threshold: float = Z_THRESHOLD,
) -> list:
    """Attach a z-score and verdict to each window, against the sample's own baseline.

    Note:
        A single criterion is correct here now, where it was not before. The earlier
        version needed a separate effect-size floor because the statistic was being read
        as if it were the true shift, so statistical significance and biological meaning
        had come apart. With the noise model in place, ``z >= z_threshold`` corresponds
        to a specific, computable mosaic fraction -- :func:`detection_limit` reports
        which -- so the threshold *is* the biological floor rather than standing next to
        one.

        The multiple-comparisons discipline still lives in :func:`rollup_contig`: window
        scoring is deliberately permissive, and contiguity plus a median statistic decide
        what becomes a chromosome-level finding.
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

        z = ((w.mean_abs_dev - baseline_mean) / baseline_sd) if baseline_sd > 0 else None
        verdict = "unscored" if z is None else ("imbalanced" if z >= z_threshold else "baseline")

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


def rollup_contig(
    windows: list,
    baseline_mean: float,
    sigma: float | None = None,
    limit: float | None = None,
) -> dict:
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
    # Deconvolve the noise floor out of the median statistic before converting to a cell
    # fraction. Reading `excess` as the shift under-reports by 2-3x; see `folded_mean`.
    if sigma is None:
        dps = [w.median_dp for w in informative if w.median_dp]
        sigma = sigma_at_depth(statistics.median(dps)) if dps else None

    delta = excess_shift(mid, baseline_mean, sigma) if sigma else None
    f_gain = mosaic_fraction(delta, "gain") if delta is not None else None
    f_loss = mosaic_fraction(delta, "loss") if delta is not None else None

    # A chromosome is only called when its own estimated fraction clears the sensitivity
    # the method actually has. Without this, verdicts and the reported detection limit
    # contradict each other -- calling a chromosome imbalanced at f = 0.19 while stating
    # the method cannot resolve below f = 0.27 is not a finding anyone can act on.
    if run < MIN_SUPPORTING_WINDOWS:
        verdict = "baseline"
    elif limit is not None and (f_gain is None or f_gain < limit):
        verdict = "below_detection_limit"
    else:
        verdict = "imbalanced"

    return {
        "verdict": verdict,
        "n_windows": len(windows),
        "n_informative_windows": len(informative),
        "n_flagged_windows": n_flagged,
        "longest_contiguous_flagged": run,
        "median_abs_dev": mid,
        "excess_over_baseline": excess,
        "deconvolved_shift": delta,
        "mosaic_fraction_if_gain": f_gain,
        "mosaic_fraction_if_loss": f_loss,
    }
