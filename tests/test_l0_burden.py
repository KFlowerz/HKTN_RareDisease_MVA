"""Tests for the L0 aneuploidy-burden analysis.

Every fixture here is synthetic and invented. No patient data is read, and no value from
the real dataset appears -- including as an expected result. A test that encoded the
subject's actual burden would put a clinical finding about a child into the repository.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from src.l0_genomics import burden
from src.l0_genomics.scan import scan_vcf


# --------------------------------------------------------------------------- maths


def test_baf_handles_zero_coverage() -> None:
    assert burden.baf(0, 0) is None
    assert burden.baf(10, 10) == pytest.approx(0.5)
    assert burden.baf(0, 20) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "fraction",
    [0.1, 0.25, 0.5, 0.75, 1.0],
)
def test_mosaic_fraction_gain_round_trips(fraction: float) -> None:
    """dev = f / (2(2+f)) must invert back to f."""
    dev = fraction / (2 * (2 + fraction))
    assert burden.mosaic_fraction(dev, "gain") == pytest.approx(fraction, abs=1e-9)


@pytest.mark.parametrize("fraction", [0.1, 0.25, 0.5, 0.75, 1.0])
def test_mosaic_fraction_loss_round_trips(fraction: float) -> None:
    """dev = f / (2(2-f)) must invert back to f."""
    dev = fraction / (2 * (2 - fraction))
    assert burden.mosaic_fraction(dev, "loss") == pytest.approx(fraction, abs=1e-9)


def test_pure_trisomy_sits_at_the_ceiling() -> None:
    """A fully clonal gain produces dev = 1/6; nothing above that is a pure gain."""
    assert burden.mosaic_fraction(1 / 6, "gain") == pytest.approx(1.0)
    assert burden.mosaic_fraction(0.30, "gain") is None


def test_gain_and_loss_disagree_for_the_same_deviation() -> None:
    """The models are not interchangeable -- which is why both are reported.

    Reporting a single fraction would hide a real ambiguity: BAF alone cannot tell a
    gain from a loss.
    """
    dev = 0.05
    assert burden.mosaic_fraction(dev, "gain") != pytest.approx(
        burden.mosaic_fraction(dev, "loss")
    )


def test_unknown_model_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown model"):
        burden.mosaic_fraction(0.05, "duplication")


# --------------------------------------------------------------------------- sex


def test_hemizygous_x_infers_xy() -> None:
    assert burden.infer_sex(0.06, 50_000) == "XY"


def test_diploid_x_infers_xx() -> None:
    assert burden.infer_sex(0.60, 50_000) == "XX"


def test_sex_call_declined_on_thin_evidence() -> None:
    """Sex drives which model is applied to X, so a thin call is declined, not guessed."""
    assert burden.infer_sex(0.06, 10) == "undetermined"


# --------------------------------------------------------------------------- baseline


def _win(contig: str, idx: int, dev: float, n_het: int = 5000) -> burden.WindowStat:
    return burden.WindowStat(contig, idx * 10, (idx + 1) * 10, n_het, n_het, dev, 30.0)


def test_baseline_is_robust_to_aneuploid_windows() -> None:
    """An aneuploid chromosome must not inflate the null it is scored against."""
    normal = [_win("1", i, 0.07) for i in range(20)]
    aneuploid = [_win("21", i, 0.15) for i in range(5)]
    mean, sd = burden.summarize(normal + aneuploid)
    assert mean == pytest.approx(0.07, abs=0.005)
    assert sd >= 0.0


def test_thin_windows_are_insufficient_not_zero() -> None:
    """A window with too few het sites says so; it never reports a burden of zero.

    Zero reads as 'no aneuploidy detected' when the truth is 'not measurable here'.
    """
    thin = burden.WindowStat("7", 0, 10, 5, 5, 0.0, 30.0)
    scored = burden.score_windows([thin], 0.07, 0.01)
    assert scored[0].verdict == "insufficient_sites"
    assert scored[0].z is None


def test_elevated_window_is_flagged() -> None:
    scored = burden.score_windows([_win("21", 0, 0.15)], 0.07, 0.01)
    assert scored[0].verdict == "imbalanced"
    assert scored[0].z == pytest.approx(8.0)


# --------------------------------------------------------------------------- screen


def test_globally_elevated_baseline_reads_as_contamination() -> None:
    """Contamination lifts every chromosome; mosaicism lifts some and leaves the rest."""
    assert burden.contamination_indicator(0.20, 0.01) == "suspect"
    assert burden.contamination_indicator(0.07, 0.01) == "clean"
    assert burden.contamination_indicator(0.0, 0.0) == "unknown"


# --------------------------------------------------------------------------- scanning

VCF_HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##contig=<ID=1,length=1000000>
##contig=<ID=X,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""


def _write_vcf(tmp_path: Path, rows: list) -> Path:
    """Write and index a tiny synthetic VCF. Invented data only."""
    pysam = pytest.importorskip("pysam")
    raw = tmp_path / "synthetic.vcf"
    raw.write_text(VCF_HEADER + "".join(rows), encoding="utf-8")
    gz = str(tmp_path / "synthetic.vcf.gz")
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)
    return Path(gz)


def _row(contig: str, pos: int, gt: str, ref_d: int, alt_d: int, flt: str = "PASS") -> str:
    return f"{contig}\t{pos}\t.\tA\tG\t100\t{flt}\t.\tGT:AD:DP\t{gt}:{ref_d},{alt_d}:{ref_d + alt_d}\n"


def test_scan_counts_het_sites_and_ignores_low_depth(tmp_path: Path) -> None:
    rows = [_row("1", 100 + i, "0/1", 25, 25) for i in range(10)]
    rows += [_row("1", 500 + i, "0/1", 2, 2) for i in range(10)]  # below min_dp
    windows, het = scan_vcf(_write_vcf(tmp_path, rows), min_dp=10, window=1_000_000)
    assert sum(w.n_het for w in windows) == 10
    assert het["1"] == pytest.approx(1.0)


def test_scan_skips_non_pass_and_non_snv(tmp_path: Path) -> None:
    rows = [_row("1", 100, "0/1", 25, 25, flt="LowQual")]
    rows += ["1\t200\t.\tAT\tG\t100\tPASS\t.\tGT:AD:DP\t0/1:25,25:50\n"]  # not an SNV
    rows += [_row("1", 300, "0/1", 25, 25)]
    windows, _ = scan_vcf(_write_vcf(tmp_path, rows), min_dp=10, window=1_000_000)
    assert sum(w.n_het for w in windows) == 1


def test_scan_measures_a_planted_imbalance(tmp_path: Path) -> None:
    """A window built at a known BAF must report the deviation that implies."""
    rows = [_row("1", 100 + i, "0/1", 20, 30) for i in range(500)]  # BAF 0.6, dev 0.1
    windows, _ = scan_vcf(_write_vcf(tmp_path, rows), min_dp=10, window=1_000_000)
    assert windows[0].mean_abs_dev == pytest.approx(0.1)
    assert burden.mosaic_fraction(0.1, "gain") == pytest.approx(0.5)


def test_scan_rejects_contig_name_mismatch(tmp_path: Path) -> None:
    """A `chr` prefix mismatch must fail loudly, not return an empty genome."""
    rows = [_row("1", 100, "0/1", 25, 25)]
    with pytest.raises(ValueError, match="chr` prefix mismatch|none of the requested contigs"):
        scan_vcf(_write_vcf(tmp_path, rows), contigs=("chr1", "chr2"))


def test_scan_rejects_multisample_vcf(tmp_path: Path) -> None:
    """Allele ratios are meaningless across mixed samples."""
    header = VCF_HEADER.replace("\tSAMPLE\n", "\tS1\tS2\n")
    raw = tmp_path / "multi.vcf"
    raw.write_text(header + "1\t100\t.\tA\tG\t100\tPASS\t.\tGT:AD:DP\t0/1:25,25:50\t0/1:25,25:50\n", encoding="utf-8")
    pysam = pytest.importorskip("pysam")
    gz = str(tmp_path / "multi.vcf.gz")
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)
    with pytest.raises(ValueError, match="exactly 1 sample"):
        scan_vcf(Path(gz))


# --------------------------------------------------------- thresholding discipline


def test_flat_genome_produces_no_calls() -> None:
    """A genome at uniform baseline must yield zero imbalanced windows.

    This is the regression that matters. The first implementation scored on z alone;
    because each window averages thousands of sites the baseline SD is ~0.002, so
    biologically trivial scatter cleared any z threshold and 20 of 22 autosomes were
    flagged. An effect-size floor is what separates significance from meaning here.
    """
    flat = [_win("1", i, 0.0700 + (i % 3) * 0.0005) for i in range(14)]
    mean, sd = burden.summarize(flat)
    scored = burden.score_windows(flat, mean, sd)
    assert not any(w.verdict == "imbalanced" for w in scored)


def test_effect_below_reportable_fraction_is_not_called() -> None:
    """A 5% mosaic is below what this method claims to resolve, so it is not claimed."""
    dev_5pct = 0.05 / (2 * (2 + 0.05))
    w = _win("1", 0, 0.07 + dev_5pct)
    assert burden.score_windows([w], 0.07, 0.0005)[0].verdict != "imbalanced"


def test_real_effect_is_still_detected() -> None:
    """The floor must not blind the method to an event it should find."""
    dev_30pct = 0.30 / (2 * (2 + 0.30))
    w = _win("1", 0, 0.07 + dev_30pct)
    assert burden.score_windows([w], 0.07, 0.0005)[0].verdict == "imbalanced"


def test_contig_rollup_uses_median_not_max() -> None:
    """One spike in an otherwise flat chromosome is baseline, never a call.

    Taking the max over ~14 windows per chromosome across 22 chromosomes is an
    uncorrected multiple-comparisons procedure and flags nearly everything. The lone
    window stays visible in ``n_flagged_windows`` without being promoted to a finding.
    """
    dev_50pct = 0.50 / (2 * (2 + 0.50))
    ws = [_win("5", i, 0.07) for i in range(13)] + [_win("5", 13, 0.07 + dev_50pct)]
    scored = burden.score_windows(ws, 0.07, 0.0005)
    roll = burden.rollup_contig(scored, 0.07)
    assert roll["verdict"] == "baseline"
    assert roll["longest_contiguous_flagged"] == 1
    assert roll["n_flagged_windows"] == 1


def test_contiguous_windows_make_a_call() -> None:
    """A genuine segmental event spans neighbouring windows and is called."""
    dev_40pct = 0.40 / (2 * (2 + 0.40))
    ws = [_win("21", i, 0.07 + dev_40pct) for i in range(8)]
    scored = burden.score_windows(ws, 0.07, 0.0005)
    roll = burden.rollup_contig(scored, 0.07)
    assert roll["verdict"] == "imbalanced"
    assert roll["longest_contiguous_flagged"] >= burden.MIN_SUPPORTING_WINDOWS
    assert roll["mosaic_fraction_if_gain"] == pytest.approx(0.40, abs=0.02)
