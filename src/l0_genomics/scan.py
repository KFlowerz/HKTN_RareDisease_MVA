"""L0 -- VCF scanning: turn variant records into windowed allelic statistics.

Purpose
    Walk a single-sample VCF once, keeping only sites that can carry allelic signal, and
    accumulate per-window statistics for :mod:`src.l0_genomics.burden`.

Inputs
    A bgzipped, tabix-indexed single-sample VCF with ``FORMAT/AD`` and ``FORMAT/DP``.

Outputs
    :class:`~src.l0_genomics.burden.WindowStat` records, and the per-contig heterozygous
    fractions used to infer sex.

Guardrail
    Nothing leaves this module at variant resolution. Statistics are accumulated in
    windows and the per-site values are discarded as they are read -- no coordinate,
    genotype, or allele depth is retained, logged, or returned. In a ~50-patient
    population a single exact variant coordinate can be identifying.

    Sites are counted, never listed. If a future change needs per-site output, that is a
    decision to take deliberately against ``COMPLIANCE.md``, not a refactor.
"""

from __future__ import annotations

import logging
import statistics
from dataclasses import dataclass, field

from .burden import (
    DEFAULT_MIN_DP,
    PRIMARY_CONTIGS,
    WindowStat,
    baf,
)

LOGGER = logging.getLogger(__name__)

#: Window size for segmentation. Mosaic events are frequently segmental, so a
#: chromosome-wide mean dilutes a strong local signal into surrounding diploid sequence.
#: 10 Mb is a compromise: small enough to localize an arm-level event, large enough that
#: a window still clears MIN_INFORMATIVE_SITES at this dataset's het density.
DEFAULT_WINDOW = 10_000_000

#: Heterozygous genotype, normalized. Phase is irrelevant to an allele ratio, so
#: 0/1, 0|1 and 1|0 all reduce to this.
HET_NORMALIZED = (0, 1)


@dataclass
class _Acc:
    """Mutable accumulator for one window. Holds no per-site data."""

    n_sites: int = 0
    n_het: int = 0
    dev_sum: float = 0.0
    dps: list = field(default_factory=list)


def _is_biallelic_snv(rec) -> bool:
    """True for a single-nucleotide, single-ALT record."""
    if rec.alts is None or len(rec.alts) != 1:
        return False
    return len(rec.ref) == 1 and len(rec.alts[0]) == 1


def scan_vcf(
    vcf_path,
    *,
    min_dp: int = DEFAULT_MIN_DP,
    window: int = DEFAULT_WINDOW,
    contigs=PRIMARY_CONTIGS,
    pass_only: bool = True,
):
    """Accumulate windowed allelic statistics from a VCF.

    Args:
        vcf_path: Path to a bgzipped, indexed single-sample VCF.
        min_dp: Minimum ``FORMAT/DP`` for a site to count as informative.
        window: Window size in bases.
        contigs: Contigs to scan. Defaults to the primary assembly only.
        pass_only: Keep only records whose FILTER is PASS.

    Returns:
        ``(windows, het_fraction_by_contig)``.

    Raises:
        ValueError: If the VCF holds anything other than exactly one sample -- every
            statistic here assumes a single sample, and a multi-sample file would
            silently mix individuals.
    """
    import pysam

    # htslib writes warnings to stderr that quote the VCF's path, and this dataset's
    # filename embeds a lab accession and a sequencer flowcell id, which must never
    # propagate into a log or an artifact (see run.py's guardrail). Silence them.
    pysam.set_verbosity(0)

    vcf = pysam.VariantFile(str(vcf_path))
    if len(vcf.header.samples) != 1:
        raise ValueError(
            f"expected exactly 1 sample, found {len(vcf.header.samples)}: "
            "per-site allele ratios are meaningless across mixed samples"
        )

    present = set(vcf.header.contigs)
    wanted = [c for c in contigs if c in present]
    if not wanted:
        raise ValueError(
            f"none of the requested contigs are in the VCF header. Expected names like "
            f"{contigs[0]!r}; the header uses e.g. {sorted(present)[:3]!r}. "
            "A `chr` prefix mismatch is the usual cause."
        )
    if len(wanted) < len(contigs):
        LOGGER.warning("contigs absent from VCF header, skipped: %s", sorted(set(contigs) - present))

    acc: dict = {}
    per_contig = {c: [0, 0] for c in wanted}  # contig -> [informative, het]

    for contig in wanted:
        for rec in vcf.fetch(contig):
            if pass_only and rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue
            if not _is_biallelic_snv(rec):
                continue

            sample = rec.samples[0]
            ad = sample.get("AD")
            if not ad or len(ad) < 2 or ad[0] is None or ad[1] is None:
                continue
            ref_d, alt_d = int(ad[0]), int(ad[1])
            if ref_d + alt_d < min_dp:
                continue

            key = (contig, rec.pos // window)
            a = acc.setdefault(key, _Acc())
            a.n_sites += 1
            per_contig[contig][0] += 1

            gt = sample.get("GT")
            if gt is None or tuple(sorted(g for g in gt if g is not None)) != HET_NORMALIZED:
                continue

            value = baf(ref_d, alt_d)
            if value is None:
                continue
            a.n_het += 1
            a.dev_sum += abs(value - 0.5)
            a.dps.append(ref_d + alt_d)
            per_contig[contig][1] += 1

    windows = []
    for (contig, idx), a in sorted(acc.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        windows.append(
            WindowStat(
                contig=contig,
                start=idx * window,
                end=(idx + 1) * window,
                n_sites=a.n_sites,
                n_het=a.n_het,
                mean_abs_dev=(a.dev_sum / a.n_het) if a.n_het else None,
                median_dp=statistics.median(a.dps) if a.dps else None,
            )
        )

    het_fraction = {
        c: (h / n if n else 0.0) for c, (n, h) in per_contig.items()
    }
    return windows, het_fraction
