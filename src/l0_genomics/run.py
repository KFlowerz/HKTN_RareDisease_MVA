"""L0 -- Genomics.

Purpose
    Turn the subject's whole-genome sequencing data into the scientific premise of the
    run: which SAC gene is biallelically disrupted, what the variant does to the
    protein, and which pathway module is affected. Additionally quantify **aneuploidy
    burden** per chromosome -- MVA's phenotype literally is mosaic aneuploidy, so the
    degree of missegregation is itself a feature, not just context.

Inputs
    A single-sample WGS **VCF** (bgzipped, with tabix index), read from
    ``config["data_dir"]``.

    **There is no BAM.** The gated dataset ships raw reads (``.fastq.gz``) and a called
    VCF, not alignments -- verified against dataset revision
    ``59e322d27f399006b398d366d33e703e48a29914`` on 2026-09-02; see ``DATA.md``. Any
    depth-based method requiring alignments would first have to align ~85 GB of FASTQ,
    which is out of scope. This layer is therefore **VCF-only**.

    **There is no Track-1 answer to reconcile against.** The dataset's phenotype document
    names no causal gene (verified 2026-09-07), so the reconciliation step this layer was
    originally specified around has no counterpart. L0 must make the call independently
    and report it as a *finding with its supporting evidence*, not as agreement with an
    external source. If a Track-1 result is published later, reconcile then and fail
    loudly on disagreement -- but do not block on a source that does not exist.

VCF facts this layer must accommodate (verified 2026-09-07)
    - **GRCh38**, ``GCA_000001405.15_..._plus_hs38d1``, GATK ``VariantFiltration``-hardened.
    - **Contig names carry no ``chr`` prefix** (``1``, ``2``, ... ``X``, ``Y``). Any
      annotation resource that uses ``chr`` must be translated, in one documented place.
    - **2,580 contigs**, including decoys and alts from ``hs38d1``. Restrict every
      per-chromosome statistic to the primary set or the burden vector fills with noise
      from contigs that carry no meaningful allelic signal.

Outputs
    Written under ``config["results_dir"]/l0_genomics/``:
      - ``aneuploidy_burden.json`` -- the sample's own diploid baseline, inferred sex, a
        contamination screen, and per-contig verdicts with mosaic-fraction estimates
        under both the gain and loss models
      - ``windows.json`` -- per-window statistics with z-scores against that baseline
      - ``variant_calls.json`` -- every non-reference allele the sample carries in the
        panel regions, with its genotype record and its decision-D5
        :class:`~src.l0_genomics.transcripts.VariantCall` (consequence, transcript, LoF
        tier)
      - ``causal_gene_call.json`` -- per-gene biallelic configurations, the overall
        verdict and candidate gene(s) for gate G1, provenance, and caveats

    Burden is written **before** the causal-gene half runs, so a failure there does not
    discard a completed scan of a 5-million-record VCF. The affected pathway module is
    derived from the causal gene in L1, whose purpose that is.

Method note -- aneuploidy burden without alignments
    Burden is derived from **B-allele frequency (BAF) at heterozygous sites**, using
    ``FORMAT/AD`` (per-allele read depths) in the VCF. In a disomic region heterozygous
    sites cluster near an allele ratio of 0.5; a mosaic gain or loss splits that cluster,
    and the *magnitude* of the split scales with the fraction of cells carrying the
    alteration -- so the same statistic both detects the event and estimates its mosaic
    fraction [conlin2010] doi:10.1093/hmg/ddq003. Phase-aware treatment of the same signal
    detects mosaic chromosomal alterations down to low cell fractions
    [loh2018] doi:10.1038/s41586-018-0321-x.

    This is preferable here on its own merits, not merely as a substitute. An allele ratio
    is measured between two alleles at a single locus, so coverage depth, GC content, and
    mappability affect numerator and denominator alike and cancel -- the corrections a
    depth-based caller needs do not arise. Median ``FORMAT/DP`` per chromosome, normalized
    to the autosomal median, is retained only as a weaker secondary signal; it is
    confounded by variant density and capture behaviour and must never be reported alone.

    **Calibrate the null from the sample itself.** The autosomes supply tens of thousands
    of informative heterozygous sites each at ``DP >= 10``, so the per-sample distribution
    of |BAF - 0.5| over the bulk of the genome *is* the diploid baseline. Score each
    chromosome against that internal baseline rather than against a literature constant:
    it absorbs this sample's coverage, contamination, and reference-bias behaviour, which
    an external threshold cannot.

    **Sex chromosomes need their own model, and Y needs masking.** A hemizygous X produces
    almost no heterozygous calls, so the autosomal BAF model does not apply to it and must
    not be run there by default. Y is worse: apparent heterozygosity on a hemizygous contig
    is mismapping, concentrated in pseudoautosomal and repetitive regions, and an unmasked
    Y will report the largest apparent burden in the genome -- a pure artifact presented as
    the headline finding. Mask PAR and low-complexity regions, infer sex from X
    heterozygosity, and handle X and Y explicitly or exclude them with the reason recorded.

    **Segment; do not average whole chromosomes.** Mosaic events are frequently segmental,
    and a whole-chromosome mean dilutes a strong local signal into the surrounding diploid
    genome. Segment along each chromosome and report the segment, not just the chromosome.

Method note -- the causal gene
    Panel alleles are annotated **locally** by snpEff in the two passes decision D5
    prescribes (MANE Select, then all transcripts), classified LoF or not, and assembled
    into biallelic configurations from zygosity and GATK physical phasing. The result is a
    *candidate* for gate G1 shipped with its evidence, never a setting. See ``causal.py``
    for the method and ``annotate.py`` for why a remote annotation API is forbidden.

Guardrail
    Reads **only** from ``config["data_dir"]`` and writes **only** under
    ``config["results_dir"]`` -- both outside the repository and registered in
    ``docs/data_custody.md``. No patient-derived bytes (variant coordinates, allele
    depths, sample identifiers) may be written anywhere else, logged at INFO, or embedded
    in a report figure at a resolution that permits re-identification. This is a real
    minor's genome, and the source filenames themselves embed a lab accession and a
    sequencer flowcell identifier -- never propagate them into an artifact.

    The causal gene is a *finding*, never an assumption. With no Track-1 answer available,
    it must ship with the evidence that produced it -- the variants, their consequences,
    their zygosity, and what was considered and rejected -- so a reader can disagree with
    the call rather than take it on trust.

    **Phenotype is patient data and is read at runtime only.** The dataset's phenotype
    document is HPO-coded; parse it from ``config["data_dir"]`` when needed. Never write
    HPO term sets into config, source, a test fixture, or any committed file: a specific
    combination of features is identifying in a population of roughly 50 patients, and the
    family's own publications set the boundary for what is public about them
    (``COMPLIANCE.md``).

    BAF is uninformative where a chromosome has too few heterozygous sites to populate the
    distribution. Emit an explicit ``insufficient_sites`` verdict for such a chromosome --
    never a burden of zero, which would read as "no aneuploidy detected" when the truth is
    "not measurable here".
"""

from __future__ import annotations

import json
import logging
import statistics
from dataclasses import asdict
from pathlib import Path

from . import burden as _burden
from . import causal as _causal
from . import transcripts as _transcripts
from .scan import DEFAULT_WINDOW, scan_vcf

LOGGER = logging.getLogger(__name__)

#: Spindle-assembly-checkpoint panel. Symbols only -- coordinates are resolved from a
#: local annotation source at analysis time. Symbols are public identifiers and may be
#: sent to a public service to resolve; the subject's variants may not.
SAC_PANEL = ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192")


def _find_vcf(data_dir: Path) -> Path:
    """Locate the single-sample VCF in the data directory."""
    candidates = sorted(p for p in data_dir.rglob("*.vcf.gz") if not p.name.endswith(".tbi"))
    if not candidates:
        raise FileNotFoundError(f"no .vcf.gz found under {data_dir}")
    if len(candidates) > 1:
        raise ValueError(
            f"expected exactly one VCF under {data_dir}, found {len(candidates)}. "
            "Name the intended one in config rather than guessing."
        )
    return candidates[0]


def compute_burden(vcf_path: Path, *, min_dp: int, window: int) -> _burden.BurdenResult:
    """Run the aneuploidy-burden analysis over a VCF.

    Args:
        vcf_path: The single-sample VCF.
        min_dp: Minimum depth for an informative site.
        window: Segmentation window size in bases.

    Returns:
        A :class:`~src.l0_genomics.burden.BurdenResult`.
    """
    windows, het_fraction = scan_vcf(vcf_path, min_dp=min_dp, window=window)

    baseline_mean, baseline_sd = _burden.summarize(windows)
    scored = _burden.score_windows(windows, baseline_mean, baseline_sd)
    contamination = _burden.contamination_indicator(baseline_mean, baseline_sd)

    # Per-site BAF noise at this sample's own depth. Everything downstream -- the
    # fraction estimates and the sensitivity claim -- is calibrated against it.
    depths = [w.median_dp for w in scored if w.median_dp]
    median_depth = statistics.median(depths) if depths else None
    sigma = _burden.sigma_at_depth(median_depth) if median_depth else None
    limit = _burden.detection_limit(baseline_mean, baseline_sd, sigma) if sigma else None

    x_sites = sum(w.n_sites for w in scored if w.contig == "X")
    sex = _burden.infer_sex(het_fraction.get("X", 0.0), x_sites)

    notes = []
    if contamination == "suspect":
        notes.append(
            "Genome-wide baseline deviation is elevated on every chromosome, which is the "
            "signature of sample contamination rather than mosaicism. No window may be "
            "called mosaic until this is resolved."
        )

    # Per-contig rollup. X is scored only when diploid; Y never is -- apparent
    # heterozygosity on a hemizygous contig is mismapping in PAR and repetitive regions,
    # and an unmasked Y reports the largest apparent burden in the genome.
    contigs: dict = {}
    for contig in sorted({w.contig for w in scored}):
        ws = [w for w in scored if w.contig == contig]
        informative = [w for w in ws if w.informative and w.mean_abs_dev is not None]

        if contig == "Y":
            contigs[contig] = {"verdict": "excluded", "reason": "hemizygous_mismapping"}
            continue
        if contig == "X" and sex != "XX":
            contigs[contig] = {
                "verdict": "excluded",
                "reason": f"not_diploid_under_inferred_sex_{sex}",
            }
            continue
        contigs[contig] = _burden.rollup_contig(
            sorted(ws, key=lambda w: w.start), baseline_mean, sigma, limit
        )

    result = _burden.BurdenResult(
        baseline_mean_dev=baseline_mean,
        baseline_sd=baseline_sd,
        inferred_sex=sex,
        windows=scored,
        contigs=contigs,
        contamination=contamination,
        notes=notes,
    )
    result.median_depth = median_depth
    result.sigma = sigma
    result.detection_limit = limit
    return result


def run(config: dict) -> None:
    """Execute layer L0.

    Args:
        config: Parsed pipeline configuration. Uses ``data_dir``, ``results_dir``,
            ``seed`` and ``annotator``; optionally ``l0_min_dp``, ``l0_window``,
            ``l0_panel_extension`` (extra gene symbols to search alongside the SAC panel)
            and ``l0_region_pad``.
    """
    data_dir = Path(config["data_dir"])
    out_dir = Path(config["results_dir"]) / "l0_genomics"
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = _find_vcf(data_dir)
    # Fail before the ~70 s burden scan if the transcript policy is missing or refused --
    # but after _find_vcf, so a missing dataset is still reported as a missing dataset,
    # which is the more useful of the two errors.
    _transcripts.transcript_policy(config)

    LOGGER.info("scanning VCF for aneuploidy burden (this reads ~5M records)")

    result = compute_burden(
        vcf_path,
        min_dp=int(config.get("l0_min_dp", _burden.DEFAULT_MIN_DP)),
        window=int(config.get("l0_window", DEFAULT_WINDOW)),
    )

    payload = {
        "baseline_mean_abs_dev": result.baseline_mean_dev,
        "baseline_sd": result.baseline_sd,
        "inferred_sex": result.inferred_sex,
        "contamination_screen": result.contamination,
        "notes": result.notes,
        "per_contig": result.contigs,
        "sensitivity": {
            "median_depth": result.median_depth,
            "per_site_baf_sd": result.sigma,
            "min_detectable_mosaic_fraction": result.detection_limit,
            "min_detectable_mosaic_fraction_conservative": (
                None if result.detection_limit is None else result.detection_limit * 2.5
            ),
            "z_threshold": _burden.Z_THRESHOLD,
            "interpretation": (
                "Measured from this sample's own depth and window-to-window scatter, not "
                "configured. A negative result means no event above this fraction was "
                "detected -- it does not mean no event is present."
            ),
            "caveat": (
                "Model-dependent. The optimistic figure assumes the sample's systematic "
                "baseline excess and a mosaic shift add in shift space; dropping that "
                "assumption gives roughly the conservative figure, and a spike-in "
                "simulation at this depth put the crossing near 0.25. Calibrate against "
                "simulated spike-ins before quoting a single number."
            ),
        },
        "method": {
            "statistic": "mean |BAF - 0.5| at heterozygous biallelic SNVs",
            "baseline": "sample's own autosomal windows, median + scaled MAD",
            "fraction_estimate": (
                "noise-deconvolved: the folded-normal mean is inverted for the true BAF "
                "shift before conversion, because the raw statistic sits at the noise "
                "floor when the shift is zero and responds quadratically to small shifts"
            ),
            "window_bp": int(config.get("l0_window", DEFAULT_WINDOW)),
            "min_dp": int(config.get("l0_min_dp", _burden.DEFAULT_MIN_DP)),
            "citations": ["conlin2010 doi:10.1093/hmg/ddq003", "loh2018 doi:10.1038/s41586-018-0321-x"],
        },
        "seed": config["seed"],
    }
    (out_dir / "aneuploidy_burden.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    (out_dir / "windows.json").write_text(
        json.dumps([asdict(w) for w in result.windows], indent=2), encoding="utf-8"
    )
    LOGGER.info(
        "burden written: sex=%s contamination=%s baseline=%.4f",
        result.inferred_sex, result.contamination, result.baseline_mean_dev,
    )

    # Causal-gene half: panel alleles -> the two decision-D5 annotation passes -> zygosity
    # and phase -> biallelic configurations -> a candidate for gate G1. The panel is the
    # SAC set plus any configured extension; never filter to one gene a priori. Nothing
    # here writes config["causal_gene"] -- that is a human decision at G1.
    panel = tuple(dict.fromkeys((*SAC_PANEL, *config.get("l0_panel_extension", ()))))
    summary = _causal.call(vcf_path, config, out_dir, panel=panel)
    LOGGER.info(
        "causal-gene call written: verdict=%s, %d panel allele(s), %d candidate gene(s)",
        summary["verdict"], summary["n_alleles"], len(summary["candidate_genes"]),
    )
