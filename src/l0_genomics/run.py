"""L0 -- Genomics.

Purpose
    Turn the subject's whole-genome sequencing data into the scientific premise of the
    run: which SAC gene is biallelically disrupted, what the variant does to the
    protein, and which pathway module is affected. Additionally quantify **aneuploidy
    burden** per chromosome -- MVA's phenotype literally is mosaic aneuploidy, so the
    degree of missegregation is itself a feature, not just context.

Inputs
    A single-sample WGS **VCF** (bgzipped, with tabix index), read from
    ``config["data_dir"]``. The Track-1 validated causal variant, for reconciliation.

    **There is no BAM.** The gated dataset ships raw reads (``.fastq.gz``) and a called
    VCF, not alignments -- verified against dataset revision
    ``59e322d27f399006b398d366d33e703e48a29914`` on 2026-09-02; see ``DATA.md``. Any
    depth-based method requiring alignments would first have to align ~85 GB of FASTQ,
    which is out of scope. This layer is therefore **VCF-only**.

Outputs
    Written under ``config["results_dir"]``:
      - causal gene symbol + the biallelic variant pair supporting the call
      - predicted variant effect (VEP / OpenCRAVAT-style consequence, LoF confidence)
      - the affected pathway module identifier handed to L1
      - a per-chromosome aneuploidy burden vector + a scalar summary

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

Guardrail
    Reads **only** from ``config["data_dir"]`` and writes **only** under
    ``config["results_dir"]`` -- both outside the repository and registered in
    ``docs/data_custody.md``. No patient-derived bytes (variant coordinates, allele
    depths, sample identifiers) may be written anywhere else, logged at INFO, or embedded
    in a report figure at a resolution that permits re-identification. This is a real
    minor's genome, and the source filenames themselves embed a lab accession and a
    sequencer flowcell identifier -- never propagate them into an artifact.

    The causal gene is a *finding*, never an assumption: if reconciliation with Track 1
    disagrees, fail loudly rather than proceeding on the config value.

    BAF is uninformative where a chromosome has too few heterozygous sites to populate the
    distribution. Emit an explicit ``insufficient_sites`` verdict for such a chromosome --
    never a burden of zero, which would read as "no aneuploidy detected" when the truth is
    "not measurable here".
"""

from __future__ import annotations


def run(config: dict) -> None:
    """Execute layer L0.

    Args:
        config: Parsed pipeline configuration. Uses ``data_dir``, ``results_dir``,
            ``seed``, and (for cross-checking only) ``causal_gene``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: implement L0 as four steps.
    #  1. Ingest: open the VCF with `pysam.VariantFile` (bgzipped + .tbi, so region
    #     queries work); restrict to the SAC gene panel (BUB1B, CEP57, TRIP13, BUB1,
    #     BUB3, CEP192) plus a configurable extension set. Do not filter to a single
    #     gene a priori.
    #  2. Annotate: run VEP / OpenCRAVAT-style consequence prediction; keep LoF calls
    #     (nonsense, frameshift, canonical splice) and rank by predicted severity.
    #  3. Zygosity: identify *biallelic* configurations (hom-alt, or comp-het phased or
    #     inferred). MVA is autosomal recessive -- a single het LoF is not a causal call.
    #     Note the dataset is single-sample: no parental data, so comp-het phasing must
    #     rely on read-backed phasing or be reported as unphased with that caveat.
    #  4. Aneuploidy burden (VCF-only, see the method note above): select high-confidence
    #     biallelic SNVs with adequate `FORMAT/DP`; compute BAF = AD[alt]/(AD[ref]+AD[alt])
    #     at heterozygous sites; per chromosome, characterize the BAF distribution
    #     (deviation from 0.5, modality, dispersion) and convert the deviation into a
    #     mosaic-fraction estimate. Emit the per-chromosome vector plus a scalar summary
    #     as an L2/L4 feature. Carry median normalized `FORMAT/DP` per chromosome as a
    #     secondary signal only. Chromosomes with too few informative het sites are
    #     `insufficient_sites`, never zero.
    # Finally: reconcile the called gene against the Track-1 validated variant and raise
    # on disagreement rather than silently preferring one source.
    raise NotImplementedError("l0_genomics.run is a scaffold stub")
