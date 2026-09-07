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
    #  4. Aneuploidy burden (VCF-only, see the method note above):
    #     a. Restrict to the primary contigs -- 1-22, X, Y with no `chr` prefix. The VCF
    #        declares 2,580 contigs; decoys and alts must not enter the burden vector.
    #     b. Keep FILTER=PASS biallelic SNVs with adequate `FORMAT/DP` (>=10 is workable;
    #        make the floor configurable and record it). Exclude PAR and low-complexity
    #        regions before anything is scored.
    #     c. BAF = AD[alt] / (AD[ref] + AD[alt]) at heterozygous sites.
    #     d. Build the per-sample diploid baseline from the autosomal bulk, then score
    #        each segment against it -- an internal null, not a literature constant.
    #     e. Segment along each chromosome; report segments, not chromosome-wide means.
    #        Convert deviation magnitude into a mosaic-fraction estimate with an interval,
    #        not a point value: at low fractions the estimate is not well determined.
    #     f. Infer sex from X heterozygosity and branch. Do not apply the autosomal model
    #        to a hemizygous contig; record the exclusion rather than emitting a number.
    #     g. Rule out contamination before calling anything mosaic -- low-level sample
    #        contamination and low-fraction mosaicism look alike in BAF, and calling the
    #        wrong one would misstate the central feature of this child's disease.
    #     Carry median normalized `FORMAT/DP` per chromosome as a secondary signal only.
    #     Segments with too few informative het sites are `insufficient_sites`, never zero.
    # Finally: emit the causal-gene call with its supporting and rejected evidence. There
    # is no Track-1 result to reconcile against; if one is published later, reconcile then
    # and raise on disagreement rather than silently preferring either source.
    raise NotImplementedError("l0_genomics.run is a scaffold stub")
