"""L0 -- Genomics.

Purpose
    Turn the subject's whole-genome sequencing data into the scientific premise of the
    run: which SAC gene is biallelically disrupted, what the variant does to the
    protein, and which pathway module is affected. Additionally quantify **aneuploidy
    burden** from per-chromosome BAM read depth -- MVA's phenotype literally is mosaic
    aneuploidy, so the degree of missegregation is itself a feature, not just context.

Inputs
    WGS **VCF** (small variants, incl. gVCF) and **BAM** (alignments), read from
    ``config["data_dir"]``. The Track-1 validated causal variant, for reconciliation.

Outputs
    Written under ``config["results_dir"]``:
      - causal gene symbol + the biallelic variant pair supporting the call
      - predicted variant effect (VEP / OpenCRAVAT-style consequence, LoF confidence)
      - the affected pathway module identifier handed to L1
      - a per-chromosome aneuploidy burden vector + a scalar summary

Guardrail
    Reads **only** from ``config["data_dir"]`` and writes **only** under
    ``config["results_dir"]`` -- both gitignored. No patient-derived bytes (variant
    coordinates, read depths, sample identifiers) may be written anywhere else, logged
    at INFO, or embedded in a report figure at a resolution that permits
    re-identification. This is a real minor's genome.

    The causal gene is a *finding*, never an assumption: if reconciliation with Track 1
    disagrees, fail loudly rather than proceeding on the config value.
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
    #  1. Ingest: open the VCF with `pysam.VariantFile`; restrict to the SAC gene panel
    #     (BUB1B, CEP57, TRIP13, BUB1, BUB3, CEP192) plus a configurable extension set.
    #     Do not filter to a single gene a priori.
    #  2. Annotate: run VEP / OpenCRAVAT-style consequence prediction; keep LoF calls
    #     (nonsense, frameshift, canonical splice) and rank by predicted severity.
    #  3. Zygosity: identify *biallelic* configurations (hom-alt, or comp-het phased or
    #     inferred from parental data if present). MVA is autosomal recessive -- a single
    #     het LoF is not a causal call.
    #  4. Aneuploidy burden: from the BAM, compute per-chromosome mean read depth in
    #     fixed bins, normalize to the autosomal median, and derive a mosaicism estimate
    #     per chromosome. Correct for GC bias and mappability before calling gains/losses.
    #     Emit the burden vector as an L2/L4 feature.
    # Finally: reconcile the called gene against the Track-1 validated variant and raise
    # on disagreement rather than silently preferring one source.
    raise NotImplementedError("l0_genomics.run is a scaffold stub")
