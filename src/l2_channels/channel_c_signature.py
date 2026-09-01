"""L2 Channel C -- Signature reversion (on a PROXY signature).

Purpose
    Nominate drugs whose transcriptional consequences *reverse* the disease signature --
    the classic connectivity-map logic.

    **There is no patient RNA-seq.** The subject's dataset is WGS (VCF/BAM) plus
    phenotype; no expression data exists to build a real disease signature from. This
    channel therefore runs on a **PROXY** signature, one of:

      1. LINCS L1000 **knockdown of the causal gene** (e.g. shRNA/CRISPR BUB1B) in the
         most phenotypically relevant available cell line, or
      2. a **curated aneuploidy-response gene set** (aneuploidy/CIN signatures from the
         literature) scored as an up/down list.

Inputs
    L0 causal gene; LINCS L1000 (or CLUE) signatures; the curated aneuploidy gene set.

Outputs
    Ranked candidates with a reversal/connectivity score, the identity of the proxy
    signature used, and the cell line and perturbation it came from.

Guardrail
    **The proxy substitution must travel with the result.** Every row this channel emits,
    and every figure or table built from it downstream, must name the proxy and its
    limits. A knockdown in an immortalized cell line is not a child's tissue; a curated
    aneuploidy set is not this child's aneuploidy. Reporting a connectivity score without
    that caveat would overstate the evidence -- this is the single most over-claimable
    channel in the pipeline.

    Never label proxy-derived output as a patient signature, and never impute patient
    expression from genotype to manufacture one.
"""

from __future__ import annotations


def generate(config: dict) -> None:
    """Generate candidates by signature reversion against a proxy signature.

    Args:
        config: Parsed pipeline configuration; uses the L0 causal gene, ``seed``, and
            ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: select the proxy signature -- prefer a LINCS L1000 knockdown of the causal
    # gene if one exists at adequate quality (check signature strength / replicate
    # correlation), else fall back to the curated aneuploidy gene set; record WHICH was
    # used in the output, per row.
    #
    # Then score every LINCS compound perturbagen for reversal (weighted connectivity
    # score / reverse-GSEA), aggregate across cell lines and doses with an explicit
    # summarization rule, and rank by the strength of the reversal.
    #
    # Emit `proxy_signature_id`, `proxy_type`, `cell_line`, and a `proxy_caveat` string
    # on every row so L3/L5 cannot drop the caveat by accident.
    raise NotImplementedError("channel_c_signature.generate is a scaffold stub")
