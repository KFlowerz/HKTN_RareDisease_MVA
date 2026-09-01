"""L3 -- Integration.

Purpose
    Turn five heterogeneous channel outputs into one ranked, reasoned candidate list:
    harmonize drug identities, aggregate ranks, then reason over the result.

Inputs
    The per-channel ranked tables from L2.

Outputs
    A single ranked candidate table under ``config["results_dir"]/l3/``, each row
    carrying its RxCUI, per-channel ranks, the aggregated score, the number of channels
    supporting it, a synthesized rationale, contradicting evidence, and a calibrated
    confidence.

Guardrail
    Identity harmonization happens **before** aggregation. Aggregating on drug names
    would split the same compound across channels under different synonyms and destroy
    exactly the cross-channel convergence signal the design depends on.

    No NC / ShareAlike-derived field may reach a redistributed output from this layer.
    Enrichment joins are by RxCUI into a segregated zone -- see
    ``src/l4_validate/sources.md``.
"""

from __future__ import annotations


def run(config: dict) -> None:
    """Execute layer L3.

    Args:
        config: Parsed pipeline configuration. Uses ``channels``, ``results_dir``, and
            ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: load each enabled channel's table; call `harmonize.to_rxcui` on every
    # candidate; call `aggregate.rank_aggregate` over the harmonized tables; then call
    # `claude_reasoning.reason_over_candidates` on the top-N. Persist the intermediate
    # after each step -- the reasoning step costs API calls and must not be re-run
    # because a later step failed.
    raise NotImplementedError("l3_integrate.run is a scaffold stub")
