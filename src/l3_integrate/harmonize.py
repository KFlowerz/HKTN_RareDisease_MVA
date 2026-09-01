"""L3 -- Drug identity harmonization on the RxNorm RxCUI backbone.

Purpose
    Collapse the many identifier spaces the five channels speak -- brand names, INNs,
    ChEMBL IDs, DrugBank IDs, PubChem CIDs, LINCS perturbagen names -- onto a single
    canonical key: the **RxNorm RxCUI**. Then cross-map outward to DrugBank / ChEMBL /
    UNII / ATC / PubChem / NDC via UniChem for downstream annotation.

Inputs
    Per-channel candidate tables with heterogeneous identifiers.

Outputs
    The same rows keyed by RxCUI, plus a crosswalk table and an explicit
    ``unresolved`` list.

Guardrail
    RxCUI is the **backbone**, chosen because RxNorm is US-Government public domain and
    therefore safe in a CC-BY-4.0 output. Cross-mapped identifiers from NC/ShareAlike
    sources (DrugBank, ChEMBL) are for *joining into the segregated enrichment zone
    only* -- do not redistribute their annotation content.

    Never guess a mapping. Salt forms, racemates/enantiomers, and combination products
    are genuinely distinct entities; silently merging them would corrupt the convergence
    count that drives the ranking. Unresolved candidates go to the ``unresolved`` list
    for human review, not to a fuzzy-matched RxCUI.
"""

from __future__ import annotations


def to_rxcui(config: dict) -> None:
    """Resolve candidate drug identifiers to RxNorm RxCUIs.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: resolve via the RxNorm REST API (`/rxcui?name=`, `/approximateTerm` only with
    # a strict score threshold and a recorded match score); normalize to ingredient level
    # (TTY=IN/PIN) so brand and generic collapse but salt forms stay distinguishable;
    # cross-map with UniChem for DrugBank/ChEMBL/PubChem/UNII; and cache every lookup
    # keyed by input string so a re-run is deterministic and offline-replayable.
    #
    # Emit `rxcui`, `rxcui_name`, `match_method`, `match_score`, and the source
    # identifier on every row. Collect failures in `unresolved.csv` -- do not drop them
    # silently, since a systematically unresolvable channel is a bug, not a null result.
    raise NotImplementedError("harmonize.to_rxcui is a scaffold stub")
