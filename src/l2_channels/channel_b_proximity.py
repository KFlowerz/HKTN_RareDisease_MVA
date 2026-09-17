"""L2 Channel B -- Network proximity.

Purpose
    Score drugs by the shortest-path proximity of their target sets to the L1 disease
    module in the human interactome (Guney/Barabasi). Drugs whose targets sit inside or
    adjacent to the module are nominated.

Inputs
    L1 disease module (node set); a drug-target bipartite set (targets resolved to the
    same protein identifier space as the interactome); the interactome itself.

Outputs
    Ranked candidates with the raw proximity ``d_c``, a degree-preserving null
    distribution, and the resulting z-score.

Guardrail
    This channel exists precisely because it **does not require the disease to appear in
    any knowledge graph** -- it needs only a target set and a module. That property is
    what makes it trustworthy for a hyper-rare disease, so it must never be
    reimplemented in terms of a disease node.

    Proximity must always be reported against a **degree-preserving** random null.
    Well-studied hub proteins are close to everything; a raw distance without the
    degree-matched z-score is a popularity measurement, not evidence.

    **The drug-target data is licence-restricted and stays segregated.** Targets come from
    :mod:`src.l2_channels.enrichment`, whose content may not be redistributed; only the
    fields ``enrichment.PUBLISHABLE_FIELDS`` names reach the table written here. Every row
    is built through :func:`~src.l2_channels.enrichment.publishable`, and a test asserts
    the output carries no target identifier.

    **This channel annotates; it does not exclude.** A drug scored here is nominated for
    L3, not recommended: L4 owns every safety and eligibility decision, and nothing in
    this file may filter on toxicity, age, or indication.

    No patient data is read or written here -- the channel works from the L1 module and
    public resources, and never opens ``data_dir``.
"""

from __future__ import annotations

import json
import logging
import random
from pathlib import Path

from ..l1_target import module as l1_module, sources
from . import enrichment, proximity

LOGGER = logging.getLogger(__name__)

CAVEATS = (
    "A proximity z-score measures network position, not pharmacology. It says a drug's "
    "targets sit closer to the module than degree-matched proteins do; it says nothing "
    "about whether engaging them helps, and nothing about safety.",
    "Drug-target coverage is uneven. A well-studied drug has more annotated targets than "
    "a poorly-studied one, and a drug with no annotated target cannot be scored at all -- "
    "so absence from this ranking is absence of annotation, not absence of effect.",
    "Targets are annotated per drug, not per dose or tissue. A target the drug reaches "
    "only at concentrations never achieved clinically counts the same as its primary one.",
    "The module is L1's top-N cut. Drugs are scored against that boundary, so a drug just "
    "outside a larger module would score differently; the module parameters travel in "
    "l1_target/module.json.",
    "The interactome is restricted to its largest connected component, because distance "
    "between components is undefined. Targets outside it are reported per drug and not "
    "scored.",
)


def _module_nodes(results_dir: Path) -> list:
    """The L1 module's proteins, read from L1's artifact.

    Raises:
        FileNotFoundError: If L1 has not run. The channel needs a module; inventing one
            would make the score meaningless while still producing a ranking.
        ValueError: If the module is empty.
    """
    path = results_dir / "l1_target" / "module_nodes.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found -- Channel B scores proximity to L1's disease module, so "
            "L1 must run first (python -m src.pipeline --only l1_target)"
        )
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        column = header.index("protein")
        nodes = [line.rstrip("\n").split("\t")[column] for line in handle if line.strip()]
    if not nodes:
        raise ValueError(f"{path} lists no module protein -- every distance would be "
                         "undefined, and the ranking would be over nothing")
    return nodes


def _module_score_min(results_dir: Path, config: dict) -> int:
    """The interaction cutoff L1 actually used, read from the module it wrote.

    Deriving it from config instead would let the module and the graph it is scored
    against disagree after an edit to ``l1.string_score_min`` -- silently, because both
    values are individually valid. The recorded one is the only one that matches the
    module's own boundary.
    """
    path = results_dir / "l1_target" / "module.json"
    if path.is_file():
        recorded = json.loads(path.read_text(encoding="utf-8")).get("parameters", {})
        if "string_score_min" in recorded:
            return int(recorded["string_score_min"])
    LOGGER.warning("%s records no string_score_min; falling back to config", path)
    return int((config.get("l1") or {}).get("string_score_min", l1_module.DEFAULT_SCORE_MIN))


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def generate(config: dict) -> None:
    """Generate candidates by network proximity to the disease module.

    Args:
        config: Parsed pipeline configuration; uses the L1 disease module, ``seed``, and
            ``results_dir``.
    """
    settings = ((config.get("l2") or {}).get("channel_b") or {})
    permutations = int(settings.get("permutations", proximity.DEFAULT_PERMUTATIONS))
    min_bin = int(settings.get("min_bin_size", proximity.DEFAULT_MIN_BIN_SIZE))

    results_dir = Path(config["results_dir"])
    out_dir = results_dir / "l2" / "channel_b_proximity"
    out_dir.mkdir(parents=True, exist_ok=True)
    module_nodes = _module_nodes(results_dir)
    score_min = _module_score_min(results_dir, config)

    paths, interactome_provenance = sources.ensure_sources(config)
    LOGGER.info("building interactome at confidence >= %d", score_min)
    graph = proximity.largest_component(l1_module.build_graph(
        sources.load_edges(paths["links"], score_min=score_min), score_min=score_min))

    # The crosswalk is built from STRING's own aliases (CC BY 4.0), so drug identity is
    # the only thing taken from the restricted zone.
    records, drug_provenance, stats = enrichment.drug_records(
        config, enrichment.ensembl_index(paths["aliases"]))
    # An empty channel is a scientific claim -- "no drug is near this module" -- and every
    # way of reaching zero here is a broken input instead: a crosswalk that matched
    # nothing, a release missing the clinical-stage field under approved_only, a module
    # whose proteins are absent from the graph. Say so rather than writing an empty table
    # that reads like a finding (see l2_channels/run.py).
    if not records:
        raise ValueError(
            "no drug has a target in the interactome's identifier space. That is a broken "
            f"crosswalk, not a finding: {stats['drugs_with_a_mechanism']} drug(s) carry a "
            f"mechanism and {stats['drugs_with_unmappable_targets_only']} lost every "
            "target in mapping. Check the Ensembl-to-protein index and "
            "l2.channel_b.approved_only."
        )
    in_graph = sum(1 for n in module_nodes if n in graph)
    if not in_graph:
        raise ValueError(
            f"none of the {len(module_nodes)} module protein(s) is in the interactome's "
            "largest component, so every distance is undefined. The module and the graph "
            "were probably built from different releases or identifier spaces."
        )
    LOGGER.info("scoring %d drug(s) against a %d-node module", len(records),
                len(module_nodes))

    distance = proximity.distances_from_module(graph, module_nodes)
    pools = proximity.degree_bins(graph, min_bin_size=min_bin)
    rng = random.Random(config["seed"])
    scored = {r.chembl_id: proximity.proximity(distance, pools, r.targets, rng=rng,
                                               permutations=permutations)
              for r in records}
    by_id = {r.chembl_id: r for r in records}

    rows = []
    for position, (chembl_id, result) in enumerate(proximity.rank(scored), start=1):
        # Every published field passes through the enrichment zone's whitelist; the target
        # set that produced the score never appears here.
        public = enrichment.publishable(by_id[chembl_id])
        rows.append([
            position, public["chembl_id"], public["name"], public["drug_type"],
            public["clinical_stage"],
            "" if result.z is None else f"{result.z:.4f}",
            "" if result.d_c is None else f"{result.d_c:.4f}",
            "" if result.null_mean is None else f"{result.null_mean:.4f}",
            "" if result.null_sd is None else f"{result.null_sd:.4f}",
            result.n_targets, result.n_usable, result.n_in_module, result.note,
        ])
    _write_tsv(out_dir / "candidates.tsv",
               ["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "z",
                "d_c", "null_mean", "null_sd", "n_targets", "n_targets_scored",
                "n_targets_in_module", "note"],
               rows)

    scored_count = sum(1 for r in scored.values() if r.z is not None)
    if not scored_count:
        raise ValueError(
            f"{len(scored)} drug(s) were scored and none received a z-score. A ranking "
            "without a null is a popularity measurement, so the channel produces nothing "
            "rather than a table of raw distances."
        )
    payload = {
        "channel": "channel_b_proximity",
        "method": {
            "measure": "closest: mean shortest-path distance from each target to the "
                       "nearest module protein [guney2016] doi:10.1038/ncomms10331",
            "distance": "unweighted hops on the largest connected component",
            "null": "each target resampled from proteins of comparable degree; the module "
                    "is held fixed because it is the same for every drug",
            "ranking": "ascending z -- more negative is closer than degree predicts",
        },
        "parameters": {
            "permutations": permutations,
            "min_bin_size": min_bin,
            "string_score_min": score_min,
            "approved_only": stats["approved_only"],
            "open_targets_release": str(settings.get("open_targets_release", "26.06")),
        },
        "counts": {
            "graph_nodes": graph.number_of_nodes(),
            "graph_edges": graph.number_of_edges(),
            "module_nodes": len(module_nodes),
            "module_nodes_in_graph": in_graph,
            "drugs_scored": scored_count,
            "drugs_without_a_z_score": len(scored) - scored_count,
            **stats,
        },
        "sources": {
            "interactome": interactome_provenance,
            "drug_targets": drug_provenance,
        },
        "licensing": {
            "published_fields": list(enrichment.PUBLISHABLE_FIELDS),
            "note": (
                "Drug-target content is licence-restricted and is not redistributed: it "
                "stays in the enrichment zone and only the fields above reach this "
                "output. See src/l4_validate/sources.md, Table 5."
            ),
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }
    (out_dir / "channel.json").write_text(json.dumps(payload, indent=2, sort_keys=True),
                                          encoding="utf-8")
    LOGGER.info("channel B written: %d drug(s) scored, %d without a z-score",
                scored_count, len(scored) - scored_count)
