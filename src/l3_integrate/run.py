"""L3 -- Integration.

Purpose
    Turn the channel outputs into one ranked, reasoned candidate list: harmonize drug
    identities, aggregate ranks, then reason over the result.

Inputs
    ``results/l2/channels.json`` -- the authority for which channels produced evidence in
    *this* run -- and the ``candidates.tsv`` of each channel it records as complete.

Outputs
    Under ``config["results_dir"]/l3/``:
      - ``candidates.tsv`` -- one row per harmonised molecule, with its aggregated score,
        per-channel ranks, supporting channels, convergence kind, and identifiers
      - ``integration.json`` -- method, parameters, counts, source provenance, caveats

    These rows must satisfy the **per-candidate field contract** in
    ``mngmt/decisions.md`` (decision D2), because L5 renders the candidate dossier from
    them and is forbidden to compute anything itself. The reasoning fields --
    ``rationale``, ``contradicting_evidence``, ``confidence`` -- are **not** emitted yet;
    they come from the Claude step, which is deliberately a separate stage so that a
    failure downstream never forces an expensive re-run. Until that stage exists this
    layer's output is explicitly incomplete against D2, and says so in ``integration.json``
    rather than presenting a partial contract as a finished one.

Guardrail
    Identity harmonization happens **before** aggregation. Aggregating on drug names would
    split the same compound across channels under different synonyms and destroy exactly
    the cross-channel convergence signal the design depends on. On the real data this is
    not hypothetical: channel B's top five were three molecules and two of their own salt
    forms.

    **``channels.json`` is the authority, never a directory listing.** A table on disk from
    an earlier run is not evidence from this one, and a channel that failed must not
    contribute silently.

    **No patient-derived file is read here.** Channel D's ``evidence.json`` states which
    diseases resemble the subject; this layer reads only ``candidates.tsv``, which carries
    identities and counts.

    No NC / ShareAlike-derived field may reach a redistributed output from this layer.
    Enrichment joins are by identifier into a segregated zone -- see
    ``src/l4_validate/sources.md``.
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

# Imported by name, not as `from ..l2_channels import run`: that package re-exports its
# own run() function under the same name, which shadows the submodule.
from ..l2_channels.run import STATUS_FILE as L2_STATUS_FILE
from ..l2_channels.run import output_root as l2_output_root
from . import aggregate as aggregate_mod
from . import harmonize

LOGGER = logging.getLogger(__name__)

CANDIDATES_FILE = "candidates.tsv"
INTEGRATION_FILE = "integration.json"
#: The column every channel writes, and the only one this layer needs from them.
IDENTITY_COLUMN = "chembl_id"

DEFAULTS = {"broad_channel_fraction": aggregate_mod.BROAD_CHANNEL_FRACTION}

CAVEATS = (
    "This layer ranks; it does not judge. Safety, paediatric suitability and the exclusion "
    "of genotoxic agents are L4's, and nothing here has been through them. A high rank in "
    "this table is not a recommendation and several of the top rows are expected to be "
    "excluded outright.",
    "Convergence between channels is reported by kind, not as a bare count. A channel that "
    "ranks most of the candidate pool cannot discriminate, so its agreement is weak "
    "evidence; the convergence column says whether a row was found by two discriminating "
    "channels, by one plus a broad one, or by broad channels only.",
    "The channels are not fully independent. All of them reach drug identity, and two of "
    "them reach drug indications or targets, through Open Targets and ChEMBL. Agreement "
    "between channels sharing a substrate is worth less than agreement between "
    "independent ones, and no count here corrects for that.",
    "Ranks are normalised within each channel, so a rank means 'in the top fraction of "
    "what this channel ranked'. Channels rank very different numbers of drugs, and a "
    "candidate absent from a channel contributes nothing rather than a penalty -- absence "
    "is missing evidence, not evidence against.",
    "The RRA score is a p-value against random ranking, Bonferroni-corrected across "
    "channels. It is not a probability that the drug works, and it inherits every bias in "
    "the channels that produced the ranks.",
    "Roughly 45% of candidates carry no RxCUI, because the public sources used here cover "
    "drugs marketed in the United States. An unresolved RxCUI means the crosswalk did not "
    "reach it, never that the drug is not approved somewhere.",
    "A drug's formulations are collapsed onto the parent molecule Open Targets names. A "
    "formulation whose parent is not recorded stays separate and may occupy its own rank.",
)


def output_dir(config: dict) -> Path:
    return Path(config["results_dir"]) / "l3"


def _settings(config: dict) -> dict:
    return {**DEFAULTS, **((config.get("l3") or {}))}


def load_channel_tables(config: dict) -> tuple:
    """Read the ranked table of every channel ``channels.json`` records as complete.

    Returns ``({channel: [chembl_id, ...]}, status)`` with each list in the channel's own
    rank order.

    Raises:
        FileNotFoundError: If no status record exists -- L2 has not run in this results
            directory, and guessing from directory contents is exactly what the record
            exists to prevent.
        ValueError: If a channel recorded complete has no readable table.
    """
    status_path = l2_output_root(config) / L2_STATUS_FILE
    if not status_path.is_file():
        raise FileNotFoundError(
            f"{status_path} does not exist, so there is no record of which channels "
            "produced evidence in this run. Run the l2_channels layer first; L3 will not "
            "infer channel status from whatever tables happen to be on disk.")
    status = json.loads(status_path.read_text(encoding="utf-8"))

    tables = {}
    for name in status["summary"]["complete"]:
        record = status["channels"][name]
        # output_dir is recorded relative to results_dir (see l2_channels/run.py), so the
        # path travels with the results directory rather than assuming where it sits.
        table = Path(config["results_dir"]) / record["output_dir"] / CANDIDATES_FILE
        if not table.is_file():
            raise ValueError(
                f"channel {name} is recorded complete in {status_path} but {table} is "
                "missing. The status record and the results directory disagree; re-run L2 "
                "rather than proceeding on a stale table.")
        with open(table, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows or IDENTITY_COLUMN not in rows[0]:
            raise ValueError(f"channel {name}'s table has no {IDENTITY_COLUMN!r} column; "
                             "every channel must emit one so identities can be harmonised")
        tables[name] = [r[IDENTITY_COLUMN] for r in rows if r.get(IDENTITY_COLUMN)]
    return tables, status


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def run(config: dict) -> None:
    """Execute layer L3.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir``, ``reference_dir``,
            ``enrichment_dir``, ``seed`` and ``l3``.

    Raises:
        FileNotFoundError: If L2 left no status record.
        ValueError: If fewer than two channels completed, or a recorded table is missing.
    """
    settings = _settings(config)
    out_dir = output_dir(config)
    out_dir.mkdir(parents=True, exist_ok=True)

    tables, status = load_channel_tables(config)
    if len(tables) < 2:
        raise ValueError(
            f"only {len(tables)} channel(s) completed, so there is nothing to aggregate "
            "across. Rank aggregation over one channel is that channel's ranking with "
            "extra steps; report it as such rather than presenting it as consensus.")
    LOGGER.info("aggregating %d channel(s): %s", len(tables), ", ".join(sorted(tables)))

    every_id = {cid for rows in tables.values() for cid in rows}
    identities, harmonise_stats, provenance = harmonize.build(config, every_id)

    harmonised = {channel: [identities[cid].key for cid in rows if cid in identities]
                  for channel, rows in tables.items()}
    results, aggregate_stats = aggregate_mod.aggregate(
        harmonised, broad_fraction=float(settings["broad_channel_fraction"]))

    channels = aggregate_stats["channels"]
    header = (["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "rxcui",
               "unii", "pharm_class_epc", "pharm_class_moa", "rxcui_resolution",
               "n_chembl_ids_collapsed", "rra_score", "rra_p_value",
               "n_channels_supporting", "supporting_channels", "convergence"]
              + [f"rank_{c}" for c in channels])
    by_key = {identity.key: identity for identity in identities.values()}
    rows = []
    for position, item in enumerate(results, start=1):
        public = by_key[item.key].publishable()
        rows.append([position, public["chembl_id"], public["drug_name"],
                     public["drug_type"], public["clinical_stage"], public["rxcui"],
                     public["unii"], public["pharm_class_epc"], public["pharm_class_moa"],
                     public["rxcui_resolution"], public["n_chembl_ids_collapsed"],
                     f"{item.score:.6g}", f"{item.p_value:.6g}",
                     item.n_channels_supporting, "|".join(item.channels),
                     item.convergence]
                    + [item.ranks.get(c, "") for c in channels])
    _write_tsv(out_dir / CANDIDATES_FILE, header, rows)

    (out_dir / INTEGRATION_FILE).write_text(json.dumps({
        "layer": "l3_integrate",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": config["seed"],
        "causal_gene": config.get("causal_gene"),
        "therapeutic_endpoint": config.get("therapeutic_endpoint"),
        "method": {
            "identity": "ChEMBL molecules collapsed onto the parent Open Targets names, "
                        "then cross-referenced to UNII (UniChem) and RxCUI (openFDA)",
            "aggregation": "Robust Rank Aggregation over per-channel normalised ranks "
                           "[kolde2012] doi:10.1093/bioinformatics/btr709; score is the "
                           "minimum Beta CDF over rank-order statistics, Bonferroni-"
                           "corrected by the number of channels",
            "absent_rank": "a channel that did not rank a candidate contributes nothing; "
                           "absence is missing evidence, not evidence against",
            "convergence": "reported by kind, not as a count -- a channel ranking more "
                           f"than {settings['broad_channel_fraction']:.0%} of the pool is "
                           "'broad' and its agreement alone is not convergence",
            "ranking": "ascending RRA score, then convergence kind, then support count, "
                       "then ChEMBL id",
        },
        "parameters": {k: settings[k] for k in DEFAULTS},
        "channels_used": {"from": str((l2_output_root(config)
                                       / L2_STATUS_FILE)),
                          "complete": status["summary"]["complete"],
                          "not_implemented": status["summary"]["not_implemented"],
                          "l2_finished_utc": status.get("finished_utc")},
        "counts": {**harmonise_stats,
                   **{k: v for k, v in aggregate_stats.items() if k != "channels"}},
        "field_contract": {
            "satisfied": sorted(set(header)),
            "outstanding": ["rationale", "contradicting_evidence", "confidence",
                            "confidence_basis", "evidence_grade", "references",
                            "safety_verdict"],
            "note": ("Decision D2's per-candidate contract is NOT yet satisfied. The "
                     "reasoning fields come from the Claude step and the safety fields "
                     "from L4; both are separate stages so that a failure in one never "
                     "forces an expensive re-run of the other. This file records the gap "
                     "rather than letting L5 discover it at render time."),
        },
        "sources": provenance,
        "caveats": list(CAVEATS),
    }, indent=2, sort_keys=True), encoding="utf-8")

    LOGGER.info("L3 written: %d candidate(s) from %d channel(s); convergence %s",
                len(results), len(tables), aggregate_stats["by_convergence"])
