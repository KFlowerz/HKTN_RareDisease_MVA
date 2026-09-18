"""L4 -- Validation and safety.

Purpose
    Take the ranked candidate list from L3 and decide what is defensible to publish: a
    hard paediatric safety triage first, then the recorded criteria, then a blinded
    internal benchmark of the pipeline itself.

Inputs
    ``results/l3/candidates.tsv`` -- harmonised, RxCUI-bearing, ranked. Drug annotations
    per ``sources.md``: openFDA label prose and Drugs@FDA product status, both
    US-Government public domain, for anything redistributed.

Outputs
    Under ``config["results_dir"]/l4/``:
      - ``survivors.tsv`` -- candidates no hard rule excluded, in L3's rank order
      - ``excluded.tsv`` -- everything dropped, with the rule that dropped it, the source
        field, and a quotable snippet. **This ships with the results.**
      - ``verdicts.json`` -- every rule's verdict for every candidate, survivors and
        excluded alike
      - ``validation.json`` -- method, counts per rule, provenance, caveats

    Verdict rows satisfy the **per-candidate field contract** in ``mngmt/decisions.md``
    (decision D2) for the safety half: one row per candidate per rule, carrying the rule
    name, the verdict, the source field, and a quotable snippet. The reasoning fields
    (``rationale``, ``contradicting_evidence``, ``confidence``) remain outstanding -- they
    come from L3's Claude step, which is not built.

Guardrail
    Exclusions are **hard gates, not weights**. A candidate that fails the genotoxic,
    paediatric or marketing rule is dropped, never merely down-ranked -- MVA is
    cancer-predisposing, and a ranked list is read as a recommendation no matter how it is
    captioned.

    The excluded table ships with the results. A safety filter whose decisions are
    invisible cannot be audited, and on this data the exclusions are the more informative
    half of the output.

    **Fail closed.** A candidate whose status cannot be determined is excluded as
    ``insufficient_evidence``. Most candidates are, because openFDA describes drugs
    marketed in the United States and the channels nominate from a much wider pool. That
    is reported as a coverage fact, not disguised as a safety finding.

    No NC/ShareAlike-derived annotation reaches a redistributed output path from this
    layer. Everything quoted here is openFDA's, which is public domain.
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from . import benchmark
from . import labels as labels_mod
from . import safety_triage

LOGGER = logging.getLogger(__name__)

SURVIVORS_FILE = "survivors.tsv"
EXCLUDED_FILE = "excluded.tsv"
VERDICTS_FILE = "verdicts.json"
VALIDATION_FILE = "validation.json"

#: Columns carried through from L3 so a survivor row stands on its own.
CARRIED = ("rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "rxcui",
           "unii", "pharm_class_epc", "pharm_class_moa", "rra_score",
           "n_channels_supporting", "supporting_channels", "convergence")

CAVEATS = (
    "Exclusion here is a statement about the evidence this pipeline could reach, not a "
    "clinical judgement about a drug. Most exclusions are 'insufficient evidence', which "
    "means openFDA carries no label this candidate could be matched to -- usually because "
    "the drug is not marketed in the United States.",
    "The rules are lexical. They match fixed patterns against label prose and cannot read "
    "a sentence the way a pharmacologist does. Section 13.1 is written in careful "
    "negatives, so negation is handled explicitly and conservatively: an affirmative "
    "finding of harm beats a negative finding in the same section.",
    "Surviving the triage is not a safety claim. It means no rule here found a reason to "
    "exclude, on a label that openFDA happens to carry, using patterns that can miss. "
    "Nothing in this pipeline is a clinical recommendation.",
    "openFDA states its data is not for clinical use and may be incomplete or inaccurate. "
    "Every verdict drawn from it inherits that.",
    "Paediatric use is read from Section 8.4 as written. A label establishing use in one "
    "age band and denying it in another is treated as an exclusion, because the denial is "
    "what binds for a child outside the approved band.",
    "Blood-brain-barrier penetration is recorded as not applicable, not scored: the "
    "documented phenotype does not establish a CNS requirement (decision D4).",
    "Drug interactions are not assessed. The NLM RxNav interaction API was discontinued "
    "on 2024-01-02 and no structured, licence-clean source with severity exists, so this "
    "layer makes no interaction claim rather than an unsourced one.",
)


def output_dir(config: dict) -> Path:
    return Path(config["results_dir"]) / "l4"


def load_candidates(config: dict) -> list:
    """L3's ranked table.

    Raises:
        FileNotFoundError: If L3 has not run. Triaging a candidate list that does not
            exist would produce an empty exclusions table, which reads as "nothing was
            excluded".
    """
    path = Path(config["results_dir"]) / "l3" / "candidates.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} does not exist. Run the l3_integrate layer first; L4 triages L3's "
            "output and has nothing to say without it.")
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _split(value: str) -> tuple:
    return tuple(v for v in (value or "").split("|") if v)


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        writer.writerows(rows)


def run(config: dict) -> None:
    """Execute layer L4.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir``, ``reference_dir``,
            ``seed``, and ``therapeutic_endpoint`` (for endpoint-specific criteria).

    Raises:
        FileNotFoundError: If L3 has not run.
        ValueError: If every candidate is excluded -- a result worth reporting, but not
            silently as an empty survivors table.
    """
    out_dir = output_dir(config)
    out_dir.mkdir(parents=True, exist_ok=True)
    endpoint = config.get("therapeutic_endpoint") or ""

    candidates = load_candidates(config)
    rxcuis = {r for c in candidates for r in _split(c.get("rxcui"))}
    uniis = {u for c in candidates for u in _split(c.get("unii"))}
    LOGGER.info("triaging %d candidate(s) against %d RxCUI(s) and %d UNII(s)",
                len(candidates), len(rxcuis), len(uniis))

    paths, provenance = labels_mod.ensure_datasets(config)
    label_index, label_stats = labels_mod.load_labels(paths["label"], rxcuis, uniis)
    marketing_index, marketing_stats = labels_mod.load_marketing(
        paths["drugsfda"], rxcuis, uniis)

    survivors, excluded, verdict_records = [], [], []
    rule_counts = {rule: {} for rule in safety_triage.TRIAGE_RULES}
    for candidate in candidates:
        rxcui, unii = _split(candidate.get("rxcui")), _split(candidate.get("unii"))
        label = labels_mod.for_identity(label_index, rxcui, unii)
        marketing = labels_mod.for_identity(marketing_index, rxcui, unii)
        classes = _split(candidate.get("pharm_class_epc")) + _split(
            candidate.get("pharm_class_moa"))

        verdicts = safety_triage.triage_candidate(
            label=label, marketing=marketing, pharm_classes=classes, endpoint=endpoint)
        for verdict in verdicts:
            counts = rule_counts[verdict.rule]
            counts[verdict.verdict] = counts.get(verdict.verdict, 0) + 1

        carried = {k: candidate.get(k, "") for k in CARRIED}
        verdict_records.append({**carried, "survived": safety_triage.survives(verdicts),
                                "verdicts": [v.as_row() for v in verdicts]})
        if safety_triage.survives(verdicts):
            survivors.append(carried)
        else:
            dropped = safety_triage.first_exclusion(verdicts)
            excluded.append({**carried, "excluded_by": dropped.rule,
                             "reason": dropped.reason,
                             "source_field": dropped.source_field,
                             "source_section": dropped.source_section,
                             "snippet": dropped.snippet})

    _write_tsv(out_dir / SURVIVORS_FILE, list(CARRIED),
               [[s[k] for k in CARRIED] for s in survivors])
    excluded_header = list(CARRIED) + ["excluded_by", "reason", "source_field",
                                       "source_section", "snippet"]
    _write_tsv(out_dir / EXCLUDED_FILE, excluded_header,
               [[e[k] for k in excluded_header] for e in excluded])
    (out_dir / VERDICTS_FILE).write_text(
        json.dumps({"candidates": verdict_records, "seed": config["seed"]},
                   indent=2, sort_keys=True), encoding="utf-8")

    (out_dir / VALIDATION_FILE).write_text(json.dumps({
        "layer": "l4_validate",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": config["seed"],
        "therapeutic_endpoint": endpoint,
        "method": {
            "hard_rules": list(safety_triage.HARD_RULES),
            "recorded_rules": [r for r in safety_triage.TRIAGE_RULES
                               if r not in safety_triage.HARD_RULES],
            "fail_closed": "a candidate whose status could not be determined is excluded "
                           "as insufficient_evidence; absence of a label is not evidence "
                           "of safety",
            "genotoxic_signals": "FDA Established Pharmacologic Class, an affirmative "
                                 "finding in SPL section 13.1, or an increased-malignancy "
                                 "statement in a boxed warning or warnings section; any "
                                 "one excludes",
            "negation": "negative genotoxicity findings are consulted only when no "
                        "affirmative finding is present in the same section",
            "marketing": "a molecule with no prescription or over-the-counter product on "
                         "Drugs@FDA is excluded; Open Targets' APPROVAL stage means 'was "
                         "approved', not 'is marketed'",
        },
        "counts": {
            "candidates": len(candidates),
            "survivors": len(survivors),
            "excluded": len(excluded),
            "by_rule": rule_counts,
            "excluded_by_rule": {rule: sum(1 for e in excluded if e["excluded_by"] == rule)
                                 for rule in safety_triage.HARD_RULES},
            "excluded_by_reason": {reason: sum(1 for e in excluded
                                               if e["reason"] == reason)
                                   for reason in sorted({e["reason"] for e in excluded})},
            **label_stats, **marketing_stats,
        },
        "sources": provenance,
        "caveats": list(CAVEATS),
    }, indent=2, sort_keys=True), encoding="utf-8")

    # The benchmark scores L3's ranking, not the survivors, so it runs whatever the triage
    # decided -- including the case where the triage excluded everything. What the
    # pipeline surfaced and what is safe to propose are separate questions.
    report = benchmark.run_blinded(config)
    LOGGER.info("L4 written: %d survivor(s), %d excluded; benchmark AUROC %.3f over %d "
                "positive(s)", len(survivors), len(excluded),
                report["combined"]["auroc"], report["combined"]["n_positives_found"])
    if not survivors:
        raise ValueError(
            f"every one of the {len(candidates)} candidate(s) was excluded. That is a "
            f"reportable result -- see {out_dir / EXCLUDED_FILE} for the rule that "
            "dropped each -- but it is not a shortlist, and the report must say so rather "
            "than showing an empty table.")
