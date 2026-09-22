"""Check every quantitative claim in docs/report_track2.md against the artifacts it cites.

The report asserts that each number in it came from a named artifact under ``results/``.
This checks that, so the claim is enforced rather than promised. It catches the failure
that matters: a pipeline re-run moves a number and the prose silently keeps the old one.

Each entry pairs a value read from an artifact with the string the report must contain. A
claim fails if the artifact disagrees **or** if the text has drifted away from it.

Usage::

    python scripts/verify_report_claims.py     # exit 0 if every claim holds

Requires a completed run under ``results/``. Run it before submitting, and after any
re-run that could move a number.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "results"
REPORT = REPO_ROOT / "docs" / "report_track2.md"


def load(relative: str):
    path = RESULTS / relative
    if not path.is_file():
        raise SystemExit(
            f"{path} not found. This checks a completed run; execute the pipeline first "
            "(PYTHONHASHSEED=42 python -m src.pipeline).")
    return json.loads(path.read_text(encoding="utf-8"))


def build_claims() -> list:
    """``(label, text the report must contain, value read from the artifact)``."""
    validation = load("l4/validation.json")["counts"]
    integration = load("l3/integration.json")["counts"]
    benchmark = load("l4/benchmark/recovery.json")
    dossier = load("l5/report.json")["counts"]
    burden = load("l0_genomics/aneuploidy_burden.json")
    rationales = load("l3/reasoning/rationales.json")
    mva = load("l1_target/module.json")["counts"]
    cftr = load("scalability_cftr/l1_target/module.json")["counts"]
    tier1 = [r for r in rationales if r["evidence_tier"] == "literature"][0]
    combined = benchmark["combined"]
    sensitivity = burden["sensitivity"]

    return [
        ("candidates nominated", "1,963", validation["candidates"] == 1963),
        ("candidates excluded", "1,880", validation["excluded"] == 1880),
        ("survivors", "83", validation["survivors"] == 83),
        ("tier 1 count", "Tier 1 — literature-supported: 1", dossier["tier1_literature"] == 1),
        ("tier 2 count", "82 candidates", dossier["tier2_network_only"] == 82),
        ("candidate pages", "83", dossier["candidate_pages"] == 83),
        ("genotoxic rule", "1,623",
         validation["excluded_by_rule"]["genotoxic_or_cancer_risk"] == 1623),
        ("paediatric rule", "| 256 |", validation["excluded_by_rule"]["pediatric_use"] == 256),
        ("insufficient evidence", "1,231",
         validation["excluded_by_reason"]["insufficient_evidence"] == 1231),
        ("aneugenic finding", "| 14 |",
         validation["excluded_by_reason"]["aneugenic_finding"] == 14),
        ("label records", "262,883", validation["label_records"] == 262883),
        ("label records matched", "54,930", validation["label_records_matched"] == 54930),
        ("proximity ranked", "1,956", integration["per_channel_ranked"]["proximity"] == 1956),
        ("phenotype ranked", "| 25 |", integration["per_channel_ranked"]["phenotype"] == 25),
        ("prior ranked", "| 8 |", integration["per_channel_ranked"]["prior"] == 8),
        ("discriminating convergence", "**0**",
         integration["by_convergence"]["discriminating"] == 0),
        ("mixed convergence", "| 26 |", integration["by_convergence"]["mixed"] == 26),
        ("no convergence", "1,937", integration["by_convergence"]["none"] == 1937),
        ("benchmark AUROC", "0.8147", round(combined["auroc"], 4) == 0.8147),
        ("benchmark CI low", "0.5809", combined["auroc_ci"]["ci_low"] == 0.5809),
        ("benchmark CI high", "0.9533", combined["auroc_ci"]["ci_high"] == 0.9533),
        ("benchmark positives", "8 of 1,963", combined["n_positives_found"] == 8),
        ("prior is circular", "excluded as circular",
         benchmark["circular_channels"] == ["prior"]),
        ("tier 1 confidence", "0.132", tier1["confidence"] == 0.132),
        ("tier 1 grade", "*in vitro*", tier1["evidence_grade_supplied"] == "in_vitro"),
        ("tier 1 citation", "10.1038/s41467-024-52176-x",
         tier1["citation_keys"] == ["DOI:10.1038/s41467-024-52176-x"]),
        ("burden window", "10 Mb", burden["method"]["window_bp"] == 10_000_000),
        ("burden min depth", "minimum depth 10", burden["method"]["min_dp"] == 10),
        ("detectable fraction", "0.098",
         round(sensitivity["min_detectable_mosaic_fraction"], 3) == 0.098),
        ("conservative bound", "0.245",
         round(sensitivity["min_detectable_mosaic_fraction_conservative"], 3) == 0.245),
        ("z threshold", "z = 5", sensitivity["z_threshold"] == 5.0),
        ("MVA module seeds", "| 32 |", mva["seeds"] == 32),
        ("CFTR module seeds", "| 9 |", cftr["seeds"] == 9),
        ("MVA target split", "84 / 116", (mva["upstream"], mva["downstream"]) == (84, 116)),
        ("CFTR target split", "54 / 146", (cftr["upstream"], cftr["downstream"]) == (54, 146)),
    ]


def main() -> None:
    if not REPORT.is_file():
        raise SystemExit(f"{REPORT} not found")
    report = REPORT.read_text(encoding="utf-8")

    failures = []
    for label, needle, artifact_agrees in build_claims():
        in_text = needle in report
        ok = artifact_agrees and in_text
        print(f"  {'OK  ' if ok else 'FAIL'} {label}")
        if not ok:
            failures.append(
                f"{label}: artifact agrees={artifact_agrees}, text present={in_text} "
                f"(looked for {needle!r})")

    print()
    if failures:
        print("The report and the artifacts disagree:")
        for line in failures:
            print("  -", line)
        raise SystemExit(1)
    print(f"all {len(build_claims())} quantitative claims verified against results/")


if __name__ == "__main__":
    main()
