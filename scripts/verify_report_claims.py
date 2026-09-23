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
#: The plain-language explainer restates the same figures for a non-technical reader.
#: A number that drifts there is exactly as wrong as one that drifts in the report, and
#: it is the document most likely to be forwarded on its own.
LAY = REPO_ROOT / "docs" / "how_it_works.md"


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
        # The "others" row is a remainder, so it is the one cell no per-reason check
        # covers -- and it was wrong by one until 2026-09-22, because nothing made the
        # column add up. This asserts the arithmetic, not a transcribed number.
        ("others row closes the table", "| others | 22 |",
         sum(validation["excluded_by_reason"].values())
         - sum(validation["excluded_by_reason"][reason] for reason in (
             "insufficient_evidence", "safety_not_established", "clastogenic_finding",
             "increased_tumours", "mutagenic_finding", "malignancy_risk",
             "positive_genotoxicity_assay", "cytotoxic_class", "secondary_malignancy",
             "aneugenic_finding")) == 22),
        ("reason counts sum to the exclusion total", "**1,880 of 1,963 candidates",
         sum(validation["excluded_by_reason"].values()) == validation["excluded"]),
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


def build_lay_claims() -> list:
    """Claims the plain-language explainer makes, in the words it uses.

    Deliberately a smaller set: the explainer quotes headline figures and rounds the
    sensitivity to whole percentages for a lay reader. What is checked is that each
    rounded statement still matches the artifact it came from.
    """
    validation = load("l4/validation.json")["counts"]
    integration = load("l3/integration.json")["counts"]
    dossier = load("l5/report.json")["counts"]
    sensitivity = load("l0_genomics/aneuploidy_burden.json")["sensitivity"]
    mva = load("l1_target/module.json")["counts"]
    call = load("l0_genomics/causal_gene_call.json")

    proximity_share = round(
        integration["per_channel_ranked"]["proximity"] / validation["candidates"] * 100, 1)
    silent_genes = [gene for gene, per in call["per_gene"].items()
                    if not per["configurations"]]

    return [
        ("nominated", "1,963", validation["candidates"] == 1963),
        ("excluded", "1,880", validation["excluded"] == 1880),
        ("genotoxic rule", "1,623",
         validation["excluded_by_rule"]["genotoxic_or_cancer_risk"] == 1623),
        ("paediatric rule", "256", validation["excluded_by_rule"]["pediatric_use"] == 256),
        ("no-label bucket", "**1,231**",
         validation["excluded_by_reason"]["insufficient_evidence"] == 1231),
        # Bound to the sentence, not the bare digits: "14" alone would match almost
        # anything and would pass while saying nothing.
        ("aneugenic count", "Fourteen drugs were removed",
         validation["excluded_by_reason"]["aneugenic_finding"] == 14),
        ("proximity ranked", "1,956",
         integration["per_channel_ranked"]["proximity"] == 1956),
        ("proximity share", "99.6%", proximity_share == 99.6),
        ("no convergence", "| **0** |", integration["by_convergence"]["discriminating"] == 0),
        ("tier 1", "| **1** |", dossier["tier1_literature"] == 1),
        ("tier 2", "| **82** |", dossier["tier2_network_only"] == 82),
        ("module size", "200 proteins", mva["module"] == 200),
        ("one gene survives", "exactly **one**", call["candidate_genes"] == ["BUB1B"]),
        ("five genes silent", "The other\nfive produced nothing", len(silent_genes) == 5),
        ("sensitivity rounds to 10%", "roughly **10%**",
         round(sensitivity["min_detectable_mosaic_fraction"] * 100) == 10),
        ("conservative bound rounds to 25%", "**25%**",
         round(sensitivity["min_detectable_mosaic_fraction_conservative"] * 100) in (24, 25)),
    ]


def _check(document: Path, claims: list, failures: list) -> None:
    if not document.is_file():
        raise SystemExit(f"{document} not found")
    text = document.read_text(encoding="utf-8")
    print(f"{document.relative_to(REPO_ROOT)}")
    for label, needle, artifact_agrees in claims:
        # Prose wraps, so a needle spanning a line break is matched on collapsed space.
        in_text = needle in text or " ".join(needle.split()) in " ".join(text.split())
        ok = artifact_agrees and in_text
        print(f"  {'OK  ' if ok else 'FAIL'} {label}")
        if not ok:
            failures.append(
                f"{document.name} / {label}: artifact agrees={artifact_agrees}, "
                f"text present={in_text} (looked for {needle!r})")
    print()


def main() -> None:
    failures: list = []
    report_claims, lay_claims = build_claims(), build_lay_claims()
    _check(REPORT, report_claims, failures)
    _check(LAY, lay_claims, failures)

    if failures:
        print("The documents and the artifacts disagree:")
        for line in failures:
            print("  -", line)
        raise SystemExit(1)
    print(f"all {len(report_claims) + len(lay_claims)} quantitative claims "
          "verified against results/")


if __name__ == "__main__":
    main()
