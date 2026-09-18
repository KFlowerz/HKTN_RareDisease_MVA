"""Tests for the L4 layer and its blinded benchmark.

Invented drugs and identifiers throughout. The openFDA downloads are stubbed, so nothing
here touches the network.
"""

from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path

import pytest

from src.l4_validate import benchmark
from src.l4_validate import labels as labels_mod

l4 = importlib.import_module("src.l4_validate.run")

CARC = "carcinogenesis_and_mutagenesis_and_impairment_of_fertility"
L3_HEADER = ["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "rxcui",
             "unii", "pharm_class_epc", "pharm_class_moa", "rxcui_resolution",
             "n_chembl_ids_collapsed", "rra_score", "rra_p_value",
             "n_channels_supporting", "supporting_channels", "convergence",
             "rank_phenotype", "rank_prior", "rank_proximity"]


def _l3_table(results: Path, rows) -> None:
    directory = results / "l3"
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "candidates.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(L3_HEADER)
        for position, row in enumerate(rows, start=1):
            full = {"rank": position, "rxcui_resolution": "unii", "rra_score": "0.1",
                    "rra_p_value": "0.1", "n_chembl_ids_collapsed": "0",
                    "n_channels_supporting": "1", "supporting_channels": "proximity",
                    "convergence": "none", "rank_proximity": str(position / 10)}
            full.update(row)
            writer.writerow([full.get(c, "") for c in L3_HEADER])


def _config(tmp_path: Path) -> dict:
    return {"seed": 42, "results_dir": str(tmp_path / "results"),
            "therapeutic_endpoint": "chemoprevention"}


def _stub_sources(monkeypatch, labels: dict, marketing: dict) -> None:
    monkeypatch.setattr(labels_mod, "ensure_datasets",
                        lambda config, datasets=None: ({"label": [], "drugsfda": []}, []))
    monkeypatch.setattr(labels_mod, "load_labels",
                        lambda paths, rx, un: (labels, {"label_records": 0}))
    monkeypatch.setattr(labels_mod, "load_marketing",
                        lambda paths, rx, un: (marketing, {"drugsfda_records": 0}))


def _clean_label():
    return labels_mod.Label(sections={
        CARC: "Inventib was not genotoxic in any assay.",
        "pediatric_use": "Safety and effectiveness have been established in pediatric "
                         "patients 2 years of age and older."})


def _marketed():
    return labels_mod.Marketing(statuses=frozenset({"Prescription"}))


# ---------------------------------------------------------------------- the layer


def test_refuses_to_run_without_l3(tmp_path, monkeypatch):
    """An empty exclusions table would read as 'nothing was excluded'."""
    _stub_sources(monkeypatch, {}, {})
    with pytest.raises(FileNotFoundError, match="l3"):
        l4.run(_config(tmp_path))


def test_writes_survivors_and_exclusions_with_provenance(tmp_path, monkeypatch):
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]), [
        {"chembl_id": "CHEMBL9000001", "drug_name": "INVENTIB", "rxcui": "999001"},
        {"chembl_id": "CHEMBL9000002", "drug_name": "TOXICAINE", "rxcui": "999002"},
    ])
    _stub_sources(monkeypatch,
                  {"rxcui:999001": _clean_label(),
                   "rxcui:999002": labels_mod.Label(sections={
                       CARC: "Toxicaine was clastogenic in human lymphocytes.",
                       "pediatric_use": "Safety and effectiveness have been established "
                                        "in pediatric patients 2 years and older."})},
                  {"rxcui:999001": _marketed(), "rxcui:999002": _marketed()})
    l4.run(config)

    out = Path(config["results_dir"]) / "l4"
    survivors = list(csv.DictReader(open(out / "survivors.tsv"), delimiter="\t"))
    excluded = list(csv.DictReader(open(out / "excluded.tsv"), delimiter="\t"))
    assert [s["drug_name"] for s in survivors] == ["INVENTIB"]
    assert [e["drug_name"] for e in excluded] == ["TOXICAINE"]
    assert excluded[0]["excluded_by"] == "genotoxic_or_cancer_risk"
    assert excluded[0]["reason"] == "clastogenic_finding"
    assert "clastogenic" in excluded[0]["snippet"]
    assert excluded[0]["source_section"].startswith("13.1")


def test_every_candidate_appears_in_the_verdict_record(tmp_path, monkeypatch):
    """Survivors and excluded alike -- a filter whose decisions are invisible is unauditable."""
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]), [
        {"chembl_id": "C1", "drug_name": "INVENTIB", "rxcui": "999001"},
        {"chembl_id": "C2", "drug_name": "NOLABEL", "rxcui": "999009"},
    ])
    _stub_sources(monkeypatch, {"rxcui:999001": _clean_label()},
                  {"rxcui:999001": _marketed()})
    l4.run(config)
    record = json.loads((Path(config["results_dir"]) / "l4" / "verdicts.json")
                        .read_text(encoding="utf-8"))
    assert len(record["candidates"]) == 2
    for candidate in record["candidates"]:
        assert {v["rule"] for v in candidate["verdicts"]}
        assert isinstance(candidate["survived"], bool)


def test_a_candidate_with_no_label_is_excluded_not_passed(tmp_path, monkeypatch):
    """Fail closed: absence of a label is not evidence of safety."""
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]),
              [{"chembl_id": "C1", "drug_name": "NOLABEL", "rxcui": "999009"}])
    _stub_sources(monkeypatch, {}, {})
    with pytest.raises(ValueError, match="every one of the"):
        l4.run(config)
    excluded = list(csv.DictReader(
        open(Path(config["results_dir"]) / "l4" / "excluded.tsv"), delimiter="\t"))
    assert excluded[0]["reason"] == "insufficient_evidence"


def test_excluding_everything_raises_rather_than_shipping_an_empty_table(tmp_path,
                                                                         monkeypatch):
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]),
              [{"chembl_id": "C1", "drug_name": "NOLABEL", "rxcui": "999009"}])
    _stub_sources(monkeypatch, {}, {})
    with pytest.raises(ValueError, match="reportable result"):
        l4.run(config)


def test_the_record_carries_counts_and_caveats(tmp_path, monkeypatch):
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]),
              [{"chembl_id": "C1", "drug_name": "INVENTIB", "rxcui": "999001"}])
    _stub_sources(monkeypatch, {"rxcui:999001": _clean_label()},
                  {"rxcui:999001": _marketed()})
    l4.run(config)
    meta = json.loads((Path(config["results_dir"]) / "l4" / "validation.json")
                      .read_text(encoding="utf-8"))
    assert meta["counts"]["survivors"] == 1
    assert "genotoxic_or_cancer_risk" in meta["method"]["hard_rules"]
    assert any("not a safety claim" in c for c in meta["caveats"])
    assert any("not assessed" in c for c in meta["caveats"]), "DDI limitation must be stated"


# ------------------------------------------------------------------- the benchmark


def test_the_circular_channel_is_excluded_from_the_headline(tmp_path, monkeypatch):
    """Channel E is built from the positive set; scoring it on that set measures nothing."""
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]),
              [{"chembl_id": "C1", "drug_name": "INVENTIB", "rank_prior": "0.5"}])
    monkeypatch.setattr(benchmark, "positive_set", lambda path=None: (["inventib"], "abc"))
    report = benchmark.run_blinded(config)

    assert "prior" in report["circular_channels"]
    assert "prior" not in report["headline_channels"]
    assert report["per_channel"]["prior"]["circular"] is True
    assert "measures nothing" in report["per_channel"]["prior"]["note"]
    assert report["combined"]["partly_circular"] is True


def test_the_positive_set_is_fingerprinted():
    """A benchmark whose contents can drift unnoticed is not blinded."""
    names, digest = benchmark.positive_set()
    assert names and len(digest) == 64
    again = benchmark.positive_set()[1]
    assert digest == again


def test_auroc_is_chance_when_positives_are_spread_evenly():
    positions = [2, 4, 6, 8]
    value = benchmark.auroc(positions, n_ranked=8, n_positives=4)
    assert 0.3 < value < 0.7


def test_auroc_is_one_when_positives_lead_the_ranking():
    assert benchmark.auroc([1, 2, 3], n_ranked=10, n_positives=3) == pytest.approx(1.0)


def test_auroc_is_zero_when_positives_trail_the_ranking():
    assert benchmark.auroc([8, 9, 10], n_ranked=10, n_positives=3) == pytest.approx(0.0)


def test_auroc_is_undefined_without_negatives():
    value = benchmark.auroc([1, 2], n_ranked=2, n_positives=2)
    assert value != value  # NaN


def test_the_bootstrap_interval_is_reported_for_a_small_set():
    ci = benchmark.bootstrap_auroc([1, 5, 9], n_ranked=20, n_positives=3, seed=42,
                                   samples=200)
    assert ci["ci_low"] is not None and ci["ci_low"] <= ci["ci_high"]


def test_the_bootstrap_declines_to_resample_a_single_positive():
    ci = benchmark.bootstrap_auroc([1], n_ranked=10, n_positives=1, seed=42)
    assert ci["ci_low"] is None and "too few" in ci["note"]


def test_the_benchmark_records_that_it_is_not_an_efficacy_claim(tmp_path, monkeypatch):
    config = _config(tmp_path)
    _l3_table(Path(config["results_dir"]), [{"chembl_id": "C1", "drug_name": "INVENTIB"}])
    monkeypatch.setattr(benchmark, "positive_set", lambda path=None: (["inventib"], "abc"))
    report = benchmark.run_blinded(config)
    assert any("not clinical efficacy" in c for c in report["caveats"])
    assert any("re-freezing" in c for c in report["caveats"])
