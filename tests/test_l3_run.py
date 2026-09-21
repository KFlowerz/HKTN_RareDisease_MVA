"""Tests for the L3 layer: what it reads, what it refuses, and what it writes.

Invented ChEMBL ids throughout. The harmonisation and crosswalk downloads are stubbed, so
nothing here reaches the network or the enrichment cache.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

# The package re-exports its own run() function, so both `from src.l3_integrate import
# run` and `import src.l3_integrate.run as l3` hand back the function rather than the
# module. import_module reads sys.modules and is unambiguous.
import importlib

from src.l3_integrate import harmonize

l3 = importlib.import_module("src.l3_integrate.run")


def _channel_table(results: Path, channel: str, chembl_ids) -> None:
    directory = results / "l2" / channel
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "candidates.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("rank\tchembl_id\tdrug_name\n")
        for position, chembl_id in enumerate(chembl_ids, start=1):
            handle.write(f"{position}\t{chembl_id}\tINVENTED {chembl_id}\n")


def _status(results: Path, complete: dict, not_implemented=()) -> None:
    (results / "l2").mkdir(parents=True, exist_ok=True)
    channels = {name: {"module": f"src.l2_channels.channel_{name}", "status": "complete",
                       "output_dir": f"l2/{name}", "n_candidates": len(ids)}
                for name, ids in complete.items()}
    channels.update({name: {"module": "x", "status": "not_implemented"}
                     for name in not_implemented})
    (results / "l2" / "channels.json").write_text(json.dumps({
        "finished_utc": "2026-09-18T00:00:00+00:00", "seed": 42, "causal_gene": "INVGENE",
        "channels": channels,
        "summary": {"complete": sorted(complete), "failed": [],
                    "not_implemented": sorted(not_implemented), "disabled": []},
        "n_complete": len(complete),
    }, indent=2), encoding="utf-8")


def _config(tmp_path: Path) -> dict:
    # Reasoning off: these tests cover harmonisation and rank aggregation, and should not
    # depend on a model server being up. The reasoning step has its own tests, and L3
    # records its absence in integration.json rather than failing -- which is asserted in
    # test_field_contract_reports_reasoning_absence below.
    return {"seed": 42, "results_dir": str(tmp_path / "results"),
            "causal_gene": "INVGENE", "therapeutic_endpoint": "chemoprevention",
            "l3": {"reasoning": False}}


def _stub_harmonize(monkeypatch) -> None:
    """Identity harmonisation that maps each id to itself, so L3's own logic is tested."""
    def build(config, chembl_ids):
        identities = {cid: harmonize.Identity(key=cid, name=f"INVENTED {cid}",
                                              chembl_ids=frozenset({cid}),
                                              rxcui=(), resolution="unresolved")
                      for cid in chembl_ids}
        return identities, {"unii": 0, "name": 0, "unresolved": len(identities)}, []

    monkeypatch.setattr(l3.harmonize, "build", build)


def _setup(tmp_path, monkeypatch, complete, not_implemented=()):
    results = Path(tmp_path) / "results"
    for name, ids in complete.items():
        _channel_table(results, name, ids)
    _status(results, complete, not_implemented)
    _stub_harmonize(monkeypatch)
    return _config(tmp_path)


# ----------------------------------------------------------- channels.json rules


def test_refuses_to_run_without_an_l2_status_record(tmp_path, monkeypatch):
    """Inferring channel status from whatever is on disk is what the record prevents."""
    _stub_harmonize(monkeypatch)
    with pytest.raises(FileNotFoundError, match="no record of which channels"):
        l3.run(_config(tmp_path))


def test_reads_only_the_channels_recorded_complete(tmp_path, monkeypatch):
    """A failed channel's table on disk must not contribute silently."""
    results = Path(tmp_path) / "results"
    _channel_table(results, "alpha", ["C1", "C2"])
    _channel_table(results, "beta", ["C2", "C3"])
    _channel_table(results, "ghost", ["C9"])      # on disk, absent from the record
    _status(results, {"alpha": ["C1", "C2"], "beta": ["C2", "C3"]})
    _stub_harmonize(monkeypatch)

    config = _config(tmp_path)
    l3.run(config)
    rows = list(csv.DictReader(
        open(Path(config["results_dir"]) / "l3" / "candidates.tsv"), delimiter="\t"))
    assert {r["chembl_id"] for r in rows} == {"C1", "C2", "C3"}
    assert "C9" not in {r["chembl_id"] for r in rows}


def test_a_complete_channel_with_a_missing_table_is_an_error(tmp_path, monkeypatch):
    """The record and the directory disagreeing means a stale state, not a smaller run."""
    results = Path(tmp_path) / "results"
    _channel_table(results, "alpha", ["C1"])
    _status(results, {"alpha": ["C1"], "beta": ["C2"]})
    _stub_harmonize(monkeypatch)
    with pytest.raises(ValueError, match="recorded complete"):
        l3.run(_config(tmp_path))


def test_a_table_without_an_identity_column_is_refused(tmp_path, monkeypatch):
    results = Path(tmp_path) / "results"
    directory = results / "l2" / "alpha"
    directory.mkdir(parents=True)
    (directory / "candidates.tsv").write_text("rank\tdrug_name\n1\tINVENTIB\n",
                                              encoding="utf-8")
    _channel_table(results, "beta", ["C1"])
    _status(results, {"alpha": ["C1"], "beta": ["C1"]})
    _stub_harmonize(monkeypatch)
    with pytest.raises(ValueError, match="chembl_id"):
        l3.run(_config(tmp_path))


def test_one_channel_is_not_aggregation(tmp_path, monkeypatch):
    """Rank aggregation over one list is that list with extra steps -- say so."""
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1", "C2"]})
    with pytest.raises(ValueError, match="nothing to aggregate"):
        l3.run(config)


# ------------------------------------------------------------------- the output


def test_writes_a_ranked_table_and_a_record(tmp_path, monkeypatch):
    config = _setup(tmp_path, monkeypatch,
                    {"alpha": ["C1", "C2"], "beta": ["C2", "C3"]})
    l3.run(config)

    out = Path(config["results_dir"]) / "l3"
    rows = list(csv.DictReader(open(out / "candidates.tsv"), delimiter="\t"))
    assert [r["chembl_id"] for r in rows][0] == "C2", "the only drug both channels ranked"
    assert rows[0]["n_channels_supporting"] == "2"
    assert set(rows[0]["supporting_channels"].split("|")) == {"alpha", "beta"}

    meta = json.loads((out / "integration.json").read_text(encoding="utf-8"))
    assert meta["channels_used"]["complete"] == ["alpha", "beta"]
    assert meta["seed"] == 42 and meta["causal_gene"] == "INVGENE"


def test_each_channel_gets_its_own_rank_column(tmp_path, monkeypatch):
    """D2's contract: one column per channel, and absent is blank rather than zero."""
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1"], "beta": ["C2"]})
    l3.run(config)
    rows = list(csv.DictReader(
        open(Path(config["results_dir"]) / "l3" / "candidates.tsv"), delimiter="\t"))
    by_id = {r["chembl_id"]: r for r in rows}
    assert by_id["C1"]["rank_alpha"] and by_id["C1"]["rank_beta"] == ""
    assert by_id["C2"]["rank_beta"] and by_id["C2"]["rank_alpha"] == ""


def test_the_record_states_which_fields_are_missing(tmp_path, monkeypatch):
    """L5 must not discover missing reasoning fields at render time (D2)."""
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1"], "beta": ["C1"]})
    l3.run(config)
    meta = json.loads((Path(config["results_dir"]) / "l3" / "integration.json")
                      .read_text(encoding="utf-8"))
    outstanding = set(meta["field_contract"]["outstanding"])
    assert {"rationale", "contradicting_evidence", "confidence"} <= outstanding
    assert meta["field_contract"]["reasoning"] == "disabled"


def test_field_contract_reports_reasoning_absence(tmp_path, monkeypatch):
    """With reasoning on but no server, L3 still completes and says what is missing.

    The ranking costs minutes of harmonisation and aggregation. Losing it because a
    model server was not running would be the wrong trade, and silently omitting the
    reasoning fields would leave L5 to discover the gap at render time. So the layer
    succeeds, and the gap is a recorded fact.
    """
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1"], "beta": ["C1"]})
    # Point at a loopback port with nothing on it, and keep the probe short: a dead
    # loopback port does not refuse under WSL mirrored networking, it hangs.
    config["reasoning_endpoint"] = "http://127.0.0.1:8099"
    config["l3"] = {"reasoning": True, "reasoning_probe_seconds": 1}

    l3.run(config)

    out = Path(config["results_dir"]) / "l3"
    assert (out / "candidates.tsv").exists(), "the ranking must survive a missing server"
    meta = json.loads((out / "integration.json").read_text(encoding="utf-8"))
    assert meta["field_contract"]["reasoning"].startswith("unavailable")
    assert "rationale" in meta["field_contract"]["outstanding"]
    assert not (out / "reasoning" / "rationales.json").exists()


def test_the_caveats_travel_with_the_result(tmp_path, monkeypatch):
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1"], "beta": ["C1"]})
    l3.run(config)
    meta = json.loads((Path(config["results_dir"]) / "l3" / "integration.json")
                      .read_text(encoding="utf-8"))
    assert any("does not judge" in c for c in meta["caveats"])
    assert any("not fully independent" in c for c in meta["caveats"])


def test_no_patient_derived_file_is_read(tmp_path, monkeypatch):
    """Channel D's evidence.json states which diseases resemble the subject."""
    config = _setup(tmp_path, monkeypatch, {"alpha": ["C1"], "beta": ["C1"]})
    marker = Path(config["results_dir"]) / "l2" / "alpha" / "evidence.json"
    marker.write_text('{"patient_derived": true}', encoding="utf-8")

    import builtins

    opened = []
    real_open = builtins.open

    def watched(path, *args, **kwargs):
        opened.append(str(path))
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr("builtins.open", watched)
    l3.run(config)
    assert not any("evidence.json" in p for p in opened)
