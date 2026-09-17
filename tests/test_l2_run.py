"""Tests for the L2 channel runner.

Channels are replaced with invented stand-ins; no real channel, drug or data is used.
"""

from __future__ import annotations

import importlib
import json
import random
import types
from pathlib import Path

import pytest

# import_module, not `from src.l2_channels import run` -- the package re-exports the
# function under that name, so the latter binds the function, not the module.
l2 = importlib.import_module("src.l2_channels.run")

ALL_OFF = {key: False for key in l2.CHANNEL_REGISTRY}


def _config(tmp_path, **channels) -> dict:
    return {"results_dir": tmp_path, "seed": 42, "causal_gene": "GENE1",
            "channels": {**ALL_OFF, **channels}}


def _fake(monkeypatch, key, generate, name=None):
    """Swap a registry entry for a stand-in module whose ``generate`` is given."""
    module = types.SimpleNamespace(__name__=f"src.l2_channels.{name or 'channel_' + key}",
                                   generate=generate)
    monkeypatch.setitem(l2.CHANNEL_REGISTRY, key, module)
    return module


def _writes(rows, key):
    def generate(config):
        out = l2.channel_dir(config, key)
        out.mkdir(parents=True, exist_ok=True)
        (out / l2.CANDIDATES_FILE).write_text("rank\tdrug\n" + "".join(
            f"{i}\tD{i}\n" for i in range(1, rows + 1)), encoding="utf-8")
    return generate


def _status(tmp_path) -> dict:
    return json.loads((tmp_path / "l2" / l2.STATUS_FILE).read_text(encoding="utf-8"))


# ------------------------------------------------------------ config checks


def test_an_unknown_channel_key_is_refused(tmp_path) -> None:
    """A typo must not silently disable a line of evidence."""
    config = _config(tmp_path)
    config["channels"]["proximty"] = True
    with pytest.raises(ValueError, match="unknown \\['proximty'\\]"):
        l2.run(config)


def test_a_missing_channel_key_is_refused(tmp_path) -> None:
    config = _config(tmp_path)
    del config["channels"]["prior"]
    with pytest.raises(ValueError, match="missing \\['prior'\\]"):
        l2.run(config)


def test_a_quoted_false_is_refused(tmp_path) -> None:
    """``"false"`` in YAML quotes is a truthy string."""
    with pytest.raises(ValueError, match="true or false"):
        l2.run(_config(tmp_path, kg="false"))


def test_enabled_channels_keep_registry_order(tmp_path) -> None:
    assert l2.enabled_channels(_config(tmp_path, prior=True, kg=True)) == ["kg", "prior"]


# ------------------------------------------------------------------ outcomes


def test_a_complete_channel_is_recorded_with_its_table(tmp_path, monkeypatch) -> None:
    _fake(monkeypatch, "proximity", _writes(3, "proximity"))
    l2.run(_config(tmp_path, proximity=True))

    status = _status(tmp_path)
    record = status["channels"]["proximity"]
    assert record["status"] == "complete"
    assert record["n_candidates"] == 3
    assert record["artifacts"] == ["l2/channel_proximity/candidates.tsv"]
    assert status["n_complete"] == 1
    assert status["channels"]["kg"]["status"] == "disabled"
    assert status["seed"] == 42 and status["causal_gene"] == "GENE1"


def test_a_failing_channel_does_not_stop_the_others(tmp_path, monkeypatch) -> None:
    def boom(config):
        raise OSError("network unreachable")

    _fake(monkeypatch, "kg", boom)
    _fake(monkeypatch, "proximity", _writes(2, "proximity"))
    with pytest.raises(RuntimeError, match="failed: kg"):
        l2.run(_config(tmp_path, kg=True, proximity=True))

    status = _status(tmp_path)
    assert status["channels"]["kg"]["status"] == "failed"
    assert "OSError: network unreachable" in status["channels"]["kg"]["detail"]
    assert status["channels"]["proximity"]["status"] == "complete"


def test_an_unimplemented_channel_makes_the_layer_not_implemented(tmp_path, monkeypatch) -> None:
    def stub(config):
        raise NotImplementedError("scaffold")

    _fake(monkeypatch, "signature", stub)
    _fake(monkeypatch, "proximity", _writes(1, "proximity"))
    with pytest.raises(NotImplementedError, match="signature.*Completed: proximity"):
        l2.run(_config(tmp_path, signature=True, proximity=True))
    assert _status(tmp_path)["summary"]["not_implemented"] == ["signature"]


def test_a_failure_outranks_a_stub(tmp_path, monkeypatch) -> None:
    """A real failure must not be reported as the milder "not implemented"."""
    def stub(config):
        raise NotImplementedError("scaffold")

    def boom(config):
        raise ValueError("broken crosswalk")

    _fake(monkeypatch, "signature", stub)
    _fake(monkeypatch, "kg", boom)
    with pytest.raises(RuntimeError):
        l2.run(_config(tmp_path, signature=True, kg=True))


@pytest.mark.parametrize("rows", [None, 0])
def test_returning_without_candidates_is_a_failure(tmp_path, monkeypatch, rows) -> None:
    """An empty list is a scientific claim; silence from a channel is not one."""
    generate = (lambda config: None) if rows is None else _writes(rows, "proximity")
    _fake(monkeypatch, "proximity", generate)
    with pytest.raises(RuntimeError):
        l2.run(_config(tmp_path, proximity=True))
    record = _status(tmp_path)["channels"]["proximity"]
    assert record["status"] == "failed"
    assert "non-empty candidates.tsv" in record["detail"]


def test_a_stale_table_cannot_survive_a_failure(tmp_path, monkeypatch) -> None:
    """Last run's table on disk must not be read as this run's evidence."""
    _fake(monkeypatch, "proximity", _writes(5, "proximity"))
    l2.run(_config(tmp_path, proximity=True))
    stale = tmp_path / "l2" / "channel_proximity" / l2.CANDIDATES_FILE
    assert stale.exists()

    def boom(config):
        raise OSError("down")

    _fake(monkeypatch, "proximity", boom)
    with pytest.raises(RuntimeError):
        l2.run(_config(tmp_path, proximity=True))
    assert not stale.exists()
    assert _status(tmp_path)["n_complete"] == 0


def test_a_stale_status_record_cannot_outlive_a_bad_config(tmp_path, monkeypatch) -> None:
    """L3 reads channels.json; last run's copy must not stand in for a run that died."""
    _fake(monkeypatch, "proximity", _writes(1, "proximity"))
    l2.run(_config(tmp_path, proximity=True))
    assert (tmp_path / "l2" / l2.STATUS_FILE).exists()

    config = _config(tmp_path, proximity=True)
    config["channels"]["typo"] = True
    with pytest.raises(ValueError):
        l2.run(config)
    assert not (tmp_path / "l2" / l2.STATUS_FILE).exists()


def test_a_channel_failure_message_is_not_logged(tmp_path, monkeypatch, caplog) -> None:
    """A channel reading patient-derived input may raise with that content in its message."""
    def boom(config):
        raise ValueError("INVENTED-SENSITIVE-DETAIL")

    _fake(monkeypatch, "phenotype", boom)
    with caplog.at_level("DEBUG"), pytest.raises(RuntimeError):
        l2.run(_config(tmp_path, phenotype=True))
    assert "INVENTED-SENSITIVE-DETAIL" not in caplog.text
    assert "INVENTED-SENSITIVE-DETAIL" in _status(tmp_path)["channels"]["phenotype"]["detail"]


# ------------------------------------------------------------- independence


def test_a_channel_cannot_alter_the_config_another_sees(tmp_path, monkeypatch) -> None:
    seen = {}

    def meddler(config):
        config["seed"] = 7
        config["channels"]["prior"] = False
        _writes(1, "kg")(config)

    def observer(config):
        seen.update(seed=config["seed"], prior=config["channels"]["prior"])
        _writes(1, "prior")(config)

    _fake(monkeypatch, "kg", meddler)
    _fake(monkeypatch, "prior", observer)
    l2.run(_config(tmp_path, kg=True, prior=True))
    assert seen == {"seed": 42, "prior": True}


def test_rng_state_does_not_leak_between_channels(tmp_path, monkeypatch) -> None:
    """A channel's draws must not depend on whether another channel ran first."""
    draws = {}

    def consumer(config):
        [random.random() for _ in range(100)]
        _writes(1, "kg")(config)

    def sampler(config):
        draws.setdefault("run", []).append(random.random())
        _writes(1, "prior")(config)

    _fake(monkeypatch, "kg", consumer)
    _fake(monkeypatch, "prior", sampler)
    l2.run(_config(tmp_path, prior=True))
    l2.run(_config(tmp_path, kg=True, prior=True))
    assert draws["run"][0] == draws["run"][1]


# ------------------------------------------------------ against the orchestrator


def test_the_orchestrator_lists_l2_artifacts(tmp_path, monkeypatch) -> None:
    """L2 writes results/l2/, not results/l2_channels/; the manifest must still see it."""
    from src import pipeline

    _fake(monkeypatch, "proximity", _writes(2, "proximity"))
    config = _config(tmp_path, proximity=True)
    pipeline.run(config, only="l2_channels")

    manifest = json.loads((tmp_path / pipeline.MANIFEST_NAME).read_text(encoding="utf-8"))
    layer = manifest["layers"]["l2_channels"]
    assert layer["status"] == "complete"
    assert "l2/channels.json" in [a.replace("\\", "/") for a in layer["artifacts"]]


def test_the_configured_channels_are_all_valid() -> None:
    import yaml

    config = yaml.safe_load((Path(__file__).resolve().parent.parent / "config" / "pipeline.yaml")
                            .read_text(encoding="utf-8"))
    l2.enabled_channels(config)


def test_a_disabled_channels_old_table_is_not_listed_as_this_runs_output(tmp_path,
                                                                        monkeypatch) -> None:
    """channels.json is what L3 reads; the manifest must not contradict it.

    A channel switched off after an earlier run keeps its table on disk. Listing that as an
    artifact of this run would present last week's candidates as today's evidence.
    """
    from src import pipeline

    _fake(monkeypatch, "proximity", _writes(2, "proximity"))
    _fake(monkeypatch, "prior", _writes(2, "prior"))
    pipeline.run(_config(tmp_path, proximity=True, prior=True), only="l2_channels")
    assert (tmp_path / "l2" / "channel_prior" / l2.CANDIDATES_FILE).exists()

    pipeline.run(_config(tmp_path, proximity=True), only="l2_channels")
    manifest = json.loads((tmp_path / pipeline.MANIFEST_NAME).read_text(encoding="utf-8"))
    listed = [a.replace("\\", "/") for a in manifest["layers"]["l2_channels"]["artifacts"]]

    assert "l2/channel_proximity/candidates.tsv" in listed
    assert "l2/channel_prior/candidates.tsv" not in listed, "a disabled channel's stale table"
    assert "l2/channels.json" in listed
    # The file itself is left alone: only channels that ran get their directory cleared.
    assert (tmp_path / "l2" / "channel_prior" / l2.CANDIDATES_FILE).exists()


def test_an_unreadable_status_record_lists_only_itself(tmp_path, monkeypatch) -> None:
    from src import pipeline

    _fake(monkeypatch, "proximity", _writes(1, "proximity"))
    pipeline.run(_config(tmp_path, proximity=True), only="l2_channels")
    (tmp_path / "l2" / l2.STATUS_FILE).write_text("{not json", encoding="utf-8")
    assert pipeline._artifacts(tmp_path, "l2_channels") == ["l2/channels.json".replace("/", __import__("os").sep)]
