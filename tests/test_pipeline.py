"""Tests for the orchestrator.

These cover the failure modes that would otherwise be silent: a config missing a key,
a results directory nested inside the data directory, the environment override not
taking effect, and the two open decisions being quietly filled in.

No patient data is touched. Every test builds its own directories under ``tmp_path``.
"""

from __future__ import annotations

import json
import logging
import random
from pathlib import Path

import pytest
import yaml

from src import pipeline


def _write_config(tmp_path: Path, **overrides) -> Path:
    """Write a minimal valid config and return its path."""
    body = {
        "seed": 42,
        "data_dir": "./data",
        # An absolute tmp path, NOT "./results": load_config resolves a relative value
        # against the repository root, so the suite was writing into -- and reading from --
        # the real results directory. Once L4 was implemented that made this test run the
        # layer against production artifacts and the multi-gigabyte openFDA cache, taking
        # 75 seconds and reporting "complete" for a layer it meant to find unimplemented.
        "results_dir": str(tmp_path / "results"),
        "causal_gene": None,
        "therapeutic_endpoint": None,
        "model": "claude-opus-5",
        "channels": {"kg": True, "proximity": True, "signature": True, "phenotype": True, "prior": True},
    }
    body.update(overrides)
    path = tmp_path / "pipeline.yaml"
    path.write_text(yaml.safe_dump(body), encoding="utf-8")
    return path


@pytest.fixture()
def data_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A custody-root-shaped directory with a raw/ subdir, exported as MVA_DATA_ROOT."""
    root = tmp_path / "mva-data"
    (root / "raw").mkdir(parents=True)
    monkeypatch.setenv("MVA_DATA_ROOT", str(root))
    return root


def test_missing_required_key_is_rejected(tmp_path: Path, data_root: Path) -> None:
    """A config without `seed` must fail loudly, not default to something."""
    body = yaml.safe_load(_write_config(tmp_path).read_text(encoding="utf-8"))
    del body["seed"]
    path = tmp_path / "broken.yaml"
    path.write_text(yaml.safe_dump(body), encoding="utf-8")

    with pytest.raises(ValueError, match="missing required key"):
        pipeline.load_config(path)


def test_data_root_env_overrides_config(tmp_path: Path, data_root: Path) -> None:
    """MVA_DATA_ROOT wins over the committed data_dir, and resolves to raw/."""
    config = pipeline.load_config(_write_config(tmp_path))
    assert config["data_dir"] == data_root / "raw"


def test_missing_data_dir_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Pointing at a data directory that does not exist fails before any layer runs."""
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "nope"))
    with pytest.raises(FileNotFoundError, match="data_dir does not exist"):
        pipeline.load_config(_write_config(tmp_path))


def test_results_dir_may_not_nest_inside_data_dir(
    tmp_path: Path, data_root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Derived artifacts must not be writable into the source-data tree.

    Nesting them would make a derived artifact indistinguishable from source data, and
    could hide it from a purge scoped to one tree (docs/data_custody.md).
    """
    monkeypatch.setattr(pipeline, "REPO_ROOT", data_root / "raw")
    with pytest.raises(ValueError, match="must not sit inside"):
        pipeline.load_config(_write_config(tmp_path, results_dir="./sub"))


def test_open_decisions_warn_but_do_not_default(
    tmp_path: Path, data_root: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Null causal_gene / therapeutic_endpoint must warn and stay null.

    Guessing either would silently fabricate the scientific premise of the run.
    """
    with caplog.at_level(logging.WARNING):
        config = pipeline.load_config(_write_config(tmp_path))

    assert config["causal_gene"] is None
    assert config["therapeutic_endpoint"] is None
    assert "causal_gene is null" in caplog.text
    assert "therapeutic_endpoint is null" in caplog.text


def test_seed_everything_is_reproducible() -> None:
    """The same seed must produce the same stream."""
    pipeline.seed_everything(42)
    first = [random.random() for _ in range(5)]
    pipeline.seed_everything(42)
    assert [random.random() for _ in range(5)] == first


def test_unknown_layer_is_rejected(tmp_path: Path, data_root: Path) -> None:
    """A typo in --only must fail rather than silently running nothing."""
    config = pipeline.load_config(_write_config(tmp_path))
    with pytest.raises(ValueError, match="unknown layer"):
        pipeline.run(config, only="l9_nonexistent")


def test_stub_layer_records_not_implemented_and_stops(tmp_path: Path, data_root: Path) -> None:
    """A scaffold stub is recorded as `not_implemented`, not as a completed layer.

    The distinction matters: a layer that ran and produced nothing is a scientific claim;
    a layer that was never implemented is not.

    Targets L5 explicitly. Every earlier layer is implemented and fails on its missing
    input instead -- a real VCF, a causal gene, a channel map, an L2 status record, an L3
    table -- so this asserts the manifest behaviour on the one layer that genuinely is a
    scaffold. It has moved down the stack as layers were built, which is the intent.
    """
    config = pipeline.load_config(_write_config(tmp_path))
    pipeline.run(config, only="l5_report")

    manifest = json.loads((config["results_dir"] / pipeline.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["layers"]["l5_report"]["status"] == "not_implemented"
    assert manifest["seed"] == 42
    assert manifest["causal_gene"] is None


def test_layer_failure_is_recorded_before_it_propagates(tmp_path: Path, data_root: Path) -> None:
    """A layer that raises must leave a manifest entry, so a re-run can resume.

    L0 with an empty data directory raises FileNotFoundError -- a real failure, not an
    unimplemented body. The manifest must capture it rather than losing the run's state.
    """
    config = pipeline.load_config(_write_config(tmp_path))
    with pytest.raises(FileNotFoundError):
        pipeline.run(config, only="l0_genomics")

    manifest = json.loads((config["results_dir"] / pipeline.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["layers"]["l0_genomics"]["status"] == "failed"
    assert "FileNotFoundError" in manifest["layers"]["l0_genomics"]["detail"]


def test_manifest_records_model_choice(tmp_path: Path, data_root: Path) -> None:
    """The model that ran is resolvable from config, not hardcoded in source."""
    from src.l3_integrate.claude_reasoning import DEFAULT_MODEL, resolve_model

    config = pipeline.load_config(_write_config(tmp_path, model="llama-3.1-8b-instruct-q4_k_m"))
    assert resolve_model(config) == "llama-3.1-8b-instruct-q4_k_m"
    assert resolve_model({}) == DEFAULT_MODEL
