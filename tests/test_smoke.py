"""Smoke tests: every layer imports, every ``run()`` is a callable stub.

These do not test behavior -- there is none yet. They guard the scaffold's two
invariants: the package tree is importable, and the config's two open decisions
(``causal_gene``, ``therapeutic_endpoint``) have not been quietly filled in.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

LAYER_MODULES = [
    "src.l0_genomics.run",
    "src.l1_target.run",
    "src.l2_channels.run",
    "src.l3_integrate.run",
    "src.l4_validate.run",
    "src.l5_report.run",
]

CHANNEL_MODULES = [
    "src.l2_channels.channel_a_kg",
    "src.l2_channels.channel_b_proximity",
    "src.l2_channels.channel_c_signature",
    "src.l2_channels.channel_d_phenotype",
    "src.l2_channels.channel_e_prior",
]

HELPER_MODULES = [
    "src.l3_integrate.harmonize",
    "src.l3_integrate.aggregate",
    "src.l3_integrate.claude_reasoning",
    "src.l4_validate.safety_triage",
    "src.l4_validate.benchmark",
]


@pytest.mark.parametrize("module_name", LAYER_MODULES)
def test_layer_run_is_a_stub(module_name: str) -> None:
    """Each layer exposes a callable ``run`` that raises ``NotImplementedError``."""
    module = importlib.import_module(module_name)
    assert callable(module.run)
    with pytest.raises(NotImplementedError):
        module.run({})


@pytest.mark.parametrize("module_name", CHANNEL_MODULES)
def test_channel_generate_is_a_stub(module_name: str) -> None:
    """Each L2 channel exposes a callable ``generate`` that raises."""
    module = importlib.import_module(module_name)
    assert callable(module.generate)
    with pytest.raises(NotImplementedError):
        module.generate({})


@pytest.mark.parametrize("module_name", HELPER_MODULES)
def test_helper_modules_import(module_name: str) -> None:
    """L3/L4 helper modules import cleanly."""
    assert importlib.import_module(module_name) is not None


def test_pipeline_imports_and_orchestrates() -> None:
    """The orchestrator imports and knows the layer order."""
    from src import pipeline

    assert pipeline.LAYER_ORDER == (
        "l0_genomics",
        "l1_target",
        "l2_channels",
        "l3_integrate",
        "l4_validate",
        "l5_report",
    )
    for name in ("load_config", "seed_everything", "run", "main"):
        assert callable(getattr(pipeline, name))


def test_channel_registry_covers_config_channels() -> None:
    """Every channel toggle in the config maps to a registered channel module."""
    from src.l2_channels import CHANNEL_REGISTRY

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    assert set(config["channels"]) == set(CHANNEL_REGISTRY)


def test_open_decisions_are_still_open() -> None:
    """``causal_gene`` and ``therapeutic_endpoint`` must not be hardcoded.

    The causal gene is a finding from L0 / Track-1 reconciliation, and the therapeutic
    endpoint is the project's primary open scientific decision. A value here means
    someone guessed.
    """
    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    assert config["causal_gene"] is None
    assert config["therapeutic_endpoint"] is None
    assert config["seed"] == 42


def test_no_patient_data_files_present() -> None:
    """No genomic file may exist anywhere in the repo -- including outside data/."""
    patterns = ("*.vcf", "*.vcf.gz", "*.bam", "*.cram", "*.fastq", "*.fastq.gz", "*.fq")
    offenders = [
        path
        for pattern in patterns
        for path in REPO_ROOT.rglob(pattern)
        if ".git" not in path.parts
    ]
    assert offenders == [], f"genomic files found in repo: {offenders}"
