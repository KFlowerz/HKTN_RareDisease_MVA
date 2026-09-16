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


#: Layers still unimplemented end to end. L0 is absent: it is implemented, and needs a
#: real VCF and the local snpEff database, so calling it with an empty config raises
#: KeyError on the missing data_dir. Layers move off this list as they are built, one at
#: a time and deliberately.
STUB_LAYER_MODULES = [m for m in LAYER_MODULES if m != "src.l0_genomics.run"]


@pytest.mark.parametrize("module_name", STUB_LAYER_MODULES)
def test_layer_run_is_a_stub(module_name: str) -> None:
    """Each unimplemented layer exposes a callable ``run`` that raises."""
    module = importlib.import_module(module_name)
    assert callable(module.run)
    with pytest.raises(NotImplementedError):
        module.run({})


def test_l0_is_partially_implemented() -> None:
    """L0 requires real config: it is no longer a stub that ignores its input.

    The burden half reads a VCF from ``data_dir``, so an empty config fails on the
    missing key rather than on an unimplemented body. Asserting this keeps the smoke
    suite honest about which layers are built.
    """
    # import_module, not `from src.l0_genomics import run` -- the package's __init__
    # re-exports the function under that name, so the latter binds the function.
    l0 = importlib.import_module("src.l0_genomics.run")

    assert callable(l0.run)
    with pytest.raises(KeyError):
        l0.run({})


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


def test_causal_gene_is_still_a_finding() -> None:
    """``causal_gene`` must not be hardcoded.

    It is a finding from L0, not a setting. A value here before L0 has run means
    someone guessed, which would silently fabricate the scientific premise of the run.
    When L0 does produce a call, this test changes to assert the config matches L0's
    output artifact -- it does not get deleted.
    """
    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    assert config["causal_gene"] is None
    assert config["seed"] == 42


def test_therapeutic_endpoint_is_decided_and_documented() -> None:
    """``therapeutic_endpoint`` holds a valid value, recorded in the decision log.

    This replaces an assertion that it was ``None``. Gate G2 passed on 2026-09-08, so
    "still null" is no longer the invariant -- but the reason the original test existed
    has not gone away, so it is tightened rather than removed: the value must be one of
    the three admissible endpoints, and it must be traceable to a written decision.
    An endpoint nobody recorded is indistinguishable from an endpoint someone guessed.
    """
    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    endpoint = config["therapeutic_endpoint"]

    assert endpoint in {"chemoprevention", "symptomatic", "mitotic_fidelity"}, (
        f"{endpoint!r} is not an admissible therapeutic endpoint"
    )

    decisions = (REPO_ROOT / "mngmt" / "decisions.md").read_text(encoding="utf-8")
    assert endpoint in decisions, (
        f"therapeutic_endpoint is {endpoint!r} but mngmt/decisions.md does not record it"
    )


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
