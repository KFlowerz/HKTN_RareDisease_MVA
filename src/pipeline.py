"""Orchestrator: load config, seed RNG, run L0 -> L5 in order.

Purpose
    Single entry point for the whole repurposing pipeline. Owns config loading, RNG
    seeding, and layer sequencing so that no layer has to know about any other layer's
    location on disk.

Inputs
    ``config/pipeline.yaml`` (path overridable via ``--config``). Keys: ``seed``,
    ``data_dir``, ``results_dir``, ``causal_gene``, ``therapeutic_endpoint``, ``model``,
    ``channels``.

    ``MVA_DATA_ROOT`` (environment) overrides ``data_dir``, resolving to
    ``$MVA_DATA_ROOT/raw`` -- location ``C1`` in ``docs/data_custody.md``. The committed
    config cannot name that path: it is machine-specific and lives outside the repo, and
    hardcoding an absolute home directory would break the reproducibility claim for a
    judge running this on their own re-downloaded copy.

Outputs
    Artifacts under ``config["results_dir"]`` (gitignored), plus ``_manifest.json``
    recording each layer's status, wall time, and the artifacts it produced. Nothing is
    written anywhere else, and nothing patient-derived leaves that directory.

Guardrail
    ``causal_gene`` and ``therapeutic_endpoint`` are read from config and are ``None``
    until deliberately set by a human. This module must never supply a default for
    either -- guessing the causal gene would silently fabricate the scientific premise
    of the entire run. It warns loudly and continues; it does not fill them in.

    ``results_dir`` may not sit inside ``data_dir``. Layers read from one and write to
    the other, and nesting them would let a derived artifact be mistaken for source data
    -- or be missed by a purge scoped to one tree (``docs/data_custody.md``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import random
import time
from datetime import datetime, timezone
from pathlib import Path

LOGGER = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "pipeline.yaml"

#: Layer order is fixed: each layer consumes the previous layer's output.
LAYER_ORDER = (
    "l0_genomics",
    "l1_target",
    "l2_channels",
    "l3_integrate",
    "l4_validate",
    "l5_report",
)

#: Keys that must be present in the config file. The two open decisions are
#: deliberately absent -- they may be None, but they must not be missing.
REQUIRED_KEYS = ("seed", "data_dir", "results_dir", "channels")

#: Written after every layer so a failed L4 does not force a re-run of the expensive
#: L2 channels or the paid L3 reasoning step.
MANIFEST_NAME = "_manifest.json"

#: Layers whose output directory is not named after the layer. Without an entry here the
#: manifest looks in the wrong place and silently records a completed layer as producing
#: nothing, which reads as a failure that did not happen. L0 and L1 write to directories
#: named after themselves; L2-L5 write to short names.
#: ``tests/test_pipeline.py`` asserts every layer's directory matches what the
#: orchestrator looks for, because the failure is invisible -- the layer still completes.
LAYER_OUTPUT_DIRS = {"l2_channels": "l2", "l3_integrate": "l3", "l4_validate": "l4",
                     "l5_report": "l5"}


def load_config(path: Path = DEFAULT_CONFIG_PATH) -> dict:
    """Load and validate the pipeline configuration.

    Args:
        path: Path to the YAML config. Defaults to ``config/pipeline.yaml``.

    Returns:
        The parsed configuration, with ``data_dir`` and ``results_dir`` resolved to
        absolute paths and ``config_path`` recorded for provenance.

    Raises:
        FileNotFoundError: If the config file or the resolved ``data_dir`` is absent.
        ValueError: If a required key is missing, or ``results_dir`` nests inside
            ``data_dir``.
    """
    import yaml  # imported here so a missing optional dep cannot break `--help`

    if not path.is_file():
        raise FileNotFoundError(f"config not found: {path}")

    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    missing = [k for k in REQUIRED_KEYS if k not in config]
    if missing:
        raise ValueError(f"config {path} is missing required key(s): {', '.join(missing)}")

    # MVA_DATA_ROOT wins over the committed value. The dataset lives in `raw/` beneath
    # the custody root; see docs/data_custody.md location C1.
    data_root = os.environ.get("MVA_DATA_ROOT")
    if data_root:
        data_dir = Path(data_root).expanduser() / "raw"
        LOGGER.info("data_dir from MVA_DATA_ROOT: %s", data_dir)
    else:
        data_dir = (REPO_ROOT / str(config["data_dir"])).resolve()
        LOGGER.warning(
            "MVA_DATA_ROOT is not set; falling back to config data_dir %s. "
            "Load the environment first: set -a; . ~/.config/mva/env; set +a",
            data_dir,
        )

    results_dir = (REPO_ROOT / str(config["results_dir"])).resolve()

    if data_dir == results_dir or results_dir.is_relative_to(data_dir):
        raise ValueError(
            f"results_dir ({results_dir}) must not sit inside data_dir ({data_dir}): "
            "derived artifacts would be indistinguishable from source data"
        )
    if not data_dir.is_dir():
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")

    results_dir.mkdir(parents=True, exist_ok=True)

    config["data_dir"] = data_dir
    config["results_dir"] = results_dir
    config["config_path"] = path

    # Warn, never default. These are findings and decisions, not settings.
    if config.get("causal_gene") is None:
        LOGGER.warning(
            "causal_gene is null -- it is a finding from L0, never a setting. The dataset "
            "names no gene and there is no Track-1 answer to reconcile against, so L0 must "
            "call it independently and ship the evidence (gate G1)"
        )
    if config.get("therapeutic_endpoint") is None:
        LOGGER.warning(
            "therapeutic_endpoint is null -- L1 will emit BOTH target sets and L3 will "
            "carry both forward, roughly doubling L4 annotation work (gate G2). It was "
            "decided on 2026-09-08; see mngmt/decisions.md D4 before setting it back to null"
        )

    return config


def seed_everything(seed: int) -> None:
    """Seed every RNG the pipeline touches.

    Args:
        seed: The seed value from ``config["seed"]``.

    Note:
        ``PYTHONHASHSEED`` cannot be changed after the interpreter starts. Setting it
        here affects **child processes only**; if it was unset at launch, str/bytes hash
        randomization is already active in this process and any ordering derived from
        set or dict iteration over strings may vary between runs. That is warned about
        rather than papered over -- a reproducibility claim that quietly depends on an
        environment variable nobody set is worse than a documented caveat.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)

    try:
        import numpy as np

        np.random.seed(seed)  # legacy global, for libraries that still use it
    except ImportError:  # pragma: no cover - numpy is a hard dep in practice
        LOGGER.warning("numpy not importable; skipping numpy seeding")

    for mod_name, seeder in (
        ("torch", lambda m: m.manual_seed(seed)),
        ("dgl", lambda m: m.seed(seed)),
    ):
        try:
            seeder(__import__(mod_name))
            LOGGER.debug("seeded %s", mod_name)
        except Exception:  # noqa: BLE001 - absence or API drift must not break the run
            pass

    LOGGER.info("seeded all RNGs with %d", seed)


#: Libraries whose version changes the numbers, not just the wrapping. Recorded per run so
#: two artifacts that disagree can be traced to the thing that differed. snpEff and
#: bcftools are deliberately absent -- they are external executables, and L0 records its
#: own tool versions where it invokes them.
NUMERIC_DEPENDENCIES = ("numpy", "scipy", "pandas", "networkx", "sklearn", "h5py",
                        "pyarrow", "matplotlib")


def _environment() -> dict:
    """Interpreter and library versions, plus whether hash randomization was fixed.

    An unpinned environment is the failure this records: results computed with a different
    numpy are close, different, and silently so. Nothing in a TSV of floats says which
    library produced it, so the manifest has to.
    """
    import importlib
    import platform
    import sys

    versions = {}
    for name in NUMERIC_DEPENDENCIES:
        try:
            versions[name] = getattr(importlib.import_module(name), "__version__", "unknown")
        except ImportError:
            versions[name] = "not installed"
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
        # Set at launch or not at all: Python fixes the hash seed before any code runs.
        "pythonhashseed": os.environ.get("PYTHONHASHSEED", "unset"),
    }


#: Config keys that are locations rather than science. Excluded from the settings digest:
#: they are absolute and differ on every machine, so including them would make the digest
#: incomparable between two people running the identical analysis -- which is the one
#: comparison it exists to support.
MACHINE_SPECIFIC_KEYS = frozenset({
    "config_path", "data_dir", "results_dir", "reference_dir", "enrichment_dir",
})


def _config_digest(config: dict) -> dict:
    """SHA-256 of the config file's bytes, and of the settings that drove the run.

    Two digests because they answer different questions. ``file_sha256`` pins the exact
    bytes on disk, so an edited config is visible even if nothing meaningful changed.
    ``settings_sha256`` covers the resolved values with the machine-specific paths removed
    (:data:`MACHINE_SPECIFIC_KEYS`), so two judges running the same analysis on different
    machines get the *same* digest and a genuine difference stands out.
    """
    path = Path(config.get("config_path") or "")
    file_digest = None
    if path.is_file():
        file_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    settings = json.dumps({k: str(v) for k, v in sorted(config.items())
                           if k not in MACHINE_SPECIFIC_KEYS}, sort_keys=True)
    return {"file_sha256": file_digest,
            "settings_sha256": hashlib.sha256(settings.encode("utf-8")).hexdigest()}


def _read_manifest(results_dir: Path) -> dict:
    """Return the run manifest, or a fresh one if none exists."""
    path = results_dir / MANIFEST_NAME
    if path.is_file():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            LOGGER.warning("manifest at %s is unreadable; starting a new one", path)
    return {"layers": {}}


def _write_manifest(results_dir: Path, manifest: dict) -> None:
    """Persist the manifest atomically, so an interrupted write cannot corrupt it."""
    path = results_dir / MANIFEST_NAME
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _artifacts(results_dir: Path, name: str) -> list:
    """What this layer produced in this run.

    Usually every file under the layer's directory. L2 is the exception: its channels have
    their own directories, and a channel that was disabled this run may still have a table
    on disk from an earlier one. Listing that as an artifact of this run would contradict
    ``channels.json``, which is what L3 reads -- so for L2 the status record is the
    authority, and only the completed channels' files are listed.
    """
    layer_dir = results_dir / LAYER_OUTPUT_DIRS.get(name, name)
    if not layer_dir.is_dir():
        return []
    index = layer_dir / "channels.json"
    if name == "l2_channels" and index.is_file():
        try:
            status = json.loads(index.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            LOGGER.warning("%s is unreadable; listing no channel artifact", index)
            return [str(index.relative_to(results_dir))]
        listed = [a for record in status.get("channels", {}).values()
                  if record.get("status") == "complete"
                  for a in record.get("artifacts", ())]
        return sorted({str(index.relative_to(results_dir)), *listed})
    return sorted(str(p.relative_to(results_dir)) for p in layer_dir.rglob("*") if p.is_file())


def run(config: dict, *, resume: bool = False, only: str | None = None) -> None:
    """Run every layer L0 -> L5 in order.

    Args:
        config: The parsed pipeline configuration.
        resume: Skip layers already recorded ``complete`` in the manifest.
        only: Run just this one layer (a name from :data:`LAYER_ORDER`).

    Raises:
        ValueError: If ``only`` names an unknown layer.
        Exception: Whatever a layer raises. The manifest records the failure before it
            propagates, so a re-run can resume rather than restart.
    """
    import importlib

    if only is not None and only not in LAYER_ORDER:
        raise ValueError(f"unknown layer {only!r}; expected one of {', '.join(LAYER_ORDER)}")

    results_dir: Path = config["results_dir"]
    manifest = _read_manifest(results_dir)
    manifest["run_started_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["seed"] = config["seed"]
    manifest["config_path"] = str(config.get("config_path", ""))
    manifest["config"] = _config_digest(config)
    manifest["environment"] = _environment()
    manifest["causal_gene"] = config.get("causal_gene")
    manifest["therapeutic_endpoint"] = config.get("therapeutic_endpoint")

    layers = (only,) if only else LAYER_ORDER

    for name in layers:
        record = manifest["layers"].get(name, {})
        if resume and record.get("status") == "complete":
            LOGGER.info("skipping %s (already complete; --resume)", name)
            continue

        LOGGER.info("--- %s ---", name)
        started = time.monotonic()
        try:
            # Imported lazily so a missing optional dependency in one layer does not
            # prevent the others from running.
            module = importlib.import_module(f"src.{name}.run")
            module.run(config)
        except NotImplementedError as exc:
            manifest["layers"][name] = {
                "status": "not_implemented",
                "detail": str(exc),
                "seconds": round(time.monotonic() - started, 3),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
            _write_manifest(results_dir, manifest)
            # A layer may raise this *after* doing real work -- L0 writes its burden
            # artifacts, then stops because its causal-gene half needs a local annotator.
            # Report the layer's own message rather than calling everything a stub.
            LOGGER.warning("%s stopped: %s", name, exc)
            return
        except Exception as exc:
            manifest["layers"][name] = {
                "status": "failed",
                "detail": f"{type(exc).__name__}: {exc}",
                "seconds": round(time.monotonic() - started, 3),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
            _write_manifest(results_dir, manifest)
            raise

        artifacts = _artifacts(results_dir, name)
        manifest["layers"][name] = {
            "status": "complete",
            "seconds": round(time.monotonic() - started, 3),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "artifacts": artifacts,
        }
        _write_manifest(results_dir, manifest)
        LOGGER.info("%s complete in %.1fs (%d artifact(s))", name, manifest["layers"][name]["seconds"], len(artifacts))


def main(argv: list | None = None) -> int:
    """CLI entry point: ``python -m src.pipeline``.

    Args:
        argv: Argument vector; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(description="MVA Track 2 drug repurposing pipeline")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--resume", action="store_true", help="skip layers already complete")
    parser.add_argument("--only", metavar="LAYER", help=f"run one layer: {', '.join(LAYER_ORDER)}")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    LOGGER.info("Loading config from %s", args.config)

    if "PYTHONHASHSEED" not in os.environ:
        LOGGER.warning(
            "PYTHONHASHSEED was unset at launch; hash randomization is active in this "
            "process and cannot be disabled retroactively. For a bit-reproducible run, "
            "launch with PYTHONHASHSEED=<seed>."
        )

    config = load_config(args.config)
    seed_everything(config["seed"])
    run(config, resume=args.resume, only=args.only)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
