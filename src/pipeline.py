"""Orchestrator: load config, seed RNG, run L0 -> L5 in order.

Purpose
    Single entry point for the whole repurposing pipeline. Owns config loading, RNG
    seeding, and layer sequencing so that no layer has to know about any other layer's
    location on disk.

Inputs
    ``config/pipeline.yaml`` (path overridable via ``--config``). Keys: ``seed``,
    ``data_dir``, ``results_dir``, ``causal_gene``, ``therapeutic_endpoint``, ``channels``.

Outputs
    Artifacts under ``config["results_dir"]`` (gitignored). Nothing is written anywhere
    else, and nothing patient-derived leaves that directory.

Guardrail
    ``causal_gene`` and ``therapeutic_endpoint`` are read from config and are ``None``
    until deliberately set by a human. This module must never supply a default for
    either -- guessing the causal gene would silently fabricate the scientific premise
    of the entire run.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "pipeline.yaml"

#: Layer order is fixed: each layer consumes the previous layer's output.
LAYER_ORDER = (
    "l0_genomics",
    "l1_target",
    "l2_channels",
    "l3_integrate",
    "l4_validate",
    "l5_report",
)


def load_config(path: Path = DEFAULT_CONFIG_PATH) -> dict:
    """Load the pipeline configuration from a YAML file.

    Args:
        path: Path to the YAML config. Defaults to ``config/pipeline.yaml``.

    Returns:
        The parsed configuration mapping.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: read `path` with `yaml.safe_load`; validate that `seed`, `data_dir`,
    # `results_dir`, and `channels` are present; resolve `data_dir`/`results_dir` to
    # absolute paths relative to the repo root; and warn (do NOT default) when
    # `causal_gene` or `therapeutic_endpoint` is still None.
    raise NotImplementedError("load_config is a scaffold stub")


def seed_everything(seed: int) -> None:
    """Seed every RNG the pipeline touches, for reproducibility.

    Args:
        seed: The seed value from ``config["seed"]``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: seed `random`, `numpy.random`, and PYTHONHASHSEED; seed any framework RNG
    # (torch/DGL) used by the L2 knowledge-graph channel. A judge re-running this on
    # freshly downloaded data must reproduce the ranking exactly.
    raise NotImplementedError("seed_everything is a scaffold stub")


def run(config: dict) -> None:
    """Run every layer L0 -> L5 in order.

    Args:
        config: The parsed pipeline configuration.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: import each layer's `run` lazily (so a missing optional dependency in one
    # layer does not break the others), call them in LAYER_ORDER, pass the accumulated
    # artifact manifest forward, and checkpoint after each layer into `results_dir` so a
    # failed L4 does not force a re-run of the expensive L2 channels.
    raise NotImplementedError("pipeline.run is a scaffold stub")


def main(argv: list | None = None) -> int:
    """CLI entry point: ``python -m src.pipeline``.

    Args:
        argv: Argument vector; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    parser = argparse.ArgumentParser(description="MVA Track 2 drug repurposing pipeline")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    LOGGER.info("Loading config from %s", args.config)

    config = load_config(args.config)
    seed_everything(config["seed"])
    run(config)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
