"""L2 -- Channel fan-out.

Purpose
    Dispatch to each candidate-generation channel enabled in
    ``config["channels"]`` and record, per channel, whether it produced a ranked table --
    so L3 aggregates only over channels that actually ran in this invocation.

Inputs
    L1 output (disease module, upstream/downstream target sets) and
    ``config["channels"]`` -- a mapping of ``{kg, proximity, signature, phenotype,
    prior} -> bool``.

Outputs
    Under ``config["results_dir"]/l2/``:
      - ``<channel>/candidates.tsv`` (and whatever else the channel writes) -- one ranked
        table per channel that completed, each row carrying the channel's own score, its
        rank and per-candidate provenance
      - ``channels.json`` -- every registered channel's status in this run
        (``complete`` | ``failed`` | ``not_implemented`` | ``disabled``), wall time,
        failure detail, artifacts and candidate count. **L3 reads this, not the
        directory listing**: a table on disk from an earlier run is not evidence from
        this one.

Guardrail
    Channels run **independently**. No channel may read another's output or scores --
    L3's rank aggregation assumes independent evidence, and leaking scores between
    channels would manufacture false consensus. Each channel therefore gets its own deep
    copy of the config, and the global RNGs are reseeded from ``config["seed"]`` before
    each one, so a channel's result cannot depend on which channels ran before it.

    A channel that fails is recorded as *failed*, never silently treated as "returned no
    candidates": an empty list is a scientific claim, an exception is not. For the same
    reason a channel that returns without writing a non-empty ``candidates.tsv`` is
    recorded as failed, and its output directory is emptied before it runs, so a partial
    or stale table cannot survive a failure and be read as current.

    The layer does not report success over a gap. If any enabled channel failed, the
    layer raises after writing ``channels.json``; if any enabled channel is not yet
    implemented, it raises :class:`NotImplementedError`, which the orchestrator records as
    ``not_implemented``. The tables of channels that did complete stay on disk either way.
"""

from __future__ import annotations

import copy
import json
import logging
import random
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

from . import (
    channel_a_kg,
    channel_b_proximity,
    channel_c_signature,
    channel_d_phenotype,
    channel_e_prior,
)

LOGGER = logging.getLogger(__name__)

#: Maps the config key under ``channels:`` to the module implementing that channel.
CHANNEL_REGISTRY = {
    "kg": channel_a_kg,
    "proximity": channel_b_proximity,
    "signature": channel_c_signature,
    "phenotype": channel_d_phenotype,
    "prior": channel_e_prior,
}

#: The ranked table every channel must write into its own directory.
CANDIDATES_FILE = "candidates.tsv"
#: The per-run status record, beside the channel directories.
STATUS_FILE = "channels.json"


def output_root(config: dict) -> Path:
    """``results_dir/l2`` -- where every channel writes its own subdirectory."""
    return Path(config["results_dir"]) / "l2"


def channel_dir(config: dict, key: str) -> Path:
    """A channel's output directory, named after its module (``channel_b_proximity``)."""
    return output_root(config) / CHANNEL_REGISTRY[key].__name__.rsplit(".", 1)[-1]


def enabled_channels(config: dict) -> list:
    """The enabled channel keys, in registry order.

    Raises:
        ValueError: On an unknown key, a registered channel missing from config, or a
            non-boolean toggle. A typo must not silently disable a line of evidence, and
            ``"false"`` in YAML quotes is a string, which is truthy.
    """
    toggles = config.get("channels")
    if not isinstance(toggles, dict):
        raise ValueError("config['channels'] must map channel names to true/false")
    unknown = sorted(set(toggles) - set(CHANNEL_REGISTRY))
    missing = sorted(set(CHANNEL_REGISTRY) - set(toggles))
    if unknown or missing:
        raise ValueError(f"config['channels'] does not match the registry: unknown {unknown}, "
                         f"missing {missing}; expected exactly {sorted(CHANNEL_REGISTRY)}")
    not_bool = sorted(k for k, v in toggles.items() if not isinstance(v, bool))
    if not_bool:
        raise ValueError(f"config['channels'] toggles must be true or false, not strings or "
                         f"numbers: {not_bool}")
    return [key for key in CHANNEL_REGISTRY if toggles[key]]


def _reseed(seed: int) -> None:
    """Reset the global RNGs, so no channel inherits another's consumption of them."""
    random.seed(seed)
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - numpy is a hard dependency in practice
        return
    np.random.seed(seed)


def _count_rows(path: Path) -> int:
    """Data rows in a TSV with a header."""
    with open(path, encoding="utf-8") as handle:
        return max(0, sum(1 for line in handle if line.strip()) - 1)


def run_channel(config: dict, key: str) -> dict:
    """Run one channel in isolation and return its status record.

    A channel's own exception is caught and recorded, never propagated -- the other
    channels still run. The only raise is the path guard below, which would mean the
    output layout itself is wrong.
    """
    module = CHANNEL_REGISTRY[key]
    out_dir = channel_dir(config, key)
    root = output_root(config)
    # Emptied first: a table left by an earlier run must not survive this one's failure.
    if out_dir.exists():
        if root.resolve() not in out_dir.resolve().parents:
            raise ValueError(f"refusing to clear {out_dir}: not inside {root}")
        shutil.rmtree(out_dir)

    record = {"module": module.__name__, "output_dir": str(out_dir.relative_to(root.parent))}
    _reseed(config["seed"])
    started = time.monotonic()
    try:
        module.generate(copy.deepcopy(config))
    except NotImplementedError as exc:
        record.update(status="not_implemented", detail=str(exc))
    except Exception as exc:  # noqa: BLE001 -- recorded, then reported by run()
        # Only the type reaches the log. The message goes to channels.json under the
        # gitignored results_dir: a channel reading patient-derived input (channel D
        # parses phenotype terms) could raise with that content in its message.
        LOGGER.error("channel %s failed: %s (detail in %s)", key, type(exc).__name__,
                     STATUS_FILE)
        record.update(status="failed", detail=f"{type(exc).__name__}: {exc}")
    else:
        table = out_dir / CANDIDATES_FILE
        n_rows = _count_rows(table) if table.is_file() else 0
        if n_rows == 0:
            record.update(status="failed", detail=(
                f"returned without writing a non-empty {CANDIDATES_FILE}. A channel that "
                "means 'no candidates' must say so explicitly; silence is not a result."))
        else:
            record.update(status="complete", n_candidates=n_rows)
    record["seconds"] = round(time.monotonic() - started, 3)
    record["artifacts"] = (sorted(str(p.relative_to(root.parent)) for p in out_dir.rglob("*")
                                  if p.is_file()) if out_dir.is_dir() else [])
    return record


def run(config: dict) -> None:
    """Execute every enabled L2 channel and write ``channels.json``.

    Args:
        config: Parsed pipeline configuration. Uses ``channels``, ``results_dir``,
            ``seed``, and whatever each channel reads.

    Raises:
        ValueError: If ``config["channels"]`` does not match the registry.
        RuntimeError: If any enabled channel failed.
        NotImplementedError: If none failed but an enabled channel is not implemented.
    """
    # Removed before anything can fail -- including config validation -- so a status
    # record from an earlier run never outlives a run that did not finish.
    if "results_dir" in config:
        (output_root(config) / STATUS_FILE).unlink(missing_ok=True)
    enabled = enabled_channels(config)
    root = output_root(config)
    root.mkdir(parents=True, exist_ok=True)

    channels = {}
    for key in CHANNEL_REGISTRY:
        if key not in enabled:
            channels[key] = {"module": CHANNEL_REGISTRY[key].__name__, "status": "disabled"}
            continue
        LOGGER.info("channel %s: starting", key)
        channels[key] = run_channel(config, key)
        LOGGER.info("channel %s: %s in %.1fs", key, channels[key]["status"],
                    channels[key]["seconds"])

    by_status = {s: sorted(k for k, r in channels.items() if r["status"] == s)
                 for s in ("complete", "failed", "not_implemented", "disabled")}
    (root / STATUS_FILE).write_text(json.dumps({
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": config["seed"],
        "causal_gene": config.get("causal_gene"),
        "channels": channels,
        "summary": {status: keys for status, keys in by_status.items()},
        "n_complete": len(by_status["complete"]),
        "note": ("Only channels with status 'complete' produced evidence in this run. "
                 "Tables are never read by directory listing; a channel absent from "
                 "'complete' contributes nothing, and is not a finding of 'no candidates'."),
    }, indent=2, sort_keys=True), encoding="utf-8")

    if by_status["failed"]:
        raise RuntimeError(
            f"L2 channel(s) failed: {', '.join(by_status['failed'])} -- see "
            f"{root / STATUS_FILE}. Completed: {', '.join(by_status['complete']) or 'none'}.")
    if by_status["not_implemented"]:
        raise NotImplementedError(
            f"L2 channel(s) enabled but not implemented: "
            f"{', '.join(by_status['not_implemented'])}. Completed: "
            f"{', '.join(by_status['complete']) or 'none'}. Disable them in "
            "config['channels'] to let the layer complete without them.")
