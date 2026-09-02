"""Data custody purge and destruction attestation.

Purpose
    Destroy every location that held patient-derived bytes, verify each one is empty
    afterwards, and emit an auditable attestation of that fact. This is the machinery
    behind the deletion commitment in ``COMPLIANCE.md`` -- delete all data within 30 days
    of Hackathon close and email confirmation to the organizers.

    The register in ``docs/data_custody.md`` is the input, not this module's guesswork. A
    destruction record reconstructed from memory at deletion time omits things; custody is
    instrumented when each location is *created*, and this module only attests to what was
    registered.

Inputs
    - ``docs/data_custody.md`` -- the custody register (locations, what they hold, purge
      method, status).
    - The resolved environment: ``MVA_DATA_ROOT``, ``HF_HOME``, ``config["results_dir"]``.
    - The dataset's Hub revision SHA, recorded at download time.

Outputs
    - ``docs/purge_attestation_<UTC date>.json`` -- schema below, committed.
    - A rendered Markdown companion for the confirmation email body.
    - An updated register: every row moved to ``purged`` with its verification timestamp.

Guardrail
    **The attestation must contain no patient-derived content, by construction.** No
    filenames (dataset filenames may embed a subject identifier), no per-file digests (a
    file hash is a fingerprint that lets a holder confirm they hold the same file), no
    variant coordinates, depths, or phenotype text. The destroyed data is identified by
    the dataset's **Hub revision SHA** only -- that names the exact snapshot granted, is
    verifiable against the organizers' own records, and reveals nothing about the child.
    Aggregate file counts and total bytes per location are permitted: they evidence
    completeness without describing content.

    **Refuse to attest beyond the register.** A location absent from
    ``docs/data_custody.md`` makes the attestation false. This module fails loudly on an
    unregistered path rather than silently purging it, because a purge that succeeds
    quietly is indistinguishable from a register that was never complete.

    **Destructive by explicit consent only.** Dry-run is the default; deletion requires an
    explicit confirmation flag. Never delete a path that is not a registered location, and
    never delete the repository working tree.

    See ``docs/data_custody.md`` for the register, the two cache traps this defends
    against, and the honest limits of any destruction record.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

#: Bumped when the attestation JSON shape changes. Recorded in every attestation.
ATTESTATION_SCHEMA_VERSION = "1.0"

REPO_ROOT = Path(__file__).resolve().parent.parent

#: The custody register. Authoritative: nothing is purged or attested that is absent here.
CUSTODY_REGISTER = REPO_ROOT / "docs" / "data_custody.md"

#: The gated dataset, identified for attestation by repo id + revision SHA, never by file.
DATASET_REPO_ID = "SageBio/mva-hackathon-2026-data"
DATASET_REPO_TYPE = "dataset"

#: Register rows move through these states. ``purged`` requires post-deletion verification.
CUSTODY_STATES = ("planned", "active", "purged")

#: Checks that must all pass before an attestation may be emitted.
REQUIRED_CHECKS = (
    "custody_register_complete",     # every register row reached a terminal state
    "repo_scan_genomic_extensions",  # tests/test_smoke.py::test_no_patient_data_files_present
    "hf_cache_empty",                # `hf cache ls` reports zero repos, zero bytes
    "incomplete_blobs_removed",      # `hf cache prune` -- .incomplete partials survive `rm`
)

#: Field names that must never appear anywhere in an emitted attestation. Enforced by an
#: outbound scan, because the guardrail is only real if something checks it.
FORBIDDEN_ATTESTATION_FIELDS = (
    "file_name",
    "file_names",
    "filenames",
    "sha256",
    "md5",
    "checksum",
    "digest",
    "sample_id",
    "subject_id",
    "variant",
    "coordinates",
)


@dataclass(frozen=True)
class LocationState:
    """Aggregate state of one custody location at a point in time.

    Deliberately holds no per-file detail -- counts and bytes evidence completeness
    without describing content.
    """

    exists: bool
    file_count: int
    total_bytes: int


@dataclass(frozen=True)
class CustodyLocation:
    """One row of the custody register."""

    id: str           # e.g. "C1"
    path: Path
    holds: str        # prose description; never an enumeration of filenames
    method: str       # e.g. "rm -rf; fstrim"
    status: str       # one of CUSTODY_STATES


@dataclass(frozen=True)
class LocationVerdict:
    """Before/after evidence for one location, as it appears in the attestation."""

    location: CustodyLocation
    before: LocationState
    after: LocationState
    verified_utc: str

    @property
    def is_clean(self) -> bool:
        """True when the location is absent, or present and empty."""
        return not self.after.exists or (self.after.file_count == 0 and self.after.total_bytes == 0)


def read_register(path: Path = CUSTODY_REGISTER) -> list[CustodyLocation]:
    """Parse the custody register into structured locations.

    Args:
        path: Path to ``docs/data_custody.md``.

    Returns:
        One :class:`CustodyLocation` per register row.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: parse the Markdown table under "## Custody register", expanding
    # `$MVA_DATA_ROOT`, `$HF_HOME`, and `<repo>` against the live environment. Fail on an
    # unknown status value, on a duplicate id, and on an unexpanded placeholder -- a
    # register row that does not resolve to a real path is a row that will not be purged.
    raise NotImplementedError("purge.read_register is a scaffold stub")


def measure(location: CustodyLocation) -> LocationState:
    """Measure one location without recording anything about its contents.

    Args:
        location: The registered location to stat.

    Returns:
        Existence, file count, and total bytes.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: walk the tree accumulating a count and a byte total. Do not collect names,
    # extensions, or digests -- not even transiently into a log. Count symlinks and their
    # targets separately: the Hub cache stores real bytes in `blobs/` and only symlinks in
    # `snapshots/`, so a naive walk that follows links double-counts, and one that ignores
    # blobs/ reports an empty location that still holds 85 GB.
    raise NotImplementedError("purge.measure is a scaffold stub")


def purge(config: dict, *, confirm: bool = False) -> list[LocationVerdict]:
    """Delete every registered location and verify each is empty afterwards.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir``.
        confirm: Must be ``True`` to delete. Default performs a dry run.

    Returns:
        One verdict per registered location.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: for each registered location -- measure(), then delete by its declared method,
    # then measure() again and stamp the verification time. Order matters: purge the Hub
    # cache with `hf cache rm` followed by `hf cache prune` BEFORE removing the cache root,
    # so `.incomplete` partial blobs are accounted for rather than merely unlinked. Run
    # `fstrim -av` last, once every location is gone, so freed blocks are discarded.
    #
    # Refuse, loudly, to touch: any path absent from the register; any path inside the
    # repository other than `results/` and notebook outputs; and any path that resolves
    # outside `MVA_DATA_ROOT` unless the register explicitly names it.
    #
    # With confirm=False, do everything except the deletion and report what would change.
    raise NotImplementedError("purge.purge is a scaffold stub")


def write_attestation(verdicts: list[LocationVerdict], *, revision_sha: str) -> Path:
    """Emit the destruction attestation as JSON plus a Markdown companion.

    Args:
        verdicts: Per-location before/after evidence from :func:`purge`.
        revision_sha: The dataset's full 40-character Hub commit hash -- the only
            identifier of the destroyed data that appears in the attestation.

    Returns:
        Path to the written JSON attestation.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: assemble the schema documented in docs/data_custody.md, then -- before writing
    # -- scan the serialized payload for FORBIDDEN_ATTESTATION_FIELDS and for anything
    # resembling a genomic filename or a hex digest, and refuse to write on a hit. The
    # guardrail against leaking patient data into the destruction record is only real if
    # something enforces it at the boundary.
    #
    # Record every REQUIRED_CHECKS result, the repo's git SHA, hostname, platform, and WSL
    # distro. Refuse to emit if any verdict is not `is_clean` or any required check failed:
    # a partial purge must not be reported as a completed one.
    raise NotImplementedError("purge.write_attestation is a scaffold stub")


def main(argv: list | None = None) -> int:
    """CLI entry point: ``python -m src.purge [--confirm]``.

    Args:
        argv: Argument vector; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    parser = argparse.ArgumentParser(description="Purge patient-derived data and attest to it")
    parser.add_argument(
        "--confirm",
        action="store_true",
        help="Actually delete. Without this flag the run is a dry run and nothing is removed.",
    )
    parser.add_argument("--register", type=Path, default=CUSTODY_REGISTER)
    args = parser.parse_args(argv)

    # TODO: load config, read the register, purge (honouring args.confirm), then attest.
    # A dry run must still print the full before-state and the planned method per location,
    # so the operator can audit the plan before authorising destruction.
    raise NotImplementedError("purge.main is a scaffold stub")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
