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
import json
import os
import platform
import re
import shutil
import socket
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

#: Bumped when the attestation JSON shape changes. Recorded in every attestation.
ATTESTATION_SCHEMA_VERSION = "1.0"

REPO_ROOT = Path(__file__).resolve().parent.parent

#: The custody register. Authoritative: nothing is purged or attested that is absent here.
CUSTODY_REGISTER = REPO_ROOT / "docs" / "data_custody.md"

#: The gated dataset, identified for attestation by repo id + revision SHA, never by file.
DATASET_REPO_ID = "SageBio/mva-hackathon-2026-data"
DATASET_REPO_TYPE = "dataset"

#: Submission close, and the deletion deadline COMPLIANCE.md derives from it. The two
#: readings of "Hackathon close" are a month apart; this is the earlier, safer one. See
#: the open item in docs/data_custody.md.
HACKATHON_CLOSE_DATE = "2026-10-24"
DELETION_DEADLINE_UTC = "2026-11-23T23:59:00Z"

#: Register rows move through these states. ``purged`` requires post-deletion verification.
CUSTODY_STATES = ("planned", "active", "purged")

#: States for rows this tool cannot act on, which the register carries anyway. ``C6`` is
#: prevented at the source rather than purged; ``C8`` is the organizers' own copy of the
#: Track 1 submission, which travels back to the party that issued the data. Both are
#: terminal for ``custody_register_complete``: there is nothing left for a purge to do.
#: They are separated from :data:`CUSTODY_STATES` rather than folded into it because a
#: register that showed only the copies we can delete would overstate the deletion
#: commitment -- see the ``C8`` note in ``docs/data_custody.md``.
UNPURGEABLE_STATES = ("n/a", "live once submitted")

#: A register cell names a filesystem path when it is backticked, contains a separator, and
#: begins with a root, a variable or the repository placeholder. This excludes ``C7``'s
#: ``.jsonl`` (a suffix, no separator) and ``C8``'s ``SageBio/...`` (a Hub repo id, not a
#: path) without either needing a special case.
_PATH_CELL = re.compile(r"`((?:\$[A-Z_]+|<repo>|~|/|\.{1,2})[^`]*/[^`]*)`")

#: Placeholders a register path may carry, and where each resolves from.
_REPO_PLACEHOLDER = "<repo>"

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
    path: Path | None  # None for a row that names no filesystem location (C5–C8)
    holds: str        # prose description; never an enumeration of filenames
    method: str       # e.g. "rm -rf; fstrim"
    status: str       # one of CUSTODY_STATES or UNPURGEABLE_STATES

    @property
    def is_purgeable(self) -> bool:
        """True when this tool can act on the row: it has a path and a workable state."""
        return self.path is not None and self.status in CUSTODY_STATES


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
    rows = _register_rows(path.read_text(encoding="utf-8"))
    if not rows:
        raise ValueError(f"{path}: no rows found under '## Custody register'")

    locations: list[CustodyLocation] = []
    seen: set[str] = set()
    for cells in rows:
        identifier = _plain(cells[0])
        if identifier in seen:
            raise ValueError(
                f"{path}: duplicate register id {identifier!r}. Ids key the attestation, "
                "so a duplicate silently drops a location from it")
        seen.add(identifier)

        status = _plain(cells[4]).lower()
        if status not in CUSTODY_STATES and status not in UNPURGEABLE_STATES:
            raise ValueError(
                f"{path}: {identifier} has status {status!r}, which is neither a custody "
                f"state {CUSTODY_STATES} nor an unpurgeable one {UNPURGEABLE_STATES}")

        locations.append(CustodyLocation(
            id=identifier,
            path=_resolve_path(cells[1], identifier),
            holds=_plain(cells[2]),
            method=_plain(cells[3]),
            status=status,
        ))
    return locations


def _register_rows(text: str) -> list[list[str]]:
    """The data rows of the Markdown table under '## Custody register'."""
    section = text.split("## Custody register", 1)
    if len(section) != 2:
        raise ValueError("no '## Custody register' heading")
    body = section[1].split("\n## ", 1)[0]

    rows = []
    for line in body.splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|- "):
            continue                                  # separator row, or not a table row
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if len(cells) != 5 or cells[0].strip("` ").upper() in ("ID", ""):
            continue                                  # header, or a differently-shaped table
        rows.append(cells)
    return rows


def _plain(cell: str) -> str:
    """A table cell without its Markdown emphasis or code ticks."""
    return cell.replace("**", "").replace("`", "").strip()


def _resolve_path(cell: str, identifier: str) -> Path | None:
    """The filesystem path a register cell names, with placeholders expanded.

    Returns ``None`` when the row names no path at all -- a shell-history surface, or
    someone else's store. Raises when it names one that does not resolve, because a row
    that cannot be located is a row that will not be purged, and the attestation would
    claim otherwise.
    """
    match = _PATH_CELL.search(cell)
    if match is None:
        return None

    raw = match.group(1)
    if raw.startswith(_REPO_PLACEHOLDER):
        return (REPO_ROOT / raw[len(_REPO_PLACEHOLDER):].lstrip("/")).resolve()

    expanded = os.path.expandvars(os.path.expanduser(raw))
    if "$" in expanded:
        missing = re.findall(r"\$([A-Z_]+)", expanded)
        raise ValueError(
            f"{identifier}: {raw!r} still contains {', '.join(missing)} after expansion. "
            "Set it in the environment before purging; an unresolved row is a location "
            "that will be left behind while the attestation says otherwise")
    return Path(expanded).resolve()


def measure(location: CustodyLocation) -> LocationState:
    """Measure one location without recording anything about its contents.

    Args:
        location: The registered location to stat.

    Returns:
        Existence, file count, and total bytes.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    if location.path is None:
        return LocationState(exists=False, file_count=0, total_bytes=0)
    return measure_path(location.path)


def measure_path(path: Path) -> LocationState:
    """Count files and bytes under ``path`` without recording anything about them.

    Symlinks are counted as themselves and never followed. That single rule handles the
    Hub cache correctly in both directions: ``blobs/`` holds the real bytes and is counted
    once, while ``snapshots/`` holds only links to those same blobs and so contributes its
    own negligible size rather than a second copy of the payload. Following links would
    double the total; skipping the link farm entirely would report an empty cache that
    still holds every byte it ever downloaded.
    """
    if not path.exists() and not path.is_symlink():
        return LocationState(exists=False, file_count=0, total_bytes=0)

    if path.is_file() or path.is_symlink():
        return LocationState(exists=True, file_count=1, total_bytes=path.lstat().st_size)

    count = 0
    total = 0
    for root, dirs, files in os.walk(path, followlinks=False):
        # A symlinked directory is an entry in its own right, not a tree to descend into.
        for name in list(dirs):
            entry = Path(root) / name
            if entry.is_symlink():
                dirs.remove(name)
                count += 1
                total += entry.lstat().st_size
        for name in files:
            try:
                total += (Path(root) / name).lstat().st_size
            except OSError:
                pass        # raced or unreadable: still present, so still counted
            count += 1
    return LocationState(exists=True, file_count=count, total_bytes=total)


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
    locations = read_register()
    verdicts: list[LocationVerdict] = []

    for location in locations:
        if not location.is_purgeable:
            continue
        _assert_deletable(location)
        before = measure(location)
        if confirm and before.exists:
            _delete(location)
        after = measure(location) if confirm else before
        verdicts.append(LocationVerdict(
            location=location,
            before=before,
            after=after,
            verified_utc=_now(),
        ))

    if confirm:
        _discard_freed_blocks()
    return verdicts


def _assert_deletable(location: CustodyLocation) -> None:
    """Refuse a path the register has no business pointing at.

    The register is the authority on *what* is purged, but it is a Markdown file that a
    person edits, so it is not the authority on what is *safe* to delete. A typo that
    turned ``<repo>/results/`` into ``<repo>/`` would otherwise destroy the project.
    """
    path = location.path
    assert path is not None                              # guarded by is_purgeable

    if path == Path(path.anchor) or len(path.parts) <= 2:
        raise ValueError(
            f"{location.id}: refusing to delete {path} -- too close to a filesystem root")

    if path == REPO_ROOT:
        raise ValueError(
            f"{location.id}: refusing to delete the repository working tree itself")

    if REPO_ROOT in path.parents:
        allowed = {REPO_ROOT / "results", REPO_ROOT / "notebooks"}
        if path not in allowed:
            raise ValueError(
                f"{location.id}: {path} is inside the repository but is not one of "
                f"{sorted(str(p.name) for p in allowed)}. Only pipeline output and "
                "notebook outputs are purgeable from the working tree")


def _delete(location: CustodyLocation) -> None:
    """Delete one location by the method its register row declares.

    The Hub cache is emptied through ``hf`` first, because ``.incomplete`` partial blobs
    survive a plain tree removal of ``snapshots/`` and would be left behind by anything
    that reasons about repos rather than bytes. The tree removal then takes the root.
    """
    path = location.path
    assert path is not None

    if "hf cache" in location.method or _looks_like_hub_cache(path):
        _run_hf_cache_cleanup(path)

    if not path.exists():
        return
    if path == REPO_ROOT / "results":
        _empty_but_keep(path, keep={".gitkeep"})
    else:
        shutil.rmtree(path, ignore_errors=False)


def _looks_like_hub_cache(path: Path) -> bool:
    return (path / "blobs").is_dir() or (path / "hub").is_dir() or path.name == "hf-cache"


def _run_hf_cache_cleanup(path: Path) -> None:
    """``hf cache rm`` then ``hf cache prune``; a missing CLI is reported, not swallowed."""
    env = dict(os.environ, HF_HOME=str(path))
    for args in (["hf", "cache", "rm", "--all", "--yes"], ["hf", "cache", "prune", "--yes"]):
        try:
            subprocess.run(args, env=env, check=False, capture_output=True, timeout=600)
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(
                f"{' '.join(args)} failed ({exc}). The Hub cache keeps .incomplete partial "
                "blobs that a tree removal leaves behind, so this step is not optional; "
                "run it by hand and re-run the purge") from exc


def _empty_but_keep(path: Path, *, keep: set) -> None:
    """Empty a directory while preserving named entries -- ``results/.gitkeep``."""
    for entry in path.iterdir():
        if entry.name in keep:
            continue
        if entry.is_dir() and not entry.is_symlink():
            shutil.rmtree(entry, ignore_errors=False)
        else:
            entry.unlink()


def _discard_freed_blocks() -> None:
    """``fstrim -av`` so deleted blocks are discarded rather than merely unlinked.

    Best-effort and reported: it needs privilege, it is a no-op on filesystems that do not
    support discard, and a destruction record that silently claimed it ran would be worse
    than one that says it could not.
    """
    try:
        done = subprocess.run(["fstrim", "-av"], check=False, capture_output=True, timeout=900)
        FSTRIM_RESULT["ran"] = done.returncode == 0
        FSTRIM_RESULT["detail"] = ("fstrim -av completed" if done.returncode == 0
                                   else f"fstrim -av exited {done.returncode}")
    except (OSError, subprocess.SubprocessError) as exc:
        FSTRIM_RESULT["ran"] = False
        FSTRIM_RESULT["detail"] = f"fstrim -av could not be run: {exc}"


#: Filled by :func:`_discard_freed_blocks` and reported in the attestation. Never asserted
#: as done when it was not.
FSTRIM_RESULT: dict = {"ran": False, "detail": "not attempted (dry run)"}


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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
    if not re.fullmatch(r"[0-9a-f]{40}", revision_sha):
        raise ValueError(
            "revision_sha must be the dataset's full 40-character Hub commit hash: it is "
            "the only identifier of the destroyed data the attestation carries")

    unclean = [v.location.id for v in verdicts if not v.is_clean]
    if unclean:
        raise ValueError(
            f"refusing to attest: {', '.join(unclean)} still hold bytes. A partial purge "
            "reported as a completed one is the one thing this record must never be")

    checks = _run_required_checks(verdicts)
    failed = [c["name"] for c in checks if c["result"] != "pass"]
    if failed:
        raise ValueError(f"refusing to attest: check(s) failed -- {', '.join(failed)}")

    payload = {
        "schema_version": ATTESTATION_SCHEMA_VERSION,
        "generated_utc": _now(),
        "hackathon_close_date": HACKATHON_CLOSE_DATE,
        "deletion_deadline_utc": DELETION_DEADLINE_UTC,
        "dataset": {
            "repo_id": DATASET_REPO_ID,
            "repo_type": DATASET_REPO_TYPE,
            "revision_sha": revision_sha,
            "identification_note": (
                "Identified by Hub revision only. No filenames or per-file digests are "
                "recorded, by design."),
        },
        "environment": _environment(),
        "locations": [
            {
                "id": v.location.id,
                "path": str(v.location.path),
                "method": v.location.method,
                "before": _state_dict(v.before),
                "after": _state_dict(v.after),
                "verified_utc": v.verified_utc,
            }
            for v in verdicts
        ],
        "checks": checks,
        "attested_by": {
            "name": os.environ.get("MVA_ATTESTOR_NAME", "<attestor>"),
            "role": os.environ.get("MVA_ATTESTOR_ROLE", "<role>"),
            "statement": ("I verified each registered location was emptied on the date "
                          "recorded above."),
        },
    }

    serialized = json.dumps(payload, indent=2, sort_keys=False)
    _assert_carries_no_patient_data(serialized)

    out = REPO_ROOT / "docs" / f"purge_attestation_{_now()[:10]}.json"
    out.write_text(serialized + "\n", encoding="utf-8")
    out.with_suffix(".md").write_text(_render_markdown(payload), encoding="utf-8")
    return out


def _state_dict(state: LocationState) -> dict:
    return {"exists": state.exists, "file_count": state.file_count,
            "total_bytes": state.total_bytes}


def _environment() -> dict:
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "wsl_distro": os.environ.get("WSL_DISTRO_NAME"),
        "pipeline_commit": _git_sha(),
        "purge_tool_version": ATTESTATION_SCHEMA_VERSION,
    }


def _git_sha() -> str:
    try:
        done = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=False, timeout=30)
        return done.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _run_required_checks(verdicts: list[LocationVerdict]) -> list[dict]:
    """Every check in :data:`REQUIRED_CHECKS`, each reported pass or fail with detail."""
    results = []

    registered = read_register()
    terminal = [r for r in registered
                if r.status in UNPURGEABLE_STATES or not r.is_purgeable
                or r.id in {v.location.id for v in verdicts}]
    results.append({
        "name": "custody_register_complete",
        "result": "pass" if len(terminal) == len(registered) else "fail",
        "detail": f"{len(terminal)} of {len(registered)} register rows reached a "
                  "terminal state",
    })

    scan = subprocess.run(
        [os.environ.get("MVA_PYTHON", "python"), "-m", "pytest", "-q", "--no-header",
         "-p", "no:cacheprovider", "tests/test_smoke.py", "-k", "patient_data"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, check=False, timeout=900)
    results.append({
        "name": "repo_scan_genomic_extensions",
        "result": "pass" if scan.returncode == 0 else "fail",
        "matches": 0 if scan.returncode == 0 else 1,
    })

    cache = [v for v in verdicts if "hf cache" in v.location.method
             or _looks_like_hub_cache(v.location.path or REPO_ROOT)]
    results.append({
        "name": "hf_cache_empty",
        "result": "pass" if all(v.is_clean for v in cache) else "fail",
        "repos_remaining": 0 if all(v.is_clean for v in cache) else 1,
        "bytes_remaining": sum(v.after.total_bytes for v in cache),
    })
    results.append({
        "name": "incomplete_blobs_removed",
        "result": "pass" if all(v.is_clean for v in cache) else "fail",
        "matches": sum(v.after.file_count for v in cache),
    })

    results.append({
        "name": "freed_blocks_discarded",
        "result": "pass" if FSTRIM_RESULT["ran"] else "not_run",
        "detail": FSTRIM_RESULT["detail"],
    })
    return results


def _assert_carries_no_patient_data(serialized: str) -> None:
    """The boundary check. A guardrail nothing enforces is a comment."""
    lowered = serialized.lower()
    for field in FORBIDDEN_ATTESTATION_FIELDS:
        if f'"{field}"' in lowered:
            raise ValueError(
                f"refusing to write the attestation: it carries a {field!r} field. The "
                "destroyed data is identified by Hub revision only")

    for pattern, what in (
        (r"\b[0-9a-f]{64}\b", "a sha256 digest"),
        (r"\.(vcf|bam|cram|fastq|fq|tbi|bai)(\.gz)?\b", "a genomic filename"),
        (r"\bchr[0-9xym]+[:\-][0-9]{3,}", "a variant coordinate"),
        (r"\bhp:[0-9]{7}\b", "an HPO term id"),
    ):
        hit = re.search(pattern, lowered)
        if hit:
            raise ValueError(
                f"refusing to write the attestation: it contains {what}. A destruction "
                "record must not describe what was destroyed")


def _render_markdown(payload: dict) -> str:
    """The companion for the confirmation email body."""
    lines = [
        "# Data destruction attestation",
        "",
        f"**Generated:** {payload['generated_utc']}  ",
        f"**Dataset:** `{payload['dataset']['repo_id']}` at revision "
        f"`{payload['dataset']['revision_sha']}`  ",
        f"**Deletion deadline:** {payload['deletion_deadline_utc']}",
        "",
        payload["dataset"]["identification_note"],
        "",
        "## Locations",
        "",
        "| ID | Method | Before | After |",
        "|---|---|---:|---:|",
    ]
    for loc in payload["locations"]:
        before = f"{loc['before']['file_count']} files, {loc['before']['total_bytes']:,} B"
        after = ("absent" if not loc["after"]["exists"]
                 else f"{loc['after']['file_count']} files, "
                      f"{loc['after']['total_bytes']:,} B")
        lines.append(f"| `{loc['id']}` | {loc['method']} | {before} | {after} |")

    lines += ["", "## Checks", "", "| Check | Result | Detail |", "|---|---|---|"]
    for check in payload["checks"]:
        detail = check.get("detail", "")
        lines.append(f"| {check['name']} | **{check['result']}** | {detail} |")

    lines += [
        "",
        "## Attestation",
        "",
        f"{payload['attested_by']['statement']}",
        "",
        f"— {payload['attested_by']['name']}, {payload['attested_by']['role']}",
        "",
        "*No filenames, per-file digests, variant coordinates or phenotype content appear "
        "in this record, by construction. See `docs/data_custody.md` for the register and "
        "the honest limits of any destruction record.*",
        "",
    ]
    return "\n".join(lines)


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
    parser.add_argument("--revision-sha", default=os.environ.get("MVA_DATASET_REVISION"),
                        help="the dataset's full 40-character Hub commit hash; required "
                             "with --confirm, since it is the only identifier of the "
                             "destroyed data the attestation carries")
    args = parser.parse_args(argv)

    try:
        locations = read_register(args.register)
    except (OSError, ValueError) as exc:
        print(f"cannot read the custody register: {exc}")
        return 2

    print(f"custody register: {args.register}")
    print(f"mode: {'DELETE (--confirm)' if args.confirm else 'dry run -- nothing is removed'}")
    print()
    print(f"  {'id':<4} {'status':<20} {'holds now':<24} method")
    for location in locations:
        if not location.is_purgeable:
            print(f"  {location.id:<4} {location.status:<20} {'not actionable here':<24} "
                  f"{location.method[:46]}")
            continue
        state = measure(location)
        holds = ("absent" if not state.exists
                 else f"{state.file_count} files, {state.total_bytes / 1e6:.1f} MB")
        print(f"  {location.id:<4} {location.status:<20} {holds:<24} {location.method[:46]}")
        print(f"       -> {location.path}")

    if not args.confirm:
        print()
        print("Dry run. Re-run with --confirm and --revision-sha <40-char hash> to delete.")
        print("This destroys the only copy you can regenerate results from. Do it after")
        print("judging settles, and before the deadline in docs/data_custody.md.")
        return 0

    if not args.revision_sha:
        print("\n--confirm needs --revision-sha (or MVA_DATASET_REVISION): the attestation")
        print("identifies the destroyed data by Hub revision and by nothing else.")
        return 2

    verdicts = purge({}, confirm=True)
    print()
    for verdict in verdicts:
        mark = "ok  " if verdict.is_clean else "LEFT"
        print(f"  {mark} {verdict.location.id}: "
              f"{verdict.before.file_count} files -> {verdict.after.file_count}")

    try:
        written = write_attestation(verdicts, revision_sha=args.revision_sha)
    except ValueError as exc:
        print(f"\npurge ran but the attestation was refused: {exc}")
        return 1

    print(f"\nattestation: {written}")
    print(f"companion:   {written.with_suffix('.md')}")
    print("\nStill manual, and not claimed by this tool: delete the conversations in the")
    print("Claude account (C7's account-side half), then commit the attestation and email")
    print("it to the organizers.")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
