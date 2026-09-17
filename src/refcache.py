"""Cache for public reference files, shared by any layer that needs one.

Purpose
    Download a public reference file once, keep it outside the repository, and return a
    provenance record -- URL, retrieval date, byte size and SHA-256 -- so a result can be
    traced to the exact bytes behind it.

Inputs
    ``config["reference_dir"]`` (or ``$MVA_REF_ROOT``) and a URL.

Outputs
    The cached file's path and its provenance record.

Guardrail
    **Whole public files only, never a query.** Every caller downloads a complete public
    release and matches against it locally. A per-variant lookup against a remote API
    would put the subject's coordinates on a third-party server, which ``COMPLIANCE.md``
    prohibits -- the same reason L0 annotates with a local snpEff rather than the VEP REST
    API (``l0_genomics/annotate.py``).

    Cached files are public reference data, not patient data: they live in
    ``reference_dir``, outside the repository and outside the custody root, and must never
    enter the deletion register in ``docs/data_custody.md``. "Outside the repository" is
    enforced, not assumed -- one of these releases is a VCF, and the rule that no genomic
    file sits in the working tree stays absolute.

    This module never opens ``data_dir`` and never writes under ``results_dir``.
"""

from __future__ import annotations

import hashlib
import os
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
USER_AGENT = {"User-Agent": "mva-track2/1.0 (https://github.com/KFlowerz/HKTN_RareDisease_MVA)"}
#: Used when nothing is configured. Outside the repository -- see :func:`reference_dir`.
DEFAULT_DIR = "~/.cache/mva-track2/reference"


def cache_dir(raw, default: str, *, name: str = "cache") -> Path:
    """Resolve a cache location, refusing anywhere inside the repository.

    **A cache must not sit inside the working tree.** One cached release is ClinVar's
    GRCh38 VCF, and ``tests/test_smoke.py`` forbids a ``.vcf.gz`` anywhere in the repo --
    a rule worth keeping absolute rather than carving an exception into, since the next
    genomic file to land there may not be a public one. Another cache holds
    licence-restricted data that must never be committed. A relative path is therefore
    resolved against the user's home directory, not the repo.

    Raises:
        ValueError: If the resolved directory lies inside the repository.
    """
    path = Path(raw or default).expanduser()
    if not path.is_absolute():
        path = Path.home() / path
    path = path.resolve()
    if path == REPO_ROOT or REPO_ROOT in path.parents:
        raise ValueError(
            f"{name} {path} is inside the repository. Cached releases are never committed "
            "-- some are genomic files, some are licence-restricted -- so point it outside "
            "the working tree."
        )
    path.mkdir(parents=True, exist_ok=True)
    return path


def reference_dir(config: dict) -> Path:
    """Where cached public reference files live. ``$MVA_REF_ROOT`` overrides the config."""
    return cache_dir(os.environ.get("MVA_REF_ROOT") or config.get("reference_dir"),
                     DEFAULT_DIR, name="reference_dir")


def sha256(path: Path) -> str:
    """SHA-256 of a file, read in chunks so a 100 MB release costs no more memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, dest: Path) -> dict:
    """Download once and reuse the cached copy; return a provenance record.

    The partial file is written beside the target and renamed on completion, so an
    interrupted download cannot be mistaken for a complete one on the next run.
    """
    if not dest.exists():
        part = dest.with_name(dest.name + ".part")
        request = urllib.request.Request(url, headers=USER_AGENT)
        with urllib.request.urlopen(request, timeout=300) as response, open(part, "wb") as handle:
            while chunk := response.read(1 << 20):
                handle.write(chunk)
        part.rename(dest)
        retrieved = date.today().isoformat()
    else:
        retrieved = datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc).date().isoformat()

    return {"url": url, "file": dest.name, "bytes": dest.stat().st_size,
            "sha256": sha256(dest), "retrieved": retrieved}
