"""LINCS L1000 access -- the reference data behind Channel C.

Purpose
    Fetch, cache and read the two public LINCS L1000 releases Channel C scores against:
    the consensus gene-knockdown signatures (GEO ``GSE106127``) that supply the **proxy**
    disease signature, and the Level 5 compound signatures (GEO ``GSE70138``, LINCS Phase
    II) that supply the perturbations to score.

Inputs
    ``config["reference_dir"]`` (or ``$MVA_REF_ROOT``). Nothing else -- the release
    identity is fixed in :data:`RELEASES`, not configured, because a signature matrix and
    its own metadata file are one artifact and must not be mixed across releases.

Outputs
    Cached file paths plus a provenance record per file (URL, accession, bytes, SHA-256,
    retrieval date), and typed readers: metadata rows, gene identifiers, and blocks of
    the signature matrix.

Guardrail
    **Whole public releases only, never a query** -- the rule :mod:`src.refcache` exists
    to enforce. Nothing about the subject reaches this module: the only thing Channel C
    asks LINCS for is the *causal gene symbol*, which is a public research premise
    recorded at gate G1, not patient data. This module never opens ``data_dir`` and never
    writes under ``results_dir``.

    **The GCTX axes are transposed relative to their names, and that is the single most
    error-prone fact here.** In a ``.gctx`` file, ``0/META/ROW/id`` names the *genes* and
    ``0/META/COL/id`` names the *signatures*, but the HDF5 matrix is stored column-major,
    so ``matrix[i, j]`` reads as ``[signature i, gene j]``. Reading it the other way
    silently produces a matrix of the right dtype and the wrong meaning, which no
    downstream check would catch. :func:`read_block` is the only sanctioned way in, and
    :func:`axes` asserts the orientation every time a file is opened.

    LINCS files are public reference data, not patient data. They live in
    ``reference_dir``, outside the repository and outside the custody root, and must
    never enter the deletion register in ``docs/data_custody.md``.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import shutil
from pathlib import Path

from .. import refcache

LOGGER = logging.getLogger(__name__)

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
#: Subdirectory of ``reference_dir``; these releases are large and worth keeping apart.
CACHE_SUBDIR = "lincs"

#: The two GEO series, pinned. ``knockdown`` supplies the proxy signature and
#: ``compound`` the perturbations scored against it. Not configurable: a matrix and its
#: metadata are one artifact, and letting a config key move one without the other would
#: silently mismatch signature ids to rows.
RELEASES = {
    "knockdown": {
        "accession": "GSE106127",
        "prefix": f"{GEO_FTP}/GSE106nnn/GSE106127/suppl",
        "matrix": "GSE106127_CGS_n33839x978.gctx",
        "meta": "GSE106127_CGS_meta.txt.gz",
        "genes": "GSE106127_gene_info.txt.gz",
        "description": "LINCS L1000 consensus gene signatures (CGS) of shRNA knockdown",
    },
    "compound": {
        "accession": "GSE70138",
        "prefix": f"{GEO_FTP}/GSE70nnn/GSE70138/suppl",
        "matrix": "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
        "meta": "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "genes": "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
        # The release's pert_info table is deliberately not fetched: sig_info already
        # carries pert_iname, and a file nothing reads is a file whose licence and
        # provenance still have to be defended.
        "description": "LINCS L1000 Phase II (GSE70138) Level 5 compound signatures",
    },
}

#: LINCS' own perturbation-type vocabulary. ``trt_sh.cgs`` is a knockdown *consensus*
#: across the shRNAs targeting one gene, not a single hairpin; ``trt_cp`` is a compound
#: treatment. Controls and other perturbation classes are never scored.
KNOCKDOWN_TYPE = "trt_sh.cgs"
COMPOUND_TYPE = "trt_cp"

#: HDF5 paths inside a ``.gctx``. Fixed by the GCTX 1.0 specification.
MATRIX_PATH = "0/DATA/0/matrix"
GENE_ID_PATH = "0/META/ROW/id"
SIGNATURE_ID_PATH = "0/META/COL/id"

#: What the licence check found, recorded with the data rather than assumed. See
#: ``src/l4_validate/sources.md`` Table 10.
LICENCE = ("NCBI GEO supplementary files for GSE106127 and GSE70138 (NIH LINCS Program, "
           "Broad Institute). No machine-readable licence statement was resolvable at the "
           "LINCS project URL on 2026-09-21; treated conservatively as non-redistributable "
           "and matched locally, with only derived scores and identifiers published.")
REDISTRIBUTABLE = False


def cache_path(config: dict) -> Path:
    """``reference_dir/lincs`` -- created on demand, always outside the repository."""
    path = refcache.reference_dir(config) / CACHE_SUBDIR
    path.mkdir(parents=True, exist_ok=True)
    return path


def _checksum(path: Path) -> str:
    """SHA-256 of a cached file, with a sidecar so a 5 GB matrix is hashed once.

    The sidecar records the size and mtime it was computed for and is ignored if either
    moved, so an edited or re-downloaded file is never reported under a stale digest.
    """
    sidecar = path.with_name(path.name + ".sha256.json")
    stat = path.stat()
    if sidecar.is_file():
        try:
            cached = json.loads(sidecar.read_text(encoding="utf-8"))
            if cached.get("bytes") == stat.st_size and cached.get("mtime") == int(stat.st_mtime):
                return str(cached["sha256"])
        except (ValueError, KeyError):
            pass  # unreadable sidecar is not worth failing over; recompute below
    LOGGER.info("hashing %s (%.1f GB) -- cached for later runs", path.name,
                stat.st_size / 1e9)
    digest = refcache.sha256(path)
    sidecar.write_text(json.dumps({"sha256": digest, "bytes": stat.st_size,
                                   "mtime": int(stat.st_mtime)}), encoding="utf-8")
    return digest


def _ensure_plain(url: str, cache: Path, name: str) -> tuple:
    """Cache a file GEO serves gzipped but that has to be read uncompressed.

    A ``.gctx`` is HDF5 and cannot be read through a gzip stream, so the release is
    downloaded as ``<name>.gz`` and expanded once. The compressed copy is removed
    afterwards -- these are 5 GB files and keeping both doubles the cache for nothing --
    but its SHA-256 is recorded first, so the provenance names the bytes GEO actually
    served rather than only the local expansion.

    Expansion goes to a ``.part`` and is renamed on completion, so an interrupted run
    cannot leave a truncated matrix that looks finished.
    """
    dest = cache / name
    if dest.is_file():
        return dest, {"file": name, "bytes": dest.stat().st_size,
                      "sha256": _checksum(dest), "retrieved": _mtime_date(dest),
                      "url": url,
                      "note": "already cached; SHA-256 is of the expanded file, because "
                              "the GEO .gz it came from is not kept"}

    archive = cache / f"{name}.gz"
    record = refcache.fetch(url, archive)
    LOGGER.info("expanding %s", archive.name)
    part = dest.with_name(dest.name + ".part")
    try:
        with gzip.open(archive, "rb") as source, open(part, "wb") as handle:
            shutil.copyfileobj(source, handle, length=1 << 22)
        part.replace(dest)
    finally:
        part.unlink(missing_ok=True)
    archive.unlink(missing_ok=True)
    return dest, {**record, "file": name, "expanded_bytes": dest.stat().st_size,
                  "note": "SHA-256 and byte count are of the GEO .gz, which is removed "
                          "after expansion; the expanded file is what is read"}


def _mtime_date(path: Path) -> str:
    from datetime import datetime, timezone

    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).date().isoformat()


def ensure(config: dict) -> tuple:
    """Cache every LINCS file Channel C needs. Returns ``(paths, provenance)``.

    ``paths`` is ``{release: {role: Path}}`` over the roles in :data:`RELEASES`. Metadata
    files stay gzipped and are read through :func:`read_metadata`; the two matrices are
    expanded, because HDF5 cannot be read from a compressed stream.
    """
    cache = cache_path(config)
    paths, provenance = {}, []
    for release, spec in RELEASES.items():
        paths[release] = {}
        for role, name in spec.items():
            if role in ("accession", "prefix", "description"):
                continue
            url = f"{spec['prefix']}/{name}" if role != "matrix" else f"{spec['prefix']}/{name}.gz"
            if role == "matrix":
                path, record = _ensure_plain(url, cache, name)
            else:
                path = cache / name
                record = refcache.fetch(url, path)
            paths[release][role] = path
            provenance.append({"source": "NCBI GEO", "accession": spec["accession"],
                               "release": spec["description"], "role": role,
                               "licence": LICENCE, "redistributable": REDISTRIBUTABLE,
                               **record})
    return paths, provenance


def read_metadata(path: Path) -> list:
    """Rows of a gzipped LINCS metadata table, as dicts keyed by its header."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def landmark_ids(gene_info: list) -> list:
    """Entrez ids of the 978 landmark genes, in the order the gene table lists them.

    The landmark space is the only one both releases share and the only one L1000
    *measures*: everything else in the Phase II matrix is inferred from it, so scoring a
    query over inferred genes would count the landmarks twice under other names.

    Raises:
        ValueError: If the table carries no ``pr_is_lm`` column, or marks no landmark.
    """
    if not gene_info or "pr_is_lm" not in gene_info[0]:
        raise ValueError("gene table has no pr_is_lm column; without it the landmark "
                         "space cannot be identified and the rest is inferred expression")
    ids = [row["pr_gene_id"] for row in gene_info if str(row.get("pr_is_lm", "")).strip() == "1"]
    if not ids:
        raise ValueError("gene table marks no landmark gene (pr_is_lm == 1)")
    return ids


def axes(handle) -> tuple:
    """``(signature_ids, gene_ids)`` of an open ``.gctx``, with the orientation checked.

    Returns the ids as lists of ``str``; HDF5 stores them as fixed-width bytes.

    Raises:
        ValueError: If the matrix shape does not match ``(len(COL/id), len(ROW/id))``.
            That is the transposition trap in this module's guardrail: a file read the
            wrong way round still yields numbers, and this is the only place it can be
            caught.
    """
    gene_ids = [x.decode() if isinstance(x, bytes) else str(x)
                for x in handle[GENE_ID_PATH][:]]
    signature_ids = [x.decode() if isinstance(x, bytes) else str(x)
                     for x in handle[SIGNATURE_ID_PATH][:]]
    shape = handle[MATRIX_PATH].shape
    if shape != (len(signature_ids), len(gene_ids)):
        raise ValueError(
            f"matrix is {shape} but the file names {len(signature_ids)} signature(s) and "
            f"{len(gene_ids)} gene(s). GCTX stores the matrix transposed relative to its "
            "axis names; a mismatch here means the layout is not what this reader assumes"
        )
    return signature_ids, gene_ids


def read_block(handle, start: int, stop: int, gene_columns):
    """``(stop - start, len(gene_columns))`` z-scores -- signatures by genes.

    Args:
        handle: an open :class:`h5py.File` over a ``.gctx``.
        start: first signature index, inclusive.
        stop: last signature index, exclusive.
        gene_columns: gene indices, ascending -- here, the landmark subset.

    Deliberately a **contiguous slice** rather than a list of signature indices. The
    signatures Channel C wants are scattered across all 118,050 rows, and HDF5 would
    serve a scattered fancy index by reading everything between the first and the last
    anyway. Handing the caller slices makes that cost explicit and bounds peak memory at
    one block of every gene, instead of hiding a whole-matrix read behind an index array
    that merely looked small.
    """
    import numpy as np

    gene_columns = np.asarray(gene_columns)
    if stop <= start:
        return np.empty((0, gene_columns.size), dtype=np.float64)
    block = handle[MATRIX_PATH][start:stop, :]
    return np.asarray(block[:, gene_columns], dtype=np.float64)
