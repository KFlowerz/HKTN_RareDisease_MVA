"""L1 -- interactome and pathway sources: download, cache, parse.

Purpose
    Fetch the public resources the module is built from, cache them outside the
    repository, and parse them into the plain structures :mod:`src.l1_target.module`
    consumes.

Inputs
    ``config["l1"]`` (``string_version``, ``taxon``, ``species``, ``reference_dir``) and
    network access to STRING and Reactome.

Outputs
    Parsed edges, symbol maps, pathway memberships, and a provenance record -- URL,
    licence, retrieval date, byte size and SHA-256 -- for every file used, so a result can
    be traced to the exact bytes behind it.

Guardrail
    **Licence-clean only.** STRING [szklarczyk2022] doi:10.1093/nar/gkac1000 is CC BY 4.0
    and Reactome data [milacic2023] doi:10.1093/nar/gkad1025 is CC0, so both may feed a
    redistributed CC-BY-4.0 output path; the attribution STRING asks for travels in the
    provenance record. Nothing carrying NonCommercial or ShareAlike terms may be added
    here -- ShareAlike would relicense every derived table (see
    ``../l4_validate/sources.md``).

    **No patient data reaches this module.** It never opens ``data_dir``. Downloads are
    whole public files with no query attached, so not even a gene symbol leaves the
    machine.

    Cached files are public reference data, not patient data: they live in
    ``reference_dir``, outside the repository and outside the custody root, and must never
    enter the deletion register in ``docs/data_custody.md``.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
import urllib.request
from datetime import date, timezone, datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
USER_AGENT = {"User-Agent": "mva-track2-l1/1.0 (https://github.com/KFlowerz/HKTN_RareDisease_MVA)"}
STRING_BASE = "https://stringdb-downloads.org/download"
REACTOME_URL = "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
LICENCES = {"STRING": "CC BY 4.0", "Reactome": "CC0 1.0"}
#: UniProt accession syntax, so alias rows that are not accessions are ignored.
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def reference_dir(config: dict) -> Path:
    """Where cached source files live. ``MVA_REF_ROOT`` overrides the config."""
    raw = (os.environ.get("MVA_REF_ROOT")
           or (config.get("l1") or {}).get("reference_dir")
           or "reference")
    path = Path(raw)
    if not path.is_absolute():
        path = REPO_ROOT / path
    path.mkdir(parents=True, exist_ok=True)
    return path


def string_urls(version: str, taxon: int) -> dict:
    """The four STRING files this layer uses, by role."""
    return {
        "links": f"{STRING_BASE}/protein.links.v{version}/{taxon}.protein.links.v{version}.txt.gz",
        "physical": f"{STRING_BASE}/protein.physical.links.v{version}/{taxon}.protein.physical.links.v{version}.txt.gz",
        "info": f"{STRING_BASE}/protein.info.v{version}/{taxon}.protein.info.v{version}.txt.gz",
        "aliases": f"{STRING_BASE}/protein.aliases.v{version}/{taxon}.protein.aliases.v{version}.txt.gz",
    }


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

    digest = hashlib.sha256()
    with open(dest, "rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return {"url": url, "file": dest.name, "bytes": dest.stat().st_size,
            "sha256": digest.hexdigest(), "retrieved": retrieved}


def ensure_sources(config: dict) -> tuple:
    """Download (or reuse) every source file. Returns ``(paths, provenance)``."""
    l1 = config.get("l1") or {}
    version = str(l1.get("string_version", "12.0"))
    taxon = int(l1.get("taxon", 9606))
    cache = reference_dir(config)

    paths, provenance = {}, []
    for role, url in string_urls(version, taxon).items():
        dest = cache / url.rsplit("/", 1)[-1]
        paths[role] = dest
        provenance.append({"source": "STRING", "version": version,
                           "licence": LICENCES["STRING"], "role": role, **fetch(url, dest)})
    dest = cache / REACTOME_URL.rsplit("/", 1)[-1]
    paths["reactome"] = dest
    provenance.append({"source": "Reactome", "version": "current",
                       "licence": LICENCES["Reactome"], "role": "pathways",
                       **fetch(REACTOME_URL, dest)})
    return paths, provenance


def _open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def load_symbols(info_path: Path) -> tuple:
    """``(symbol -> protein, protein -> symbol)`` from STRING's protein-info file."""
    by_symbol, by_protein = {}, {}
    with _open(info_path) as handle:
        next(handle, None)
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            protein, name = fields[0], fields[1]
            by_symbol.setdefault(name, protein)
            by_protein[protein] = name
    return by_symbol, by_protein


def load_edges(path: Path, *, score_min: int):
    """Yield ``(a, b, score)`` for edges at or above ``score_min``."""
    with _open(path) as handle:
        next(handle, None)
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            score = int(fields[2])
            if score >= score_min:
                yield fields[0], fields[1], score


def load_uniprot(aliases_path: Path) -> dict:
    """``{UniProt accession: protein}`` from STRING's alias table.

    Isoform suffixes are stripped, because Reactome keys pathways on the accession.
    """
    out = {}
    with _open(aliases_path) as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or "UniProt" not in fields[2]:
                continue
            accession = fields[1].split("-")[0]
            if UNIPROT_RE.match(accession):
                out.setdefault(accession, fields[0])
    return out


def load_pathways(reactome_path: Path, uniprot_to_protein: dict, *, species: str) -> dict:
    """``{protein: {pathway id, ...}}`` for one species.

    Reactome ships every species in one table, so the species column is filtered here
    rather than trusted from the file name.
    """
    out = {}
    with _open(reactome_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] != species:
                continue
            protein = uniprot_to_protein.get(fields[0].split("-")[0])
            if protein:
                out.setdefault(protein, set()).add(fields[1])
    return out
