"""L1 -- interactome and pathway sources: download, cache, parse.

Purpose
    Fetch the public resources the module is built from, cache them outside the
    repository, and parse them into the plain structures :mod:`src.l1_target.module`
    consumes.

Inputs
    ``config["l1"]`` (``string_version``, ``taxon``, ``species``), ``config["reference_dir"]``
    via :mod:`src.refcache`, and network access to STRING and Reactome.

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

    Caching and its provenance record live in :mod:`src.refcache`, which states the
    whole-file-never-a-query rule every layer's downloads obey.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from ..refcache import fetch, reference_dir

STRING_BASE = "https://stringdb-downloads.org/download"
REACTOME_URL = "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
LICENCES = {"STRING": "CC BY 4.0", "Reactome": "CC0 1.0"}
#: UniProt accession syntax, so alias rows that are not accessions are ignored.
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def string_urls(version: str, taxon: int) -> dict:
    """The four STRING files this layer uses, by role."""
    return {
        "links": f"{STRING_BASE}/protein.links.v{version}/{taxon}.protein.links.v{version}.txt.gz",
        "physical": f"{STRING_BASE}/protein.physical.links.v{version}/{taxon}.protein.physical.links.v{version}.txt.gz",
        "info": f"{STRING_BASE}/protein.info.v{version}/{taxon}.protein.info.v{version}.txt.gz",
        "aliases": f"{STRING_BASE}/protein.aliases.v{version}/{taxon}.protein.aliases.v{version}.txt.gz",
    }


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
