"""L3 -- Drug identity: collapse duplicates, then attach cross-system identifiers.

Purpose
    Give every candidate one identity before anything is aggregated. Two problems have to
    be solved in that order, and they are not the same problem:

    1. **The same molecule appears more than once.** Channel B's top five on the real data
       were bortezomib, bortezomib D-mannitol, ixazomib, ixazomib citrate and carfilzomib
       -- salt and co-crystal forms competing with their own parents for rank. Aggregating
       that way splits a compound's support across its own formulations and destroys the
       cross-channel convergence signal the whole design rests on.
    2. **The candidate needs an identifier the outside world uses.** RxNorm's RxCUI is the
       backbone the drug-annotation plan is built on (``CLAUDE.md``), and L4 joins its
       label-derived annotations on it.

    Problem 1 is solved locally and completely, by Open Targets' own ``parentId``. Problem
    2 is solved for as many candidates as the public sources cover, which is not all of
    them -- see the coverage note below.

Inputs
    The L2 channel tables; Open Targets ``drug_molecule`` for the parent relation
    (restricted zone, :mod:`..l2_channels.enrichment`); UniChem's whole ChEMBL-to-UNII
    mapping; openFDA's NDC and Drugs@FDA bulk releases for UNII-to-RxCUI.

Outputs
    One :class:`Identity` per distinct molecule: the parent ChEMBL id as its key, every
    ChEMBL id that collapsed into it, its RxCUI set where resolvable, how that resolution
    was made, and the openFDA pharmacologic class fields L4 will want.

Guardrail
    **Whole public files only, never a per-drug query.** RxNav would answer this in one
    call per drug, and this module deliberately does not use it. Channel D's candidates are
    derived from the subject's phenotype, so the *set* of drugs asked about is itself a weak
    statement about the child; a bulk file asked about nothing. The same rule as everywhere
    else in this pipeline, applied to a case where the convenient route is a query.

    **This module never reads patient-derived output.** It reads ``candidates.tsv`` per
    channel, which carries identities and counts. It must not open channel D's
    ``evidence.json``, which is patient-derived by construction.

    **The licence boundary is a directory boundary** (decision D7). The ChEMBL-to-UNII
    mapping is treated as ChEMBL-derived -- UniChem states no licence this project could
    read -- so it is cached in ``enrichment_dir`` and never redistributed. RxCUI, UNII and
    the pharmacologic classes come from openFDA, are US-Government public domain, and may
    be published. Only the second kind reaches an output.

    **An unresolved RxCUI is recorded, never dropped.** Roughly 45% of ChEMBL-sourced
    candidates have no RxCUI in these sources, mostly because they are approved outside
    the United States and so have no NDC. Dropping them would silently restrict the whole
    pipeline to the US market while appearing to be a technical detail.
"""

from __future__ import annotations

import gzip
import json
import logging
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

from .. import refcache
from ..l2_channels import enrichment

LOGGER = logging.getLogger(__name__)

#: ChEMBL (src 1) to FDA/USP Substance Registration System, i.e. UNII (src 14).
UNICHEM_URL = ("https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/"
               "src_id1/src1src14.txt.gz")
UNICHEM_LICENCE = ("UniChem ChEMBL-to-UNII mapping (EMBL-EBI). No licence statement was "
                   "readable on 2026-09-18 -- the FAQ is a script-only page and the FTP "
                   "README carries none -- so the mapping is treated as ChEMBL-derived "
                   "(CC BY-SA 3.0) and kept in the enrichment zone, per decision D7.")

#: openFDA bulk releases. US-Government public domain; no key; whole files.
OPENFDA_DOWNLOADS = "https://api.fda.gov/download.json"
OPENFDA_DATASETS = ("ndc", "drugsfda")
OPENFDA_LICENCE = ("openFDA (U.S. Food and Drug Administration): US-Government work, public "
                   "domain / CC0. openFDA states its data is not for clinical use and may "
                   "be incomplete or inaccurate.")

#: Fields taken from an openFDA record. All public domain, all safe to publish.
OPENFDA_FIELDS = ("rxcui", "unii", "pharm_class_epc", "pharm_class_moa")

#: How an RxCUI was arrived at, strongest first. Recorded per candidate, because a name
#: match and a registered-substance match are not the same kind of evidence.
RESOLUTION_ORDER = ("unii", "name", "unresolved")


@dataclass(frozen=True)
class Identity:
    """One molecule, after its formulations have been collapsed onto their parent."""

    key: str                       # parent ChEMBL id -- the aggregation key
    name: str = ""
    chembl_ids: frozenset = frozenset()   # every id that collapsed into this one
    rxcui: tuple = ()
    unii: tuple = ()
    pharm_class_epc: tuple = ()
    pharm_class_moa: tuple = ()
    resolution: str = "unresolved"
    drug_type: str = ""
    clinical_stage: str = ""

    @property
    def collapsed(self) -> int:
        """How many ChEMBL ids folded into this identity beyond the first."""
        return max(0, len(self.chembl_ids) - 1)

    def publishable(self) -> dict:
        """The fields that may leave the enrichment zone.

        ChEMBL id, name, type and stage are the Table 2 whitelist
        (:data:`..l2_channels.enrichment.PUBLISHABLE_FIELDS`); RxCUI, UNII and the
        pharmacologic classes are openFDA's and are public domain.
        """
        return {"chembl_id": self.key, "drug_name": self.name,
                "drug_type": self.drug_type, "clinical_stage": self.clinical_stage,
                "rxcui": "|".join(self.rxcui), "unii": "|".join(self.unii),
                "pharm_class_epc": "|".join(self.pharm_class_epc),
                "pharm_class_moa": "|".join(self.pharm_class_moa),
                "rxcui_resolution": self.resolution,
                "n_chembl_ids_collapsed": self.collapsed}


def normalise(value: str) -> str:
    """Casefold and strip, for name matching. Deliberately conservative."""
    return " ".join(str(value or "").strip().lower().split())


def parent_index(paths) -> tuple:
    """``({chembl_id: parent_chembl_id}, identities, stats)`` from Open Targets.

    ``parentId`` is Open Targets' own statement that one molecule is a form of another --
    a salt, a hydrate, a co-crystal. Using it rather than inferring from names means the
    collapse is somebody's curated assertion, not our string matching.
    """
    parent, identity, stats = {}, {}, {"molecules": 0, "with_a_parent": 0}
    for path in paths:
        for row in enrichment._rows(path):
            chembl_id = row.get("id")
            if not chembl_id:
                continue
            stats["molecules"] += 1
            identity[chembl_id] = {"name": row.get("name") or "",
                                   "drug_type": row.get("drugType") or "",
                                   "clinical_stage": row.get("maximumClinicalStage") or ""}
            if row.get("parentId"):
                parent[chembl_id] = row["parentId"]
                stats["with_a_parent"] += 1
    return parent, identity, stats


def resolve_parent(chembl_id: str, parent: dict, *, limit: int = 10) -> str:
    """Follow ``parentId`` to the root molecule.

    ``limit`` stops a cycle in the source data from hanging the pipeline; a chain that
    long is a data problem, and returning the last id reached keeps it visible rather than
    raising in the middle of a run.
    """
    seen = {chembl_id}
    current = chembl_id
    for _ in range(limit):
        nxt = parent.get(current)
        if not nxt or nxt in seen:
            break
        seen.add(nxt)
        current = nxt
    return current


def load_chembl_to_unii(path: Path) -> dict:
    """``{chembl_id: {unii}}`` from UniChem's whole mapping file.

    Restricted content: this stays in the enrichment zone and never reaches an output.
    """
    out: dict = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if "src" not in header.lower():
            raise ValueError(f"{path} does not look like a UniChem mapping file; its first "
                             f"line was {header[:80]!r}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] and fields[1]:
                out.setdefault(fields[0], set()).add(fields[1])
    return out


def _openfda_records(path: Path):
    """Yield every record from an openFDA bulk zip."""
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if not name.endswith(".json"):
                continue
            payload = json.loads(archive.read(name))
            yield from payload.get("results", [])


def load_openfda(paths) -> tuple:
    """``(by_unii, by_name, stats)`` -- RxCUI and class fields keyed two ways.

    ``by_name`` is the fallback route and is kept separate from ``by_unii`` on purpose, so
    a candidate resolved by string match can be reported as such. A name is not a
    registered identifier: two products can share a generic name across salt forms, and a
    brand name can be reused. It resolves candidates that would otherwise be lost, and the
    output says which ones it resolved.
    """
    by_unii: dict = {}
    by_name: dict = {}
    stats = {"openfda_records": 0, "records_with_rxcui": 0}

    def absorb(target: dict, key: str, block: dict) -> None:
        entry = target.setdefault(key, {f: set() for f in OPENFDA_FIELDS})
        for field_name in OPENFDA_FIELDS:
            entry[field_name].update(block.get(field_name) or [])

    for path in paths:
        for record in _openfda_records(path):
            stats["openfda_records"] += 1
            block = record.get("openfda") or {}
            if not block.get("rxcui"):
                continue
            stats["records_with_rxcui"] += 1
            for unii in block.get("unii") or []:
                absorb(by_unii, unii, block)
            names = [record.get("generic_name"), record.get("brand_name")]
            names += [p.get("brand_name") for p in record.get("products") or []]
            for name in names:
                key = normalise(name)
                if key:
                    absorb(by_name, key, block)

    stats["unii_keys"] = len(by_unii)
    stats["name_keys"] = len(by_name)
    return by_unii, by_name, stats


def ensure_sources(config: dict) -> tuple:
    """Download (or reuse) the crosswalk files. Returns ``(paths, provenance)``.

    The two caches are the licence boundary: UniChem's ChEMBL-derived mapping goes to
    ``enrichment_dir``, openFDA's public-domain releases to ``reference_dir``.
    """
    enrichment_cache = enrichment.enrichment_dir(config)
    reference_cache = refcache.reference_dir(config)

    unichem_dest = enrichment_cache / "unichem.src1src14.txt.gz"
    provenance = [{"source": "UniChem", "dataset": "src1src14 (ChEMBL to UNII)",
                   "licence": UNICHEM_LICENCE, "redistributable": False,
                   **refcache.fetch(UNICHEM_URL, unichem_dest)}]

    # The download index names the current export and its parts, so a re-run picks up a
    # new openFDA export rather than pinning to a URL that silently starts 404ing.
    index = _download_index()
    openfda_paths = []
    for dataset in OPENFDA_DATASETS:
        entry = index["results"]["drug"][dataset]
        export = str(entry.get("export_date") or "unknown")
        for number, part in enumerate(entry["partitions"]):
            # The export date is in the filename, so a new openFDA export lands as a new
            # file. Without it the cache would keep the old bytes -- refcache reuses any
            # existing destination -- while the provenance below copied the *new* export
            # date out of the index, claiming a freshness the file did not have.
            dest = reference_cache / f"openfda.{dataset}.{export}.{number:04d}.json.zip"
            openfda_paths.append(dest)
            provenance.append({"source": "openFDA", "dataset": dataset,
                               "export_date": entry.get("export_date"),
                               "total_records": entry.get("total_records"),
                               "licence": OPENFDA_LICENCE, "redistributable": True,
                               **refcache.fetch(part["file"], dest)})
    return {"unichem": unichem_dest, "openfda": openfda_paths}, provenance


def _download_index() -> dict:
    import urllib.request

    request = urllib.request.Request(OPENFDA_DOWNLOADS, headers=refcache.USER_AGENT)
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def merge_on_rxcui(identities: dict) -> tuple:
    """Second pass: fold together identities RxNorm considers the same product.

    Open Targets' ``parentId`` is authoritative where it is populated, and it is not
    always populated. On the real data ixazomib and ixazomib citrate arrived as two
    identities holding *identical* RxCUI sets, and took ranks 2 and 3 — one molecule
    occupying two places in the ranking, which is the exact failure harmonisation exists
    to prevent.

    The rule is deliberately narrow: merge only when the whole RxCUI tuple is identical and
    non-empty. Two different drugs do not share a complete RxCUI set, while a salt and its
    parent routinely do. A subset relation would be looser and is not used, because a
    combination product legitimately contains its components' codes.

    The surviving name is the shortest among the merged members — a salt's name is its
    parent's plus a counterion — with the smallest ChEMBL id breaking ties so the choice
    never depends on iteration order. Returns ``(identities, n_merged)``.
    """
    by_rxcui: dict = {}
    for identity in {i.key: i for i in identities.values()}.values():
        if identity.rxcui:
            by_rxcui.setdefault(identity.rxcui, []).append(identity)

    remap, merged = {}, 0
    for group in by_rxcui.values():
        if len(group) < 2:
            continue
        winner = min(group, key=lambda i: (len(i.name or "~" * 99), i.key))
        combined = frozenset().union(*(i.chembl_ids for i in group))
        survivor = Identity(
            key=winner.key, name=winner.name, chembl_ids=combined,
            rxcui=winner.rxcui, unii=tuple(sorted({u for i in group for u in i.unii})),
            pharm_class_epc=tuple(sorted({c for i in group for c in i.pharm_class_epc})),
            pharm_class_moa=tuple(sorted({c for i in group for c in i.pharm_class_moa})),
            resolution=winner.resolution, drug_type=winner.drug_type,
            clinical_stage=winner.clinical_stage)
        merged += len(group) - 1
        for member in group:
            remap[member.key] = survivor

    if remap:
        identities = {cid: remap.get(identity.key, identity)
                      for cid, identity in identities.items()}
    return identities, merged


def build(config: dict, chembl_ids) -> tuple:
    """Harmonise a set of ChEMBL ids. Returns ``({chembl_id: Identity}, stats, provenance)``.

    Every input id gets an entry, including one whose RxCUI could not be resolved -- see
    the guardrail. Ids that share a parent share one :class:`Identity` object.
    """
    paths, provenance = ensure_sources(config)
    molecule_paths, molecule_provenance = enrichment.ensure_datasets(
        config, ("drug_molecule",), channel="channel_b")
    provenance.extend(molecule_provenance)

    parent, identity_rows, parent_stats = parent_index(molecule_paths["drug_molecule"])
    chembl_to_unii = load_chembl_to_unii(paths["unichem"])
    by_unii, by_name, openfda_stats = load_openfda(paths["openfda"])

    # Group the requested ids by the parent they collapse onto.
    groups: dict = {}
    for chembl_id in sorted(set(chembl_ids)):
        groups.setdefault(resolve_parent(chembl_id, parent), set()).add(chembl_id)

    identities, stats = {}, {r: 0 for r in RESOLUTION_ORDER}
    stats["collapsed_ids"] = 0
    for key, members in sorted(groups.items()):
        row = identity_rows.get(key, {})
        name = row.get("name") or next(
            (identity_rows.get(m, {}).get("name", "") for m in sorted(members)), "")

        # A parent's own UNII first, then any member's: the parent is the substance the
        # registry knows, and a salt form's UNII is a different registered substance.
        found, resolution = {}, "unresolved"
        uniis = sorted(chembl_to_unii.get(key, set()))
        for member in sorted(members):
            uniis += [u for u in sorted(chembl_to_unii.get(member, set())) if u not in uniis]
        for unii in uniis:
            if unii in by_unii:
                found, resolution = by_unii[unii], "unii"
                break
        if not found:
            hit = by_name.get(normalise(name))
            if hit:
                found, resolution = hit, "name"

        stats[resolution] += 1
        stats["collapsed_ids"] += max(0, len(members) - 1)
        identity = Identity(
            key=key, name=name, chembl_ids=frozenset(members),
            rxcui=tuple(sorted(found.get("rxcui", ()))),
            unii=tuple(sorted(found.get("unii", ())) or uniis[:1]),
            pharm_class_epc=tuple(sorted(found.get("pharm_class_epc", ()))),
            pharm_class_moa=tuple(sorted(found.get("pharm_class_moa", ()))),
            resolution=resolution,
            drug_type=row.get("drug_type", ""), clinical_stage=row.get("clinical_stage", ""))
        for member in members:
            identities[member] = identity

    identities, merged = merge_on_rxcui(identities)
    stats["merged_on_identical_rxcui"] = merged
    stats.update({"molecules_in": len(set(chembl_ids)),
                  "identities_out": len({i.key for i in identities.values()}),
                  **{f"openfda_{k}": v for k, v in openfda_stats.items()},
                  "chembl_ids_with_a_unii": len(chembl_to_unii), **parent_stats})
    # The identity count is the one after the RxCUI merge, not len(groups) before it --
    # logging the pre-merge number would disagree with the table the run then writes.
    LOGGER.info("harmonised %d ChEMBL id(s) into %d identity(ies) (%d collapsed onto a "
                "parent, %d merged on an identical RxCUI set); rxcui: %d by UNII, %d by "
                "name, %d unresolved", len(set(chembl_ids)), stats["identities_out"],
                stats["collapsed_ids"], stats["merged_on_identical_rxcui"],
                stats["unii"], stats["name"], stats["unresolved"])
    return identities, stats, provenance
