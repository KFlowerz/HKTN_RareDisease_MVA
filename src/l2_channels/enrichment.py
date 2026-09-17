"""L2 -- the non-redistributed enrichment zone: drug-target data, kept behind the counter.

Purpose
    Supply Channel B with each drug's protein targets, from a source whose licence does
    not permit redistribution, without any of that source's content reaching a published
    output.

Inputs
    ``config["l2"]["channel_b"]`` (``open_targets_release``, ``approved_only``) and
    ``config["enrichment_dir"]``; a map from Ensembl gene id to the interactome's protein
    identifiers, built from STRING's own alias table (CC BY 4.0).

Outputs
    :class:`DrugRecord` per drug -- identifiers, name, clinical stage, and the target set
    in the interactome's identifier space -- plus a provenance record per file.
    :func:`publishable` reduces a record to the fields that may leave this zone.

Guardrail
    **This is the segregated zone described in ``../l4_validate/sources.md``.** Open
    Targets marks its Platform data CC0, but its drug-target content derives from ChEMBL,
    which is CC BY-SA 3.0 -- and ShareAlike is viral: bundling that content into a derived
    table would relicense every output the project publishes under CC BY 4.0
    (``COMPLIANCE.md``). This project does not rely on one aggregator's relicensing of
    another's ShareAlike data; it takes the conservative reading and keeps the content
    here.

    **What may leave this module**: identifiers (ChEMBL id), the drug's name, its clinical
    stage, and counts -- the fields :data:`PUBLISHABLE_FIELDS` names, which Table 2 of the
    sources ledger already permits for performing a join. Names are public INN-style
    identifiers, not curated annotation.

    **What may not**: the target sets themselves, the mechanism-of-action text, and the
    action types. Those are the curated content. :func:`publishable` is the only sanctioned
    way out, and ``tests/test_l2_channel_b.py`` asserts that what Channel B writes carries
    no target identifier.

    Cached files live in ``enrichment_dir``, outside the repository and outside the custody
    root. They are not patient data and must never enter the deletion register in
    ``docs/data_custody.md``; they are also never committed, and ``src/refcache.py``
    refuses a location inside the working tree.

    **No patient data reaches this module.** It never opens ``data_dir``, and every
    download is a whole public release with no query attached.
"""

from __future__ import annotations

import ast
import gzip
import os
import re
import urllib.request
from dataclasses import dataclass

from .. import refcache

#: Separate from the public reference cache on purpose: the separation between what may
#: be redistributed and what may not is a directory boundary, not a convention.
DEFAULT_DIR = "~/.cache/mva-track2/enrichment"
OPEN_TARGETS_FTP = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
#: The two datasets Channel B needs: drug -> target, and drug -> identity/stage.
DATASETS = ("drug_mechanism_of_action", "drug_molecule")
#: Open Targets' own licence mark; the restriction below it is ChEMBL's, and is why this
#: zone exists. Both travel in the provenance record.
LICENCE = "Open Targets: CC0 1.0; drug-target content derived from ChEMBL (CC BY-SA 3.0)"
REDISTRIBUTABLE = False
#: The only fields that may leave this zone. See the guardrail.
PUBLISHABLE_FIELDS = ("chembl_id", "name", "clinical_stage", "drug_type", "n_targets")
#: Open Targets' value for "reached marketing approval".
APPROVED = "APPROVAL"


@dataclass(frozen=True)
class DrugRecord:
    """One drug, with its targets in the interactome's identifier space."""

    chembl_id: str
    name: str = ""
    drug_type: str = ""
    clinical_stage: str = ""
    targets: frozenset = frozenset()          # restricted -- never published
    unmapped_targets: frozenset = frozenset()  # Ensembl ids absent from the interactome

    @property
    def n_targets(self) -> int:
        return len(self.targets)

    @property
    def approved(self) -> bool:
        return self.clinical_stage == APPROVED


def enrichment_dir(config: dict):
    """Where restricted-licence files are cached. ``$MVA_ENRICHMENT_ROOT`` overrides config.

    Raises:
        ValueError: If it resolves to the same directory as the public reference cache.
            The whole point of decision D7 is that redistributable and non-redistributable
            data do not share a location; collapsing the two by configuration would leave
            no boundary at all, and nothing downstream would notice.
    """
    path = refcache.cache_dir(
        os.environ.get("MVA_ENRICHMENT_ROOT") or config.get("enrichment_dir"),
        DEFAULT_DIR, name="enrichment_dir")
    if path == refcache.reference_dir(config):
        raise ValueError(
            f"enrichment_dir and reference_dir both resolve to {path}. Licence-restricted "
            "data must not share a directory with redistributable data (decision D7)."
        )
    return path


def publishable(record: DrugRecord) -> dict:
    """The subset of a record that may appear in a redistributed output.

    Deliberately a whitelist, not a blacklist: a field added to :class:`DrugRecord` later
    stays inside this zone until someone decides otherwise.
    """
    values = {"chembl_id": record.chembl_id, "name": record.name,
              "clinical_stage": record.clinical_stage, "drug_type": record.drug_type,
              "n_targets": record.n_targets}
    return {field: values[field] for field in PUBLISHABLE_FIELDS}


def _get(url: str, timeout: int = 300) -> bytes:
    request = urllib.request.Request(url, headers=refcache.USER_AGENT)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def part_urls(base: str) -> list:
    """Every Parquet part in a released dataset directory, in listed order."""
    html = _get(base, timeout=120).decode("utf-8", "replace")
    names = sorted({h for h in re.findall(r'href="([^"?][^"]*)"', html)
                    if h.endswith(".parquet")})
    return [base + name for name in names]


def ensure_datasets(config: dict) -> tuple:
    """Download (or reuse) each dataset's parts. Returns ``(paths, provenance)``."""
    settings = ((config.get("l2") or {}).get("channel_b") or {})
    release = str(settings.get("open_targets_release", "26.06"))
    cache = enrichment_dir(config)

    paths, provenance = {}, []
    for dataset in DATASETS:
        base = f"{OPEN_TARGETS_FTP}/{release}/output/{dataset}/"
        urls = part_urls(base)
        if not urls:
            # An empty listing means the release, the directory layout or the link markup
            # changed. Continuing would produce zero rows, which reads downstream as
            # "no drug has a target" rather than "the download found nothing".
            raise ValueError(
                f"no Parquet part found at {base} -- check that release {release!r} exists "
                "and still publishes this dataset"
            )
        paths[dataset] = []
        for url in urls:
            dest = cache / f"{release}.{dataset}.{url.rsplit('/', 1)[-1]}"
            paths[dataset].append(dest)
            provenance.append({"source": "Open Targets Platform", "release": release,
                               "dataset": dataset, "licence": LICENCE,
                               "redistributable": REDISTRIBUTABLE,
                               **refcache.fetch(url, dest)})
    return paths, provenance


def _rows(path):
    import pyarrow.parquet as pq

    return pq.read_table(path).to_pylist()


def _as_list(value) -> list:
    """Open Targets writes list columns as lists, but a cached round-trip can stringify."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return [value]
        return list(parsed) if isinstance(parsed, (list, tuple)) else [parsed]
    return [value]


def load_mechanisms(paths, ensembl_to_protein: dict) -> tuple:
    """``({chembl_id: (targets, unmapped)}, n_rows)`` from the mechanism dataset.

    A mechanism row names one target and the drugs acting on it, so the mapping is built
    by accumulating each row's target onto every drug it lists.
    """
    targets: dict = {}
    unmapped: dict = {}
    seen = 0
    for path in paths:
        for row in _rows(path):
            seen += 1
            genes = _as_list(row.get("targets"))
            drugs = _as_list(row.get("chemblIds"))
            for drug in drugs:
                targets.setdefault(drug, set())
                for gene in genes:
                    protein = ensembl_to_protein.get(gene)
                    if protein:
                        targets[drug].add(protein)
                    else:
                        unmapped.setdefault(drug, set()).add(gene)
    # A drug whose genes all failed to map keeps an entry with an empty target set. It is
    # dropped later, by drug_records -- but only after it has been counted, because
    # "every drug lost its targets" is exactly how a broken crosswalk looks, and a
    # statistic that can never fire cannot report it.
    return {drug: (frozenset(ts), frozenset(unmapped.get(drug, ())))
            for drug, ts in targets.items()}, seen


def load_molecules(paths) -> dict:
    """``{chembl_id: {name, drug_type, clinical_stage}}`` from the molecule dataset."""
    out = {}
    for path in paths:
        for row in _rows(path):
            chembl_id = row.get("id")
            if chembl_id:
                out[chembl_id] = {
                    "name": row.get("name") or "",
                    "drug_type": row.get("drugType") or "",
                    "clinical_stage": row.get("maximumClinicalStage") or "",
                }
    return out


def ensembl_index(alias_path, source: str = "Ensembl_gene") -> dict:
    """``{Ensembl gene id: protein}`` from STRING's alias table (CC BY 4.0).

    The crosswalk comes from the interactome's own aliases, not from the restricted
    source, so nothing here depends on the enrichment zone for identity.
    """
    opener = gzip.open if str(alias_path).endswith(".gz") else open
    out = {}
    with opener(alias_path, "rt", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] == source:
                out.setdefault(fields[1], fields[0])
    return out


def drug_records(config: dict, ensembl_to_protein: dict) -> tuple:
    """Every drug with at least one mapped target. Returns ``(records, provenance, stats)``."""
    settings = ((config.get("l2") or {}).get("channel_b") or {})
    approved_only = bool(settings.get("approved_only", True))

    paths, provenance = ensure_datasets(config)
    mechanisms, n_rows = load_mechanisms(paths["drug_mechanism_of_action"], ensembl_to_protein)
    molecules = load_molecules(paths["drug_molecule"])

    records = []
    for chembl_id, (targets, unmapped) in sorted(mechanisms.items()):
        identity = molecules.get(chembl_id, {})
        record = DrugRecord(chembl_id=chembl_id, targets=targets, unmapped_targets=unmapped,
                            **{k: identity.get(k, "") for k in
                               ("name", "drug_type", "clinical_stage")})
        if targets and (record.approved or not approved_only):
            records.append(record)

    stats = {
        "mechanism_rows": n_rows,
        "drugs_with_a_mechanism": len(mechanisms),
        "molecules": len(molecules),
        "drugs_kept": len(records),
        "approved_only": approved_only,
        "drugs_with_unmappable_targets_only": sum(
            1 for _, (t, u) in mechanisms.items() if not t and u),
    }
    return records, provenance, stats
