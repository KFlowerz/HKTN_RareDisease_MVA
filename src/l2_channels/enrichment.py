"""L2 -- the non-redistributed enrichment zone: drug-target data, kept behind the counter.

Purpose
    Supply Channel B with each drug's protein targets, Channel D with each drug's approved
    indications, and Channel E with a name-to-identifier crosswalk, from a source whose
    licence does not permit redistribution, without any of that source's content reaching a
    published output.

Inputs
    ``config["l2"]["channel_b"]`` / ``["channel_d"]`` / ``["channel_e"]``
    (``open_targets_release``, ``approved_only``) and ``config["enrichment_dir"]``; for
    Channel B, a map from Ensembl gene id to the interactome's protein identifiers, built
    from STRING's own alias table (CC BY 4.0).

Outputs
    :class:`DrugRecord` per drug -- identifiers, name, clinical stage, and the target set
    in the interactome's identifier space -- plus a provenance record per file;
    :func:`load_indications` -> disease id -> approved drugs, which likewise never leaves
    the zone except as counts. :func:`publishable` reduces a record to the fields that may
    leave this zone.

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


def ensure_datasets(config: dict, datasets=DATASETS, *, channel: str = "channel_b") -> tuple:
    """Download (or reuse) each dataset's parts. Returns ``(paths, provenance)``.

    ``channel`` names the ``config["l2"]`` block whose ``open_targets_release`` applies, so
    each channel records the release it actually used.
    """
    settings = ((config.get("l2") or {}).get(channel) or {})
    release = str(settings.get("open_targets_release", "26.06"))
    cache = enrichment_dir(config)

    paths, provenance = {}, []
    for dataset in datasets:
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


def name_index(paths) -> tuple:
    """``({normalised name: chembl_id}, stats)`` -- resolve a compound name to an identity.

    Channel E starts from compound *names*, which is how the literature refers to a drug,
    and L3 aggregates on identifiers, so a name has to become a ChEMBL id somewhere. Doing
    it here keeps that identity the same one Channels B and D used.

    Preferred names win over synonyms and trade names, so ``name`` never loses to another
    molecule's synonym; among synonyms a collision keeps the first by ChEMBL id, which
    makes the result order-independent. Collisions are counted, because a name that maps to
    two molecules is a resolution this channel should not be trusted to have made.
    """
    preferred: dict = {}
    secondary: dict = {}
    collisions: set = set()
    for path in paths:
        for row in _rows(path):
            chembl_id = row.get("id")
            if not chembl_id:
                continue
            name = normalise_name(row.get("name") or "")
            if name:
                if name in preferred and preferred[name] != chembl_id:
                    collisions.add(name)
                preferred.setdefault(name, chembl_id)
            for other in _as_list(row.get("synonyms")) + _as_list(row.get("tradeNames")):
                key = normalise_name(str(other))
                if not key:
                    continue
                if key in secondary and secondary[key] != chembl_id:
                    collisions.add(key)
                    secondary[key] = min(secondary[key], chembl_id)
                else:
                    secondary.setdefault(key, chembl_id)
    index = {**secondary, **preferred}
    # A preferred name settles its own key, so only a collision among synonyms is left
    # ambiguous. The names themselves are returned, not just a count: a caller resolving a
    # handful of compounds needs to know whether *its* names are the ambiguous ones.
    ambiguous = frozenset(collisions - set(preferred))
    return index, {"molecule_names": len(preferred), "molecule_synonyms": len(secondary),
                   "ambiguous_synonyms": ambiguous}


def normalise_name(value: str) -> str:
    """Casefold and collapse punctuation, so ``17-AAG`` and ``17 AAG`` are one key."""
    return re.sub(r"[^a-z0-9]+", " ", value.strip().lower()).strip()


def load_indications(paths, *, stage: str = APPROVED) -> tuple:
    """``({disease id: {chembl_id}}, stats)`` from Open Targets' ``clinical_indication``.

    Only indications whose highest clinical stage *for that indication* equals ``stage``
    are kept: a drug approved for one disease and in a phase-2 trial for another is not an
    approved treatment of the second. Disease ids are Open Targets' (``MONDO_0012941``,
    ``HP_0001250``, ``EFO_...``). The indication pairs are ChEMBL-derived content and never
    leave this zone except as counts.
    """
    out: dict = {}
    stats = {"indication_rows": 0, "indications_at_stage": 0, "stage": stage}
    for path in paths:
        for row in _rows(path):
            stats["indication_rows"] += 1
            if row.get("maxClinicalStage") != stage:
                continue
            drug, disease = row.get("drugId"), row.get("diseaseId")
            if drug and disease:
                out.setdefault(disease, set()).add(drug)
                stats["indications_at_stage"] += 1
    return out, stats


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
