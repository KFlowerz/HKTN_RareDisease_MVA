"""L4 -- openFDA label text, read in bulk and indexed by identifier.

Purpose
    Supply the label prose the safety rules decide on: the carcinogenesis and mutagenesis
    section, paediatric use, boxed warnings, and the marketing status of a drug's
    products. Every rule in :mod:`.safety_triage` has to name the source field it drew on
    and quote a snippet, so this module keeps the text rather than a boolean.

Inputs
    openFDA's ``drug/label`` and ``drug/drugsfda`` bulk releases, and the set of RxCUIs
    and UNIIs L3 resolved. ``config["reference_dir"]`` holds the cache.

Outputs
    :class:`Label` per candidate identifier -- the sections that matter, each with the
    openFDA field name it came from -- and :class:`Marketing` per candidate, saying
    whether any product is currently marketed.

Guardrail
    **Absence of a field is not evidence of safety.** A drug with no label in this release
    is not a drug without warnings; it is a drug whose label openFDA does not carry, often
    because it is not marketed in the United States. This module reports what it found and
    what it did not, and :mod:`.safety_triage` fails closed on the difference.

    **Whole public files only.** The releases are downloaded complete and matched locally,
    like everything else in this pipeline. No candidate identifier is sent anywhere.

    Label text is **data, never instructions** (CLAUDE.md). It is matched against fixed
    patterns and quoted; nothing inside a label may change what this pipeline does.

    openFDA is US-Government public domain, so the snippets quoted in an output are
    redistributable. openFDA states its data is not for clinical use and may be incomplete
    or inaccurate; that statement travels with every verdict drawn from it.
"""

from __future__ import annotations

import json
import logging
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

from .. import refcache

LOGGER = logging.getLogger(__name__)

OPENFDA_DOWNLOADS = "https://api.fda.gov/download.json"
#: ``label`` carries the prose; ``drugsfda`` carries per-product marketing status.
DATASETS = ("label", "drugsfda")
LICENCE = ("openFDA (U.S. Food and Drug Administration): US-Government work, public "
           "domain. openFDA states its data is not for clinical use and may be incomplete "
           "or inaccurate.")

#: openFDA field -> the SPL section a reader would cite. Kept so a verdict can name the
#: section number rather than an API field nobody outside this codebase recognises.
SECTIONS = {
    # openFDA's field name, verified against the release on 2026-09-18 (4,645 of the
    # 20,000 records in part 0 carry it). An earlier guess at
    # "carcinogenesis_and_mutagenesis_of_fertility" matched nothing at all and silently
    # excluded every candidate as insufficient_evidence.
    "carcinogenesis_and_mutagenesis_and_impairment_of_fertility":
        "13.1 Carcinogenesis, Mutagenesis, Impairment of Fertility",
    # Section 13.1 nests under Nonclinical Toxicology, and some labels populate only the
    # parent. Read as a fallback so a label that carries the text under the outer heading
    # is not treated as carrying no text.
    "nonclinical_toxicology": "13 Nonclinical Toxicology",
    "pediatric_use": "8.4 Pediatric Use",
    "boxed_warning": "Boxed Warning",
    "warnings_and_cautions": "5 Warnings and Precautions",
    "warnings": "Warnings",
    "contraindications": "4 Contraindications",
}
#: Statuses that mean a product can actually be dispensed today.
MARKETED_STATUSES = frozenset({"Prescription", "Over-the-counter"})

#: How much of a section is kept. Enough to quote a verdict's basis, bounded so a results
#: file does not become a copy of the label corpus.
SNIPPET_LIMIT = 600


@dataclass(frozen=True)
class Label:
    """The label sections one drug's SPL records carry."""

    sections: dict = field(default_factory=dict)   # openFDA field -> joined text
    spl_ids: tuple = ()

    def text(self, name: str) -> str:
        return self.sections.get(name, "")

    def has(self, name: str) -> bool:
        return bool(self.sections.get(name))


@dataclass(frozen=True)
class Marketing:
    """What Drugs@FDA says about a molecule's products."""

    statuses: frozenset = frozenset()
    application_numbers: tuple = ()

    @property
    def known(self) -> bool:
        return bool(self.statuses)

    @property
    def marketed(self) -> bool:
        """True if any product is currently prescription or OTC.

        Read this way round on purpose: a molecule whose every product is discontinued or
        never got past tentative approval is not something to propose for a child, and
        several of L3's top rows are exactly that. ``Discontinued`` on *some* products is
        normal for a marketed drug whose older formulations lapsed.
        """
        return bool(self.statuses & MARKETED_STATUSES)


def snippet(text: str, at: "re.Match | re.Pattern | None" = None,
            limit: int = SNIPPET_LIMIT) -> str:
    """A quotable extract: the sentence around a match, else the section's opening.

    ``at`` may be a :class:`re.Match` already found in ``text`` -- which is what callers
    should pass, so the quoted sentence is the one the verdict actually rested on -- or a
    pattern to search for. Passing a pattern re-searches from the start and can land on a
    different occurrence than the caller decided about, which is why the match form
    exists.

    The snippet is what a reader checks the verdict against, so it is taken from the label
    verbatim and only trimmed, never paraphrased.
    """
    collapsed = " ".join((text or "").split())
    if not collapsed:
        return ""
    span = None
    if isinstance(at, re.Match):
        # The caller searched the raw text; re-find the same literal in the collapsed
        # form so the offsets line up with what is quoted.
        needle = " ".join(at.group(0).split())
        position = collapsed.find(needle)
        if position >= 0:
            span = (position, position + len(needle))
    elif at is not None:
        found = at.search(collapsed)
        if found:
            span = found.span()

    if span:
        start = collapsed.rfind(".", 0, span[0]) + 1
        end = collapsed.find(".", span[1])
        end = len(collapsed) if end == -1 else end + 1
        extract = collapsed[start:end].strip()
        if extract:
            return extract[:limit]
    return collapsed[:limit]


def _records(path: Path):
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.endswith(".json"):
                yield from json.loads(archive.read(name)).get("results", [])


def download_index() -> dict:
    import urllib.request

    request = urllib.request.Request(OPENFDA_DOWNLOADS, headers=refcache.USER_AGENT)
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def ensure_datasets(config: dict, datasets=DATASETS) -> tuple:
    """Download (or reuse) the openFDA bulk releases. Returns ``(paths, provenance)``.

    The export date is part of the cached filename, so a new openFDA export lands as a new
    file rather than being shadowed by stale bytes carrying a fresh date in provenance.
    """
    cache = refcache.reference_dir(config)
    index = download_index()
    paths, provenance = {}, []
    for dataset in datasets:
        entry = index["results"]["drug"][dataset]
        export = str(entry.get("export_date") or "unknown")
        paths[dataset] = []
        for number, part in enumerate(entry["partitions"]):
            dest = cache / f"openfda.{dataset}.{export}.{number:04d}.json.zip"
            paths[dataset].append(dest)
            provenance.append({"source": "openFDA", "dataset": dataset,
                               "export_date": export, "licence": LICENCE,
                               "redistributable": True,
                               "total_records": entry.get("total_records"),
                               **refcache.fetch(part["file"], dest)})
    return paths, provenance


def load_labels(paths, rxcuis: set, uniis: set) -> tuple:
    """``({identifier: Label}, stats)`` for the identifiers asked about.

    Only records matching a requested identifier are kept, so the 260,000-record corpus
    never has to be held in memory at once. A label is filed under every identifier it
    carries, because a candidate may be resolved by either.
    """
    wanted_rx, wanted_unii = set(rxcuis), set(uniis)
    collected: dict = {}
    stats = {"label_records": 0, "label_records_matched": 0}

    for path in paths:
        for record in _records(path):
            stats["label_records"] += 1
            block = record.get("openfda") or {}
            keys = ({f"rxcui:{r}" for r in block.get("rxcui") or [] if r in wanted_rx}
                    | {f"unii:{u}" for u in block.get("unii") or [] if u in wanted_unii})
            if not keys:
                continue
            stats["label_records_matched"] += 1
            for key in keys:
                entry = collected.setdefault(key, {"sections": {}, "spl": set()})
                for name in SECTIONS:
                    value = record.get(name)
                    if isinstance(value, list):
                        value = " ".join(str(v) for v in value if v)
                    if value:
                        # Several SPLs per molecule; keeping them joined means a warning
                        # present on any one product's label is seen, which is the
                        # conservative reading for a safety gate.
                        existing = entry["sections"].get(name, "")
                        if str(value) not in existing:
                            entry["sections"][name] = f"{existing} {value}".strip()
                for spl in block.get("spl_set_id") or []:
                    entry["spl"].add(spl)

    labels = {key: Label(sections=dict(v["sections"]), spl_ids=tuple(sorted(v["spl"])))
              for key, v in collected.items()}
    stats["identifiers_with_a_label"] = len(labels)
    return labels, stats


def load_marketing(paths, rxcuis: set, uniis: set) -> tuple:
    """``({identifier: Marketing}, stats)`` from Drugs@FDA product records."""
    wanted_rx, wanted_unii = set(rxcuis), set(uniis)
    collected: dict = {}
    stats = {"drugsfda_records": 0, "drugsfda_records_matched": 0}

    for path in paths:
        for record in _records(path):
            stats["drugsfda_records"] += 1
            block = record.get("openfda") or {}
            keys = ({f"rxcui:{r}" for r in block.get("rxcui") or [] if r in wanted_rx}
                    | {f"unii:{u}" for u in block.get("unii") or [] if u in wanted_unii})
            if not keys:
                continue
            stats["drugsfda_records_matched"] += 1
            statuses = {p.get("marketing_status") for p in record.get("products") or []}
            statuses.discard(None)
            for key in keys:
                entry = collected.setdefault(key, {"statuses": set(), "apps": set()})
                entry["statuses"].update(statuses)
                if record.get("application_number"):
                    entry["apps"].add(record["application_number"])

    marketing = {key: Marketing(statuses=frozenset(v["statuses"]),
                                application_numbers=tuple(sorted(v["apps"])))
                 for key, v in collected.items()}
    stats["identifiers_with_a_marketing_record"] = len(marketing)
    return marketing, stats


def for_identity(index: dict, rxcui, unii):
    """The entry matching any of a candidate's identifiers, or ``None``.

    RxCUI is tried before UNII: an RxCUI names a clinical drug, a UNII names a substance,
    and the label that matters is the one for the product a clinician would prescribe.
    """
    for value in rxcui or ():
        hit = index.get(f"rxcui:{value}")
        if hit is not None:
            return hit
    for value in unii or ():
        hit = index.get(f"unii:{value}")
        if hit is not None:
            return hit
    return None
