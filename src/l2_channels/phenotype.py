"""L2 channel D -- phenotype inputs: the subject's HPO terms and the public phenotype corpus.

Purpose
    Read the subject's HPO-coded phenotype from ``data_dir``, and fetch and parse the two
    public resources channel D matches it against: the Human Phenotype Ontology and
    Monarch's disease-to-phenotype associations.

Inputs
    ``config["data_dir"]`` (the phenotype ``.docx``), ``config["reference_dir"]``, and
    ``config["l2"]["channel_d"]`` for release URLs and the document pattern.

Outputs
    :func:`read_terms` -> the subject's HPO ids; :func:`parse_obo` -> the ontology's
    parents, names, and obsolete-term replacements; :func:`load_associations` -> disease
    -> annotated HPO ids and disease labels; :func:`ensure_sources` -> cached file paths
    and provenance.

Guardrail
    **The phenotype is patient data** (``COMPLIANCE.md``). :func:`read_terms` returns ids
    and nothing else; no function here logs, prints or writes a term, and error messages
    name row numbers, never content. A combination of features is identifying in a
    population of roughly 50 patients.

    **Absence of a feature is never guessed.** The document records features as present.
    A row whose cells carry a negation word ("no", "absent", "excluded", ...) is refused,
    not parsed, because reading "absent X" as X would invert the query; the fix is a
    person reading that row.

    **Nothing about the subject leaves the machine.** The two public resources are
    downloaded whole and matched locally -- never queried with the subject's terms, which
    is exactly what Monarch's web API would require (decision D10).

    Licences (``../l4_validate/sources.md``, Table 6): the HPO files are used unaltered, and
    its version is carried into every output, as the HPO licence requires; Monarch's
    association table labels diseases with Mondo (CC BY 4.0) names, so no OMIM text is
    used.
"""

from __future__ import annotations

import gzip
import re
import zipfile
from pathlib import Path

from .. import refcache

MONARCH_URL = ("https://data.monarchinitiative.org/monarch-kg/2026-09-02/tsv/"
               "disease_associations/disease_phenotype.all.tsv.gz")
HPO_URL = "http://purl.obolibrary.org/obo/hp/releases/2026-09-01/hp.obo"
DOCUMENT_GLOB = "*Phenotype*.docx"
#: Root of the HPO's phenotypic-abnormality subtree. Terms outside it -- mode of
#: inheritance, onset, frequency -- describe how a feature presents, not a feature.
PHENOTYPIC_ABNORMALITY = "HP:0000118"

_HP_ID = re.compile(r"HP:\d{7}")
_NEGATION = re.compile(r"\b(no|not|absent|absence|without|negative|excluded|denies|denied|"
                       r"ruled out|none)\b", re.IGNORECASE)
_W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


# ------------------------------------------------------------------ the document


def find_document(config: dict) -> Path:
    """The single phenotype document under ``data_dir``.

    Raises:
        FileNotFoundError / ValueError: On none, or on more than one -- choosing between
            two phenotype descriptions is a person's decision.
    """
    pattern = ((config.get("l2") or {}).get("channel_d") or {}).get("document_glob",
                                                                     DOCUMENT_GLOB)
    matches = sorted(Path(config["data_dir"]).rglob(pattern))
    if not matches:
        raise FileNotFoundError(f"no phenotype document matching {pattern!r} under data_dir")
    if len(matches) > 1:
        raise ValueError(f"{len(matches)} documents match {pattern!r} under data_dir; "
                         "set l2.channel_d.document_glob to exactly one")
    return matches[0]


def _cell_text(cell) -> str:
    return "".join(t.text or "" for t in cell.iter(f"{_W}t"))


def read_terms(path: Path) -> list:
    """The HPO ids a ``.docx`` records as present, deduplicated, in document order.

    Reads every table. A row contributes the HPO ids in its cells; a row with a negation
    word in any cell is refused. HPO ids outside tables are refused too: their status
    (present? historical? a differential?) cannot be read from structure.

    Raises:
        ValueError: On a negated row, on HPO ids outside a table, or on no id at all.
            Messages carry row numbers only.
    """
    from xml.etree import ElementTree

    with zipfile.ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))

    terms, table_ids = [], 0
    for t_index, table in enumerate(root.iter(f"{_W}tbl"), start=1):
        for r_index, row in enumerate(table.iter(f"{_W}tr"), start=1):
            cells = [_cell_text(c) for c in row.iter(f"{_W}tc")]
            ids = [m for c in cells for m in _HP_ID.findall(c)]
            if not ids:
                continue
            table_ids += len(ids)
            if any(_NEGATION.search(c) for c in cells):
                raise ValueError(
                    f"table {t_index}, row {r_index} carries an HPO id and a negation word. "
                    "Refusing to guess whether the feature is present; review that row.")
            terms.extend(ids)

    body_ids = len(_HP_ID.findall("".join(t.text or "" for t in root.iter(f"{_W}t"))))
    if body_ids != table_ids:
        raise ValueError(f"{body_ids - table_ids} HPO id(s) appear outside a table, where "
                         "their status cannot be read from structure; review the document.")
    if not terms:
        raise ValueError("the phenotype document holds no HPO id")
    return list(dict.fromkeys(terms))


# -------------------------------------------------------------------- the HPO


def parse_obo(lines) -> dict:
    """Parse an OBO file.

    Returns ``{"version", "parents", "names", "replaced_by", "alt_ids", "obsolete"}``:
    live terms' ``is_a`` parents and names; for obsolete terms their ``replaced_by``;
    ``alt_id`` -> primary id.
    """
    version = ""
    parents, names, replaced, alt_ids, obsolete = {}, {}, {}, {}, set()
    block, current = None, {}

    def flush():
        if block != "Term" or "id" not in current:
            return
        term = current["id"][0]
        for alt in current.get("alt_id", ()):
            alt_ids[alt] = term
        if current.get("is_obsolete", ["false"])[0] == "true":
            obsolete.add(term)
            if current.get("replaced_by"):
                replaced[term] = current["replaced_by"][0]
            return
        parents[term] = [p.split("!")[0].strip() for p in current.get("is_a", ())]
        names[term] = current.get("name", [""])[0]

    for raw in lines:
        line = raw.rstrip("\n")
        if line.startswith("["):
            flush()
            block, current = line.strip("[]"), {}
            continue
        if block is None and line.startswith("data-version:"):
            version = line.split(":", 1)[1].strip()
        if block and ":" in line:
            key, value = line.split(":", 1)
            current.setdefault(key.strip(), []).append(value.strip())
    flush()
    return {"version": version, "parents": parents, "names": names,
            "replaced_by": replaced, "alt_ids": alt_ids, "obsolete": obsolete}


def resolve(term: str, obo: dict):
    """A live term id for ``term`` -- itself, its primary id, or its replacement -- or None."""
    seen = set()
    while term not in obo["parents"]:
        if term in seen:
            return None
        seen.add(term)
        if term in obo["alt_ids"]:
            term = obo["alt_ids"][term]
        elif term in obo["replaced_by"]:
            term = obo["replaced_by"][term]
        else:
            return None
    return term


# ------------------------------------------------------- Monarch's associations


def load_associations(lines) -> tuple:
    """Monarch ``disease_phenotype`` rows -> ``({disease: {hpo ids}}, {disease: label}, stats)``.

    Negated associations are dropped and counted: "disease X lacks feature Y" is not an
    annotation of Y. Only Mondo-identified diseases are kept, so labels are Mondo's.
    """
    rows = iter(lines)
    header = next(rows).rstrip("\n").split("\t")
    col = {name: i for i, name in enumerate(header)}
    for needed in ("subject", "subject_label", "negated", "object"):
        if needed not in col:
            raise ValueError(f"Monarch association table has no {needed!r} column -- the "
                             "release format changed")
    annotations, labels = {}, {}
    stats = {"rows": 0, "negated": 0, "non_mondo_subject": 0, "non_hpo_object": 0}
    for line in rows:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < len(header):
            continue
        stats["rows"] += 1
        subject, obj = fields[col["subject"]], fields[col["object"]]
        if fields[col["negated"]].strip().lower() == "true":
            stats["negated"] += 1
            continue
        if not subject.startswith("MONDO:"):
            stats["non_mondo_subject"] += 1
            continue
        if not obj.startswith("HP:"):
            stats["non_hpo_object"] += 1
            continue
        annotations.setdefault(subject, set()).add(obj)
        labels.setdefault(subject, fields[col["subject_label"]])
    if not annotations:
        raise ValueError("Monarch's association table yielded no disease annotation")
    return annotations, labels, stats


def ensure_sources(config: dict) -> tuple:
    """Download (or reuse) the HPO and Monarch files. Returns ``(paths, provenance)``."""
    settings = ((config.get("l2") or {}).get("channel_d") or {})
    cache = refcache.reference_dir(config) / "phenotype"
    cache.mkdir(parents=True, exist_ok=True)
    paths, provenance = {}, {}
    for key, url, licence in (
        ("hpo", settings.get("hpo_url", HPO_URL),
         "HPO licence: free to use, cite, show the version, do not alter"),
        ("monarch", settings.get("monarch_url", MONARCH_URL),
         "Monarch KG association table; disease labels Mondo (CC BY 4.0); HPO annotations "
         "under the HPO licence"),
    ):
        # The release is in the URL; keep it in the cached name so two releases never
        # overwrite each other.
        release = next((p for p in url.split("/") if re.fullmatch(r"\d{4}-\d{2}-\d{2}", p)),
                       "unversioned")
        dest = cache / f"{release}.{url.rsplit('/', 1)[-1]}"
        paths[key] = dest
        provenance[key] = {"licence": licence, "release": release, **refcache.fetch(url, dest)}
    return paths, provenance


def open_text(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8")
