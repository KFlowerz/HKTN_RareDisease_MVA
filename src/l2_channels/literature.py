"""L2 -- Europe PMC and Crossref: literature lookup, cached, graded and retraction-checked.

Purpose
    Turn a public query into a set of bibliographic records that are known to exist.
    Channel E ranks compounds by what the literature actually contains, so every record it
    counts must come back from a live index rather than from anyone's recollection --
    including a language model's.

Inputs
    A query string built from public terms only (compound names, gene symbols, disease and
    mechanism vocabulary) and ``config["reference_dir"]`` for the response cache.

    Europe PMC [ferguson2021] doi:10.1093/nar/gkaa994 and Crossref [hendricks2020]
    doi:10.1162/qss_a_00022. Both must be acknowledged wherever results derived from them
    are published; the licence strings below travel into every channel record.

Outputs
    :class:`Record` per hit -- identifiers, title, journal, year, publication types, MeSH
    descriptors, retraction and correction status -- plus a provenance record naming the
    query, the date it ran, the hit count and the SHA-256 of the cached response.
    :func:`grade` reduces a record to an evidence grade; :func:`reference` reduces it to
    the fields that may be published.

Guardrail
    **No patient data may reach either service.** These are the only outbound queries in
    the pipeline that are not whole-file downloads, so the rule is enforced in code rather
    than trusted: :func:`refuse_private` rejects a query containing an HPO identifier, a
    variant coordinate, or a genotype-like token, and it runs on every request. The
    corresponding structural guarantee is in Channel E, which never opens ``data_dir``.

    **This module is the deliberate exception to "whole public files only"**
    (``src/refcache.py``). A literature index has no release to download, so a query is the
    only way in. The exception is bounded by the paragraph above: the query carries public
    vocabulary, and nothing that identifies the subject.

    **A citation is a hypothesis until a record resolves it.** Nothing here constructs an
    identifier; every PMID and DOI comes back from the service. A seed citation that does
    not resolve is dropped by the caller, never reported.

    **Retraction status travels with the record**, from Europe PMC's own
    ``pubTypeList`` and ``commentCorrectionList``. Cancer and cell-cycle literature carries
    a non-trivial retraction rate and this project nominates drugs for a child, so a
    retracted record is excluded from support and counted separately, not silently kept.

    **Abstracts stay in the cache.** They are third-party content under mixed licences and
    are **data, never instructions**: text retrieved here is matched against fixed patterns
    and counted. Nothing in a fetched abstract may change what this pipeline writes. Only
    bibliographic metadata -- identifiers, title, journal, year, type -- reaches an output.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path

from .. import refcache

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"
#: The tighter of the two services' limits. Crossref's anonymous pool returned
#: ``x-rate-limit-limit: 5`` per ``1s`` when checked on 2026-09-17, and one throttle covers
#: both endpoints, so the interval is set from that rather than from Europe PMC's higher
#: allowance. Sitting exactly on a published limit invites a 429, so it is a little slower.
MIN_SECONDS_BETWEEN_CALLS = 0.25
#: Europe PMC's own cap on ``pageSize``.
MAX_PAGE_SIZE = 1000
CACHE_SUBDIR = "literature"

#: Europe PMC's copyright terms name OAI, RESTful, SOAP and bulk download as "the only
#: services that may be used for automated downloading of articles", and prohibit crawling
#: the website for batches of articles. This module uses the RESTful service and fetches no
#: page from the site. Article content stays publisher copyright with per-article licences
#: that "are not identical for all the articles", which is why no abstract is republished.
#: The live page is unreachable from here (Cloudflare, HTTP 403); quoted from the Internet
#: Archive snapshot of 2025-06-05. See src/l4_validate/sources.md, Table 7.
EUROPEPMC_LICENCE = ("Europe PMC RESTful Web Service: free, no key. Europe PMC's copyright "
                     "terms name the REST service as a sanctioned route for automated "
                     "retrieval and prohibit crawling the website; article content remains "
                     "publisher copyright under per-article licences, so this pipeline "
                     "publishes bibliographic metadata only and reproduces no abstract or "
                     "full text. EMBL-EBI terms of use also apply (read 2026-09-17): "
                     "attribution expected, provided AS IS without warranty.")
#: Crossref, read 2026-09-17: "Almost all of the metadata we hold is reusable without
#: restriction, with the exception of abstracts which are subject to publisher or author
#: copyright... considered to be 'facts' which are not copyrightable and are thus in the
#: public domain (CC0)."
CROSSREF_LICENCE = ("Crossref REST API (terms read 2026-09-17): metadata is reusable "
                    "without restriction and Crossref-generated data is released as public "
                    "domain (CC0); abstracts remain under publisher or author copyright "
                    "and are not used here")

#: Patterns that must never appear in an outbound query. Deliberately broader than
#: "phenotype": the point is that nothing derived from this child leaves the machine, and a
#: check that only knew about HPO would miss a coordinate pasted into a query by mistake.
PRIVATE_PATTERNS = (
    (re.compile(r"\bHP:\d{7}\b", re.I), "an HPO term identifier"),
    (re.compile(r"\bOMIM:?\d{6}\b", re.I), "an OMIM identifier"),
    (re.compile(r"\b(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)[:-]\d{3,}\b", re.I),
     "a genomic coordinate"),
    (re.compile(r"\b[ACGT]{1,}>[ACGT]{1,}\b"), "a variant allele change"),
    (re.compile(r"\b(?:c|p|g|m|n)\.[0-9A-Za-z_*+>-]{3,}\b"), "an HGVS expression"),
    (re.compile(r"\brs\d{3,}\b", re.I), "a dbSNP identifier"),
)

#: Europe PMC publication types that mark a record as withdrawn from the literature.
RETRACTED_PUB_TYPES = {"retracted publication", "retraction of publication"}
CONCERN_PUB_TYPES = {"expression of concern"}
#: ``commentCorrectionList`` types pointing *from* this record *to* a notice about it.
RETRACTION_NOTICES = {"retraction in"}
CONCERN_NOTICES = {"expression of concern in"}
CORRECTION_NOTICES = {"erratum in", "corrected and republished in"}

#: Evidence grades, strongest first. The vocabulary is fixed by CLAUDE.md: a grade travels
#: with every claim into L3 and the report.
GRADES = ("clinical", "in_vivo", "in_vitro", "ungraded")

CLINICAL_PUB_TYPES = {
    "clinical trial", "randomized controlled trial", "controlled clinical trial",
    "clinical trial, phase i", "clinical trial, phase ii", "clinical trial, phase iii",
    "clinical trial, phase iv", "pragmatic clinical trial", "multicenter study",
    "observational study", "meta-analysis", "case reports", "clinical study",
}
CLINICAL_MESH = {
    "clinical trials as topic", "randomized controlled trials as topic",
    "treatment outcome", "antineoplastic combined chemotherapy protocols",
    "maximum tolerated dose", "drug administration schedule",
}
#: NLM assigns "Animals" to vertebrate and invertebrate animal research alike; it is the
#: single descriptor that separates whole-organism work from cell culture.
IN_VIVO_MESH = {"animals", "mice", "rats", "disease models, animal",
                "xenograft model antitumor assays", "heterografts", "zebrafish",
                "drosophila melanogaster", "caenorhabditis elegans", "mice, nude",
                "mice, inbred c57bl", "mice, scid"}
IN_VITRO_MESH = {"cell line", "cell line, tumor", "cells, cultured", "in vitro techniques",
                 "hela cells", "hek293 cells", "mcf-7 cells", "fibroblasts",
                 "saccharomyces cerevisiae", "schizosaccharomyces pombe",
                 "cell culture techniques", "tumor cells, cultured",
                 "saccharomyces cerevisiae proteins", "primary cell culture"}


@dataclass(frozen=True)
class Record:
    """One bibliographic record as Europe PMC returned it.

    ``abstract`` is held for local pattern matching only and is never published -- see the
    module guardrail and :func:`reference`.
    """

    pmid: str = ""
    pmcid: str = ""
    doi: str = ""
    title: str = ""
    journal: str = ""
    year: str = ""
    source: str = ""
    pub_types: tuple = ()
    mesh: tuple = ()
    abstract: str = ""
    retracted: bool = False
    expression_of_concern: bool = False
    corrected: bool = False

    @property
    def key(self) -> str:
        """Stable identity for de-duplication: the PMID, else the DOI, else the title."""
        return self.pmid or self.doi or self.title.lower()

    @property
    def withdrawn(self) -> bool:
        """Retracted or under an expression of concern -- excluded from support either way."""
        return self.retracted or self.expression_of_concern


def refuse_private(query: str) -> None:
    """Raise if a query contains anything derived from the subject.

    Raises:
        ValueError: On an HPO id, an OMIM id, a genomic coordinate, an allele change, an
            HGVS expression or a dbSNP id. The message names the category, never the text
            that matched -- an error string is the one place patient data would escape a
            guard designed to stop exactly that.
    """
    for pattern, description in PRIVATE_PATTERNS:
        if pattern.search(query):
            raise ValueError(
                f"refusing to send a literature query containing {description}. Queries in "
                "this channel are built from public vocabulary only -- compound names, gene "
                "symbols, mechanism terms -- and nothing derived from the subject leaves "
                "this machine (COMPLIANCE.md)."
            )


def cache_dir(config: dict) -> Path:
    """``reference_dir/literature`` -- cached API responses, outside the repository."""
    path = refcache.reference_dir(config) / CACHE_SUBDIR
    path.mkdir(parents=True, exist_ok=True)
    return path


_last_call = 0.0


def _throttle() -> None:
    global _last_call
    wait = MIN_SECONDS_BETWEEN_CALLS - (time.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.monotonic()


def _get_json(url: str, timeout: int = 60) -> dict:
    _throttle()
    request = urllib.request.Request(url, headers=refcache.USER_AGENT)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _cached(path: Path, fetch) -> tuple:
    """Return ``(payload, retrieved)``, fetching and caching on a miss.

    A cached response is reused as-is so a re-run reproduces the run that produced a
    result: the date in the provenance record is the date the bytes were retrieved, not
    the date they were read.
    """
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        return stored["payload"], stored["retrieved"]
    payload = fetch()
    retrieved = date.today().isoformat()
    path.write_text(json.dumps({"retrieved": retrieved, "fetched_utc":
                                datetime.now(timezone.utc).isoformat(),
                                "payload": payload}, sort_keys=True), encoding="utf-8")
    return payload, retrieved


def _text(value) -> str:
    return "" if value is None else str(value)


def _listed(container, key) -> list:
    """Europe PMC wraps every list as ``{"<key>List": {"<key>": [...]}}``."""
    if not isinstance(container, dict):
        return []
    inner = container.get(key)
    if isinstance(inner, list):
        return inner
    return [inner] if inner else []


def parse_record(raw: dict) -> Record:
    """Build a :class:`Record` from one Europe PMC ``resultList`` entry."""
    pub_types = tuple(sorted({_text(t).lower() for t in
                              _listed(raw.get("pubTypeList"), "pubType") if t}))
    mesh = tuple(sorted({_text(h.get("descriptorName")).lower() for h in
                         _listed(raw.get("meshHeadingList"), "meshHeading")
                         if isinstance(h, dict) and h.get("descriptorName")}))
    notices = {_text(c.get("type")).lower() for c in
               _listed(raw.get("commentCorrectionList"), "commentCorrection")
               if isinstance(c, dict)}
    journal = ""
    info = raw.get("journalInfo")
    if isinstance(info, dict) and isinstance(info.get("journal"), dict):
        journal = _text(info["journal"].get("title"))
    return Record(
        pmid=_text(raw.get("pmid")),
        pmcid=_text(raw.get("pmcid")),
        doi=_text(raw.get("doi")).lower(),
        title=_text(raw.get("title")).strip(),
        journal=journal,
        year=_text(raw.get("pubYear")),
        source=_text(raw.get("source")),
        pub_types=pub_types,
        mesh=mesh,
        abstract=_text(raw.get("abstractText")),
        retracted=bool(set(pub_types) & RETRACTED_PUB_TYPES
                       or notices & RETRACTION_NOTICES),
        expression_of_concern=bool(set(pub_types) & CONCERN_PUB_TYPES
                                   or notices & CONCERN_NOTICES),
        corrected=bool(notices & CORRECTION_NOTICES),
    )


def search(config: dict, query: str, *, max_records: int = 200,
           page_size: int = 100) -> tuple:
    """Run a Europe PMC search. Returns ``(records, provenance)``.

    Args:
        config: Pipeline configuration; uses ``reference_dir``.
        query: Europe PMC query syntax, public vocabulary only.
        max_records: Stop after this many hits, so one broad query cannot pull the index.
        page_size: Hits per request, capped at :data:`MAX_PAGE_SIZE`.

    Raises:
        ValueError: If the query would carry patient-derived content
            (:func:`refuse_private`).
    """
    refuse_private(query)
    page_size = max(1, min(int(page_size), MAX_PAGE_SIZE, int(max_records)))
    request_key = json.dumps({"endpoint": EUROPEPMC, "query": query,
                              "max_records": int(max_records), "page_size": page_size,
                              "resultType": "core"}, sort_keys=True)
    digest = hashlib.sha256(request_key.encode("utf-8")).hexdigest()
    path = cache_dir(config) / f"europepmc.{digest}.json"

    def fetch() -> dict:
        pages, cursor, hit_count = [], "*", None
        while True:
            url = f"{EUROPEPMC}?" + urllib.parse.urlencode({
                "query": query, "resultType": "core", "format": "json",
                "pageSize": page_size, "cursorMark": cursor})
            page = _get_json(url)
            pages.append(page)
            hit_count = page.get("hitCount", 0)
            results = _listed(page.get("resultList"), "result")
            nxt = page.get("nextCursorMark")
            if (not results or not nxt or nxt == cursor
                    or sum(len(_listed(p.get("resultList"), "result")) for p in pages)
                    >= min(int(max_records), hit_count)):
                break
            cursor = nxt
        return {"hit_count": hit_count, "pages": pages}

    payload, retrieved = _cached(path, fetch)

    records, seen = [], set()
    for page in payload["pages"]:
        for raw in _listed(page.get("resultList"), "result"):
            record = parse_record(raw)
            if record.key and record.key not in seen:
                seen.add(record.key)
                records.append(record)
            if len(records) >= int(max_records):
                break
    provenance = {
        "service": "Europe PMC RESTful Web Service",
        "endpoint": EUROPEPMC,
        "query": query,
        "retrieved": retrieved,
        "hit_count": payload["hit_count"],
        "records_returned": len(records),
        # The whole index moves under a fixed query, so the response is pinned by hash
        # rather than by a release number: this is the evidence the ranking was built on.
        # Hashing the payload rather than the cache file leaves out the retrieval date, so
        # two machines that got the same answer record the same hash and a reader can tell
        # "the index gave us different results" from "we ran on different days".
        "response_sha256": hashlib.sha256(
            json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest(),
        "licence": EUROPEPMC_LICENCE,
    }
    return records, provenance


def fetch_identifier(config: dict, identifier: str) -> tuple:
    """Resolve one PMID or DOI to a record. Returns ``(record | None, provenance)``.

    Used to verify a seed citation. A citation that returns nothing is unresolvable and
    the caller drops it -- a fabricated reference is worse than a missing candidate.
    """
    value = identifier.strip()
    if value.lower().startswith("pmid:"):
        query = f'EXT_ID:{value.split(":", 1)[1].strip()} AND SRC:"MED"'
    elif value.lower().startswith("doi:"):
        query = f'DOI:"{value.split(":", 1)[1].strip()}"'
    else:
        raise ValueError(f"citation {identifier!r} must be given as 'PMID:<id>' or "
                         "'DOI:<id>' so the identifier type is never guessed")
    records, provenance = search(config, query, max_records=1, page_size=1)
    return (records[0] if records else None), provenance


def crossref(config: dict, doi: str) -> dict:
    """Crossref metadata for a DOI, as an independent check that it resolves.

    Returns the fields worth cross-checking plus any ``update-to`` relation Crossref holds
    (retractions and corrections are recorded there when the publisher deposits them).
    An unresolvable DOI returns ``{"resolved": False}``; a service error is not fatal here,
    because Europe PMC has already established that the record exists.
    """
    digest = hashlib.sha256(f"crossref:{doi.lower()}".encode("utf-8")).hexdigest()
    path = cache_dir(config) / f"crossref.{digest}.json"

    def fetch() -> dict:
        try:
            message = _get_json(CROSSREF + urllib.parse.quote(doi, safe=""))["message"]
        except (urllib.error.HTTPError, urllib.error.URLError, KeyError, ValueError):
            return {"resolved": False}
        updates = [{"type": _text(u.get("type")), "doi": _text(u.get("DOI"))}
                   for u in (message.get("update-to") or []) if isinstance(u, dict)]
        titles = message.get("title") or []
        containers = message.get("container-title") or []
        parts = ((message.get("issued") or {}).get("date-parts") or [[]])[0]
        return {"resolved": True, "doi": _text(message.get("DOI")).lower(),
                "title": _text(titles[0]) if titles else "",
                "container": _text(containers[0]) if containers else "",
                "year": str(parts[0]) if parts else "",
                "type": _text(message.get("type")), "update_to": updates}

    payload, retrieved = _cached(path, fetch)
    return {**payload, "retrieved": retrieved, "licence": CROSSREF_LICENCE}


def grade(record: Record) -> str:
    """Evidence grade from publication type and MeSH indexing.

    The ladder is ``clinical`` > ``in_vivo`` > ``in_vitro`` > ``ungraded``, assigned from
    NLM's own indexing rather than from the text -- a coarse but reproducible rule that
    needs no model and no judgement call per paper.

    Its limits, which belong in every output that carries a grade:

    - MeSH indexing lags publication, so a recent paper is often ``ungraded``. That means
      *not indexed yet*, never *no evidence*.
    - A grade describes **what system was studied**, not whether the finding supports the
      compound. Direction is a separate question, and a lexical one here
      (:func:`directional`).
    - Non-animal model organisms (yeast, for instance) are graded with cell-based work. The
      distinction that matters downstream is whether a whole vertebrate or a patient was
      involved, not whether the system was a cell.
    """
    types, mesh = set(record.pub_types), set(record.mesh)
    if types & CLINICAL_PUB_TYPES or mesh & CLINICAL_MESH:
        return "clinical"
    if mesh & IN_VIVO_MESH:
        return "in_vivo"
    if mesh & IN_VITRO_MESH:
        return "in_vitro"
    return "ungraded"


#: Sentence boundary: a full stop, question or exclamation mark followed by whitespace.
#: Deliberately naive -- "et al." and "Fig. 3" split a sentence in two, which can only lose
#: a co-occurrence, never invent one.
_SENTENCE = re.compile(r"(?<=[.!?])\s+")


def term_pattern(term: str):
    """Compile a search term the way Europe PMC reads it: ``*`` is a suffix wildcard.

    Whitespace matches any run of whitespace, so a phrase broken across a line still
    matches. The same term therefore means the same thing in the query and in the local
    text check, which is the point: a record retrieved for a term should be testable
    against it.
    """
    escaped = re.escape(term.strip()).replace(r"\ ", r"\s+").replace(r"\*", r"\w*")
    return re.compile(escaped, re.I)


def cooccur(record: Record, left, right) -> bool:
    """True if some sentence of the title or abstract matches a pattern from each group.

    Retrieval only asks that both terms appear *somewhere* in a record, which is how a
    myeloma trial that names a drug in one paragraph and chromosomal instability in another
    ends up counted as evidence about the drug and that state. Requiring them in one
    sentence is a crude proximity test, but it is the difference between "this paper
    mentions both" and "this paper says something about both at once".

    The title counts as its own sentence, and an unindexed record with no abstract is
    tested on its title alone.
    """
    for sentence in [record.title] + _SENTENCE.split(record.abstract or ""):
        if not sentence:
            continue
        if any(p.search(sentence) for p in left) and any(p.search(sentence) for p in right):
            return True
    return False


def matches(record: Record, patterns) -> tuple:
    """Which of ``patterns`` (compiled) appear in the record's title or abstract.

    Returns the matching pattern names, never the matched text: the caller writes counts
    and pattern names into outputs, and third-party abstract text stays in the cache.
    """
    haystack = f"{record.title}\n{record.abstract}".lower()
    return tuple(name for name, pattern in patterns if pattern.search(haystack))


def reference(record: Record) -> dict:
    """The publishable view of a record: bibliographic metadata and status, no abstract."""
    return {"pmid": record.pmid, "doi": record.doi, "title": record.title,
            "journal": record.journal, "year": record.year, "source": record.source,
            "pub_types": list(record.pub_types), "grade": grade(record),
            "retracted": record.retracted,
            "expression_of_concern": record.expression_of_concern,
            "corrected": record.corrected}
