"""Resolve and retraction-check citations against Crossref.

Purpose
    Turn a list of loose citation strings ("Author Year, Journal Vol:Pages") into
    resolved DOIs with a retraction/correction status, so they can enter
    ``docs/references.bib`` under the evidence rule in CLAUDE.md. Also re-checks
    already-resolved DOIs, which is required before submission -- retraction status
    changes between drafting and publication.

Inputs
    A text file, one citation per line. Lines beginning with ``#`` are ignored.
    A line may be either a bibliographic string or a bare DOI (``10.xxxx/...``).

Outputs
    A table on stdout, and ``--json`` for machine-readable output. Nothing is written
    to the repository: what gets committed is a human decision, because a fuzzy
    bibliographic match is a suggestion, not a verification.

Guardrail
    **A returned record is a candidate, not a confirmation.** Crossref's bibliographic
    query is fuzzy and will happily return the closest thing it has -- frequently a
    preprint when the journal version was cited, or an unrelated paper when the query
    is thin. This tool reports a similarity score and the fields it matched on; a human
    confirms. Automating the judgement away would reintroduce exactly the failure mode
    the citation rule exists to prevent.

    No credentials and no email are sent. Crossref's polite pool asks for a mailto
    header; supplying the maintainer's address to a third-party service is not this
    tool's call to make.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

CROSSREF = "https://api.crossref.org/works"
USER_AGENT = "mva-track2-refcheck/1.0 (https://github.com/KFlowerz/HKTN_RareDisease_MVA)"
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")

#: Crossref `update` types that mean the record is compromised, not merely amended.
BLOCKING_UPDATES = {"retraction", "withdrawal", "removal"}


@dataclass
class Result:
    """One citation's resolution outcome."""

    query: str
    doi: str | None = None
    title: str | None = None
    journal: str | None = None
    volume: str | None = None
    page: str | None = None
    year: int | None = None
    first_author: str | None = None
    score: float | None = None
    updates: list = field(default_factory=list)

    @property
    def status(self) -> str:
        if self.doi is None:
            return "UNRESOLVED"
        if any(u.get("type", "").lower() in BLOCKING_UPDATES for u in self.updates):
            return "RETRACTED"
        if self.updates:
            return "CORRECTED"
        return "ok"


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _updates_for(doi: str) -> list:
    """Return records that update (retract, correct, withdraw) this DOI."""
    url = f"{CROSSREF}?{urllib.parse.urlencode({'filter': f'updates:{doi}', 'rows': 5})}"
    try:
        items = _get(url)["message"].get("items", [])
    except Exception:
        return []
    out = []
    for it in items:
        for upd in it.get("update-to", []):
            if upd.get("DOI", "").lower() == doi.lower():
                out.append({"type": upd.get("type", "?"), "doi": it.get("DOI")})
    return out


def _fill(res: Result, item: dict) -> Result:
    res.doi = item.get("DOI")
    res.title = (item.get("title") or [None])[0]
    res.journal = (item.get("container-title") or [None])[0]
    res.volume = item.get("volume")
    res.page = item.get("page") or item.get("article-number")
    parts = (item.get("issued") or {}).get("date-parts") or [[None]]
    res.year = parts[0][0]
    authors = item.get("author") or []
    res.first_author = authors[0].get("family") if authors else None
    res.score = item.get("score")
    res.updates = _updates_for(res.doi) if res.doi else []
    return res


def resolve(citation: str, *, journal_only: bool = True) -> Result:
    """Resolve one citation string or DOI against Crossref.

    Args:
        citation: A bibliographic string, or a bare DOI.
        journal_only: Restrict bibliographic search to ``type:journal-article``.

    Note:
        ``journal_only`` defaults to True because publishers register more than papers.
        AACR deposits its supplementary data as Crossref records titled "Data from
        <article title>", and those outrank the articles themselves in bibliographic
        search -- three citations in the first sweep of this project resolved to data
        deposits with no volume, page, or year. Preprints cause the same problem in the
        other direction. Without the filter the tool reports a confident match to the
        wrong object.
    """
    res = Result(query=citation)
    try:
        if DOI_RE.match(citation):
            item = _get(f"{CROSSREF}/{urllib.parse.quote(citation)}")["message"]
        else:
            params = {"query.bibliographic": citation, "rows": 1}
            if journal_only:
                params["filter"] = "type:journal-article"
            items = _get(f"{CROSSREF}?{urllib.parse.urlencode(params)}")["message"].get("items", [])
            if not items:
                return res
            item = items[0]
        return _fill(res, item)
    except Exception as exc:  # noqa: BLE001 - report, never crash mid-sweep
        res.title = f"<error: {type(exc).__name__}>"
        return res


def bibtex_for(doi: str) -> str:
    """Fetch Crossref's own BibTeX for a DOI, via content negotiation.

    Transcribing 20 author lists by hand introduces errors that a citation rule exists
    to prevent, so the publisher's registered metadata is the source. The caller still
    confirms the record is the right one -- this only removes the typing.
    """
    req = urllib.request.Request(
        f"https://api.crossref.org/works/{urllib.parse.quote(doi)}/transform/application/x-bibtex",
        headers={"User-Agent": USER_AGENT, "Accept": "application/x-bibtex"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8").strip()


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(description="Resolve and retraction-check citations")
    parser.add_argument("path", type=Path, help="text file, one citation or DOI per line")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    parser.add_argument("--bibtex", action="store_true", help="emit Crossref BibTeX for each resolved DOI")
    parser.add_argument("--any-type", action="store_true", help="do not restrict to journal articles")
    parser.add_argument("--delay", type=float, default=0.4, help="seconds between requests")
    args = parser.parse_args(argv)

    lines = [
        ln.strip()
        for ln in args.path.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]

    results = []
    for i, line in enumerate(lines, 1):
        print(f"  [{i}/{len(lines)}] {line[:70]}", file=sys.stderr)
        results.append(resolve(line, journal_only=not args.any_type))
        time.sleep(args.delay)

    if args.json:
        print(json.dumps([r.__dict__ for r in results], indent=2, default=str))
        return 0

    if args.bibtex:
        for r in results:
            if not r.doi:
                print(f"% UNRESOLVED: {r.query}")
                continue
            note = (
                f"DOI resolved via Crossref {time.strftime('%Y-%m-%d')}; "
                + ("no retraction or correction on record." if r.status == "ok" else f"STATUS: {r.status}.")
            )
            entry = bibtex_for(r.doi)
            print(entry[:-1].rstrip().rstrip(",") + f",\n  note = {{{note}}},\n}}\n")
            time.sleep(args.delay)
        return 0

    for r in results:
        flag = {"ok": "  ", "CORRECTED": "! ", "RETRACTED": "XX", "UNRESOLVED": "??"}[r.status]
        print(f"\n{flag} {r.query}")
        if r.doi:
            print(f"     -> {r.first_author} ({r.year}) {r.title}")
            print(f"        {r.journal} {r.volume or '-'}:{r.page or '-'}  doi:{r.doi}")
            print(f"        score={r.score:.1f} status={r.status}" if r.score else f"        status={r.status}")
            for u in r.updates:
                print(f"        !! {u['type']} by {u['doi']}")
        else:
            print("     -> NOT RESOLVED")

    n_ok = sum(1 for r in results if r.status == "ok")
    print(
        f"\n{n_ok}/{len(results)} resolved cleanly; "
        f"{sum(1 for r in results if r.status == 'UNRESOLVED')} unresolved, "
        f"{sum(1 for r in results if r.status in ('RETRACTED', 'CORRECTED'))} flagged.",
        file=sys.stderr,
    )
    print("Every match above is a CANDIDATE. Confirm each before it enters references.bib.", file=sys.stderr)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
