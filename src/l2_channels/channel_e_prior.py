"""L2 Channel E -- Aneuploidy-stress literature prior, verified against Europe PMC.

Purpose
    Capture what the published literature already knows about pressuring, tolerating or
    buffering aneuploidy -- biology that no structured resource encodes. Aneuploid cells
    carry characteristic proteotoxic, metabolic, mitotic and replication stress, and
    compounds tested against those states are the strongest prior this project has.

    The channel starts from a curated seed set (``aneuploidy_prior.tsv``) and ranks it by
    **verified, graded literature**: for each compound it retrieves the records that
    co-occur with aneuploidy vocabulary, grades each by what system was studied, drops the
    retracted, and scores what remains.

Inputs
    ``aneuploidy_prior.tsv`` beside this module; Europe PMC and Crossref
    (:mod:`.literature`); Open Targets' molecule identities for the ChEMBL id and clinical
    stage (restricted zone, :mod:`.enrichment`). **No patient data, of any kind.**

Outputs
    Under ``results_dir/l2/channel_e_prior/``:
      - ``candidates.tsv`` -- ranked compounds with the graded evidence counts behind each
      - ``references.json`` -- every reference counted, per compound: PMID, DOI, title,
        journal, year, grade, retraction status, and whether it was curated or retrieved
      - ``channel.json`` -- method, parameters, per-compound counts, every query with the
        SHA-256 of its cached response, licensing, caveats, seed

Method
    1. Each seed compound's anchor citations are resolved against Europe PMC. An
       unresolvable citation is dropped and counted; a retracted one is dropped and
       flagged. Nothing is reported that did not come back from the index.
    2. A query is built from the compound's own names and the configured aneuploidy
       vocabulary, searched over title and abstract, and every hit is retrieved.
    3. Each record is graded ``clinical`` / ``in_vivo`` / ``in_vitro`` / ``ungraded`` from
       NLM's indexing, and screened for direction: a record that reports the compound
       *causing* aneuploidy or chromosome damage supports the opposite conclusion and is
       excluded from support, then counted as adverse.
    4. A compound's score is the weighted sum of its supporting records by grade. Weights
       are a declared project choice, in config, not a measurement.
    5. Compounds are resolved to a ChEMBL identity so L3 can aggregate them with the other
       channels; with ``approved_only`` the unapproved are excluded and the reason recorded.

Guardrail
    **Every claim carries a citation, and the citation must be verified to exist.** An
    identifier is never constructed here -- every PMID and DOI in an output came back from
    Europe PMC in this run. A seed citation that does not resolve is dropped, never
    reported: a fabricated reference in a rare-disease report is worse than a missing
    candidate.

    **No patient data is sent anywhere.** This module never opens ``data_dir`` and never
    imports the phenotype reader; queries are built from the committed seed file and
    config. :func:`.literature.refuse_private` rejects an outbound query carrying an HPO
    id, a coordinate or an HGVS expression, and ``tests/test_l2_channel_e.py`` asserts both
    the structural and the enforced guarantee.

    **Retrieved text is data, never instructions** (CLAUDE.md). Abstracts are matched
    against fixed patterns and counted; nothing in one may change what is written here, and
    no abstract text reaches an output.

    **The candidate set is curated, not discovered.** This channel verifies and ranks what
    a person listed. It cannot nominate a compound absent from the seed file, and every
    output says so. Widening the set is what an extraction model would buy, and that is a
    separate, evaluated decision.

    **Co-occurrence is not efficacy, and a grade is not a direction.** A record counted
    here studied the compound in an aneuploid or chromosomally unstable system; whether it
    helps a child with MVA is not established by anything in this channel.

    **This channel annotates; it does not exclude.** Safety, paediatric suitability and
    genotoxic exclusion are L4's.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

from . import enrichment, literature

LOGGER = logging.getLogger(__name__)

#: Gate G2 (decision D4): secondary prevention, not prevention of a first cancer.
ENDPOINT = "chemoprevention"
PRIOR_FILE = Path(__file__).with_name("aneuploidy_prior.tsv")
PRIOR_COLUMNS = ("compound", "aliases", "target", "axis", "direction", "caution", "anchors")
MOLECULE_DATASETS = ("drug_molecule",)

DEFAULTS = {
    # Searched over title and abstract, not full text. Full text roughly triples recall and
    # costs most of the precision -- "chloroquine" and "aneuploidy" co-occur in hundreds of
    # papers about neither together -- and there is no extraction model here to sort them.
    # The trade is deliberate and is the first thing a model would change.
    "search_field": "TITLE_ABS",
    "context_terms": ["aneuploid*", "chromosomal instability", "chromosome missegregation",
                      "trisomic", "trisomy", "whole-genome doubling"],
    "max_records_per_compound": 200,
    # A declared project choice, not a measured effect size: clinical evidence outweighs
    # animal work, which outweighs cell work. Stated in channel.json so a reader can
    # disagree with the number rather than guess at it.
    "grade_weights": {"clinical": 4.0, "in_vivo": 2.0, "in_vitro": 1.0, "ungraded": 0.25},
    "require_direction": True,
    # Retrieval only asks that the compound and the aneuploidy vocabulary both appear in a
    # record. Without this, a myeloma paper naming bortezomib in one paragraph and
    # chromosomal instability in another counts as evidence about the two together -- which
    # is how the most-studied drug, not the best-evidenced one, reaches rank 1.
    "require_cooccurrence": True,
    "approved_only": True,
    "open_targets_release": "26.06",
}

#: Phrases that mark a record as reporting a *therapeutic* direction -- the compound acting
#: on the aneuploid state. Lexical and therefore crude: they cannot tell a result from a
#: proposal, and a paper that says "no selective effect" matches "selectiv" as readily as
#: one that found one. They filter; they do not interpret.
SUPPORT_PATTERNS = tuple((name, re.compile(pattern)) for name, pattern in (
    ("selectivity", r"selectiv"),
    ("sensitivity", r"\bsensitiv"),
    ("vulnerability", r"vulnerab"),
    ("synthetic_lethality", r"synthetic\s+lethal"),
    ("preferential_effect", r"preferential"),
    ("dependency", r"\bdepend(?:ence|ency|ent)\b"),
    ("antiproliferative", r"antiprolifer"),
    ("growth_inhibition", r"(?:growth|proliferation|viability)\s+\w*\s*"
                          r"(?:inhibit|suppress|arrest|reduc|impair)"),
))

#: Phrases that mark the *opposite* direction: the compound damaging chromosomes or
#: producing aneuploidy. A record matching one of these is excluded from support even if it
#: also matches a support pattern -- for a child with a cancer-predisposition syndrome, the
#: conservative reading of an ambiguous abstract is the one that does not nominate a drug.
ADVERSE_PATTERNS = tuple((name, re.compile(pattern)) for name, pattern in (
    ("aneugenic", r"aneugen"),
    ("genotoxic", r"genotoxic|clastogen"),
    ("induces_aneuploidy", r"induc\w*\s+(?:\w+\s+){0,3}(?:aneuploid|micronucle|"
                           r"chromosom\w+\s+(?:aberration|damage|loss))"),
    ("micronucleus_assay", r"micronucle(?:us|i)\s+(?:assay|test)"),
    ("chromosome_damage", r"chromosom\w+\s+(?:aberration|breakage|damage)"),
))

CAVEATS = (
    "The candidate set is curated, not discovered. This channel verifies, grades and ranks "
    "the compounds listed in aneuploidy_prior.tsv; it cannot nominate one that is absent "
    "from that file, and its ranking inherits whatever bias the curation carries.",
    "A counted record studied the compound in an aneuploid or chromosomally unstable "
    "system. That is co-occurrence of a compound and a state in one paper, not evidence of "
    "benefit, and certainly not evidence of benefit in a child with MVA.",
    "Direction is decided lexically, from fixed phrases in the title and abstract. The "
    "check cannot distinguish a finding from a proposal, nor a positive result from a "
    "reported absence of one. It is the weakest step in the channel and the one an "
    "extraction model would replace.",
    "A supporting record must name the compound and the aneuploid state in one sentence. "
    "That is proximity, not meaning: a sentence saying the two are unrelated counts the "
    "same as one saying they are. It removes records that merely mention both somewhere.",
    "clinical_stage is Open Targets' maximum clinical stage for the molecule across all "
    "indications and jurisdictions. It is not a statement that the drug is approved for "
    "this purpose, in this population, or anywhere in particular; paediatric suitability "
    "is L4's to assess and is not implied by a rank here.",
    "A curated citation is counted on the curator's assertion. This channel verifies that "
    "it exists and has not been withdrawn; it does not check that the paper says what the "
    "seed file claims, and the lexical direction and co-occurrence rules are not applied "
    "to it. The n_curated_citations_verified column says how much of a compound's score "
    "rests on that assertion.",
    "Evidence grades come from NLM indexing, not from reading the paper. MeSH indexing "
    "lags publication, so a recent record is often 'ungraded' -- which means not indexed "
    "yet, never no evidence.",
    "Search is over title and abstract only. A result reported in a paper's body but not "
    "its abstract is not counted, so a compound's count is a lower bound.",
    "Grade weights are a declared project choice, not a measured effect size. The ranking "
    "moves if they move; they are recorded in channel.json for exactly that reason.",
    "Retracted records and records under an expression of concern are excluded from "
    "support and counted separately. Retraction status is Europe PMC's, and a retraction "
    "not yet indexed there will not be caught.",
    "Most selectivity evidence behind this prior comes from transformed, often p53-null "
    "cancer cells. Constitutional MVA cells are non-transformed with impaired but present "
    "checkpoints, so the extrapolation to this proband is an assumption, not a result.",
    "Counts are not effect sizes. A compound studied often ranks above one studied rarely, "
    "whatever either study found.",
)


def _settings(config: dict) -> dict:
    return {**DEFAULTS, **((config.get("l2") or {}).get("channel_e") or {})}


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def load_prior(path: Path = PRIOR_FILE) -> list:
    """Parse the curated seed file. Returns one dict per compound.

    Raises:
        ValueError: On a missing or misspelled column, a duplicate compound, or a row with
            the wrong number of fields -- each of which would silently drop or merge a
            compound's evidence.
    """
    rows, header = [], None
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = tuple(f.strip() for f in fields)
                if header != PRIOR_COLUMNS:
                    raise ValueError(f"{path} header is {header}, expected {PRIOR_COLUMNS}")
                continue
            if len(fields) != len(PRIOR_COLUMNS):
                raise ValueError(f"{path} line {number} has {len(fields)} field(s), expected "
                                 f"{len(PRIOR_COLUMNS)}; a tab is missing or doubled")
            row = dict(zip(PRIOR_COLUMNS, (f.strip() for f in fields)))
            row["aliases"] = [a.strip() for a in row["aliases"].split(";") if a.strip()]
            row["anchors"] = [a.strip() for a in row["anchors"].split(";") if a.strip()]
            if not row["compound"]:
                raise ValueError(f"{path} line {number} has no compound name")
            rows.append(row)
    if header is None:
        raise ValueError(f"{path} has no header row")
    duplicates = sorted({r["compound"] for r in rows
                         if sum(1 for o in rows if o["compound"] == r["compound"]) > 1})
    if duplicates:
        raise ValueError(f"{path} lists {duplicates} more than once; one row per molecule, "
                         "or the same evidence is counted twice")
    return rows


def build_query(entry: dict, settings: dict) -> str:
    """The Europe PMC query for one compound: its names AND the aneuploidy vocabulary.

    Built only from the committed seed file and config. Nothing patient-derived can reach
    it, and :func:`.literature.refuse_private` re-checks that on the way out.
    """
    field = str(settings["search_field"]).strip()
    names = [entry["compound"], *entry["aliases"]]
    context = list(settings["context_terms"])
    if not context:
        raise ValueError("l2.channel_e.context_terms is empty; a compound name alone would "
                         "retrieve the whole of that drug's literature")

    def render(term: str) -> str:
        """Quote a phrase; leave a wildcard term bare.

        Europe PMC does **not** expand a wildcard inside quotes. Measured on 2026-09-17,
        ``TITLE_ABS:"aneuploid*"`` returns 9,711 hits -- exactly what ``TITLE_ABS:"aneuploid"``
        returns -- while the bare ``TITLE_ABS:aneuploid*`` returns 28,075. Quoting every term
        therefore silently deleted every wildcard in the config while the query still looked
        like it had one, and it made the retrieval mean something different from what
        :func:`.literature.term_pattern` tests locally. A phrase still needs its quotes:
        unquoted, ``chromosomal instability`` is read as two terms, not one phrase.
        """
        term = term.strip()
        if "*" in term:
            if any(c.isspace() for c in term):
                raise ValueError(
                    f"l2.channel_e term {term!r} has both a wildcard and a space. Europe PMC "
                    "expands a wildcard only on a bare single term and treats an unquoted "
                    "phrase as separate terms, so this cannot be expressed; split it into a "
                    "wildcard term and a quoted phrase.")
            return f"{field}:{term}"
        return f'{field}:"{term}"'

    def clause(terms):
        return "(" + " OR ".join(render(t) for t in terms) + ")"

    return f"{clause(names)} AND {clause(context)}"


def verify_anchors(config: dict, entry: dict) -> tuple:
    """Resolve a compound's curated citations. Returns ``(records, report, provenance)``.

    Every anchor is looked up in Europe PMC; the report says what happened to each, so a
    dropped citation is visible rather than merely absent.

    A resolved DOI is then checked against Crossref as an independent second opinion --
    that the DOI itself resolves, and whether the publisher deposited an ``update-to``
    relation marking a retraction. Checked on 2026-09-17, Crossref's ``update-to`` was
    empty for a record Europe PMC correctly reports as retracted, so Crossref is the second
    check and not the first: a retraction it reports is believed, and its silence proves
    nothing. A Crossref outage is not fatal, because Europe PMC has already established
    that the record exists.
    """
    records, report, provenance = [], [], []
    for citation in entry["anchors"]:
        record, prov = literature.fetch_identifier(config, citation)
        provenance.append({**prov, "compound": entry["compound"], "role": "curated citation"})
        if record is None:
            report.append({"citation": citation, "status": "unresolved"})
            continue

        crossref = literature.crossref(config, record.doi) if record.doi else {
            "resolved": False, "update_to": []}
        withdrawn_at_crossref = [u for u in crossref.get("update_to") or []
                                 if "retract" in u.get("type", "").lower()]
        status = "retracted" if record.retracted or withdrawn_at_crossref else (
            "expression_of_concern" if record.expression_of_concern else "verified")
        report.append({"citation": citation, "status": status, "pmid": record.pmid,
                       "doi": record.doi, "crossref_resolved": bool(crossref.get("resolved")),
                       "crossref_updates": crossref.get("update_to") or []})
        if status == "verified":
            records.append(record)
    return records, report, provenance


def classify(pooled, *, require_direction: bool, cooccurrence=None) -> tuple:
    """Split records into supporting and excluded. Returns ``(supporting, counts, notes)``.

    ``notes`` maps a record key to the pattern names that placed it, so ``references.json``
    can say why each record counted or did not.

    **Curated citations skip the lexical checks.** The direction and co-occurrence rules
    exist to triage records a keyword search returned, where nobody has looked at the paper.
    A curated anchor is the opposite case: a person read it and asserted the link, and the
    assertion does not depend on the abstract phrasing it in one sentence. Applying the
    filters to anchors discarded the best evidence in the file -- the paper that identified
    chloroquine as aneuploidy-selective never names the compound in its abstract. What
    still applies to an anchor is retraction, which is a fact about the record, and the
    adverse-direction check, because a citation that reports the compound damaging
    chromosomes is one the curator should have to defend.

    The cost is explicit and belongs in the caveats: a curated record counts on a person's
    assertion, verified only to exist and not to have been withdrawn.

    Args:
        pooled: ``(record, origin)`` pairs, ``origin`` being ``"curated"`` or
            ``"retrieved"``.
        require_direction: Drop a retrieved record with no therapeutic-direction phrase.
        cooccurrence: ``(compound patterns, context patterns)``, or ``None`` to skip the
            check. A retrieved record must name the compound and the aneuploid state in one
            sentence; retrieval only asks that both appear somewhere in the record, which
            counts a paper that mentions each in a different paragraph.
    """
    supporting, notes = [], {}
    counts = {"n_retrieved": sum(1 for _, o in pooled if o == "retrieved"),
              "n_curated": sum(1 for _, o in pooled if o == "curated"),
              "n_withdrawn": 0, "n_adverse_direction": 0, "n_no_cooccurrence": 0,
              "n_no_direction": 0, "n_supporting": 0, "n_supporting_curated": 0}
    counts["n_uncitable"] = 0
    for record, origin in pooled:
        adverse = literature.matches(record, ADVERSE_PATTERNS)
        support = literature.matches(record, SUPPORT_PATTERNS)
        notes[record.key] = {"support_patterns": list(support),
                             "adverse_patterns": list(adverse)}
        if not record.citable:
            # No PMID, no DOI, no Europe PMC id: a reader cannot look it up, so it cannot
            # carry a claim (CLAUDE.md, Evidence and citations).
            counts["n_uncitable"] += 1
            notes[record.key]["excluded"] = "no_resolvable_identifier"
            continue
        if record.withdrawn:
            counts["n_withdrawn"] += 1
            notes[record.key]["excluded"] = "withdrawn"
            continue
        if adverse:
            counts["n_adverse_direction"] += 1
            notes[record.key]["excluded"] = "adverse_direction"
            continue
        if origin == "curated":
            counts["n_supporting_curated"] += 1
        else:
            if cooccurrence and not literature.cooccur(record, *cooccurrence):
                counts["n_no_cooccurrence"] += 1
                notes[record.key]["excluded"] = "no_sentence_cooccurrence"
                continue
            if require_direction and not support:
                counts["n_no_direction"] += 1
                notes[record.key]["excluded"] = "no_direction_phrase"
                continue
        supporting.append(record)
    counts["n_supporting"] = len(supporting)
    return supporting, counts, notes


def score(records, weights: dict) -> tuple:
    """``(score, per-grade counts)`` for a compound's supporting records."""
    by_grade = {grade: 0 for grade in literature.GRADES}
    for record in records:
        by_grade[literature.grade(record)] += 1
    total = sum(float(weights.get(grade, 0.0)) * n for grade, n in by_grade.items())
    return round(total, 4), by_grade


def generate(config: dict) -> None:
    """Rank the curated prior by verified, graded literature.

    Args:
        config: Parsed pipeline configuration; uses ``reference_dir`` (response cache),
            ``enrichment_dir`` (molecule identities), ``results_dir``, ``seed`` and
            ``l2.channel_e``.

    Raises:
        ValueError: On an unusable seed file, or if no compound survives -- an unmet
            threshold reported as such, never as an empty table (``l2_channels/run.py``).
    """
    settings = _settings(config)
    weights = dict(settings["grade_weights"])
    unknown = sorted(set(weights) - set(literature.GRADES))
    if unknown:
        raise ValueError(f"l2.channel_e.grade_weights names unknown grade(s) {unknown}; "
                         f"the vocabulary is {list(literature.GRADES)}")
    out_dir = Path(config["results_dir"]) / "l2" / "channel_e_prior"
    out_dir.mkdir(parents=True, exist_ok=True)

    prior = load_prior(PRIOR_FILE)
    LOGGER.info("channel E: %d curated compound(s)", len(prior))

    ot_paths, molecule_provenance = enrichment.ensure_datasets(
        config, MOLECULE_DATASETS, channel="channel_e")
    molecules = enrichment.load_molecules(ot_paths["drug_molecule"])
    names, name_stats = enrichment.name_index(ot_paths["drug_molecule"])

    entries, excluded, examined, queries = [], [], [], []
    for row in prior:
        anchors, anchor_report, anchor_provenance = verify_anchors(config, row)
        queries.extend(anchor_provenance)

        query = build_query(row, settings)
        retrieved, provenance = literature.search(
            config, query, max_records=int(settings["max_records_per_compound"]))
        queries.append({**provenance, "compound": row["compound"]})

        # Anchors first, so a curated citation that the query also returns keeps its
        # origin. De-duplication is by PMID, then DOI, then title (Record.key).
        pooled, seen = [], set()
        for record, origin in ([(r, "curated") for r in anchors]
                              + [(r, "retrieved") for r in retrieved]):
            if record.key in seen:
                continue
            seen.add(record.key)
            pooled.append((record, origin))

        cooccurrence = None
        if bool(settings["require_cooccurrence"]):
            cooccurrence = ([literature.term_pattern(n) for n in
                             [row["compound"], *row["aliases"]]],
                            [literature.term_pattern(t) for t in settings["context_terms"]])
        supporting, counts, notes = classify(
            pooled, require_direction=bool(settings["require_direction"]),
            cooccurrence=cooccurrence)
        total, by_grade = score(supporting, weights)

        lookup = [enrichment.normalise_name(n) for n in [row["compound"], *row["aliases"]]]
        chembl_id = next((names[n] for n in lookup if n in names), "")
        # A name that two molecules share resolved to whichever the index kept. That is a
        # coin toss about drug identity, so it travels into the output rather than being
        # counted and forgotten.
        ambiguous = sorted(n for n in lookup
                           if n in names and n in name_stats["ambiguous_synonyms"])
        identity = molecules.get(chembl_id, {})
        record = enrichment.DrugRecord(chembl_id=chembl_id, **{
            k: identity.get(k, "") for k in ("name", "drug_type", "clinical_stage")})

        entry = {"row": row, "record": record, "score": total, "by_grade": by_grade,
                 "counts": counts, "notes": notes, "pooled": pooled,
                 "anchor_report": anchor_report, "query": query,
                 "ambiguous_names": ambiguous, "n_anchors_verified": len(anchors)}
        # Every compound is kept here, ranked or not: the channel-level counts below are
        # about the literature this run examined, and summing them over ranked compounds
        # only would undercount exactly the compounds whose citations were dropped.
        examined.append(entry)
        if not chembl_id:
            excluded.append({**_summary(entry), "reason": "no_chembl_identity"})
            continue
        if bool(settings["approved_only"]) and not record.approved:
            excluded.append({**_summary(entry), "reason": "not_approved",
                             "clinical_stage": record.clinical_stage})
            continue
        if not counts["n_supporting"]:
            excluded.append({**_summary(entry), "reason": "no_supporting_record"})
            continue
        if total <= 0:
            # Distinct from having no support: the records exist, but every one of them
            # landed in a grade the configured weights value at zero. That is a statement
            # about the weights, and mislabelling it as "no evidence" would hide it.
            excluded.append({**_summary(entry), "reason": "no_weighted_support"})
            continue
        entries.append(entry)

    entries.sort(key=lambda e: (-e["score"], -e["counts"]["n_supporting"],
                                e["record"].chembl_id))
    if not entries:
        raise ValueError(
            f"no compound reached channel E: {len(prior)} curated compound(s), "
            f"{len(excluded)} excluded (see channel.json). That is an unmet threshold on "
            "these inputs -- an empty prior, a licence-restricted dataset that did not "
            "load, or an approved_only filter that removed everything -- not a finding "
            "that no compound has literature behind it.")

    rows = []
    for position, entry in enumerate(entries, start=1):
        public = enrichment.publishable(entry["record"])
        by_grade = entry["by_grade"]
        rows.append([position, public["chembl_id"], public["name"], public["drug_type"],
                     public["clinical_stage"], f"{entry['score']:.4f}", entry["row"]["axis"],
                     entry["row"]["target"], entry["row"]["direction"],
                     entry["counts"]["n_supporting"], by_grade["clinical"],
                     by_grade["in_vivo"], by_grade["in_vitro"], by_grade["ungraded"],
                     entry["n_anchors_verified"], entry["counts"]["n_adverse_direction"],
                     entry["row"]["caution"], ENDPOINT])
    _write_tsv(out_dir / "candidates.tsv",
               ["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "score",
                "axis", "target", "direction", "n_supporting_records", "n_clinical",
                "n_in_vivo", "n_in_vitro", "n_ungraded", "n_curated_citations_verified",
                "n_adverse_direction_records", "caution", "endpoint"], rows)

    (out_dir / "references.json").write_text(json.dumps({
        "notice": ("Every identifier here was returned by Europe PMC in this run; none was "
                   "constructed. Records marked retracted or under an expression of concern "
                   "are listed and excluded from support. Abstract text is not reproduced."),
        "compounds": [{
            "compound": e["row"]["compound"],
            "chembl_id": e["record"].chembl_id,
            "query": e["query"],
            "curated_citations": e["anchor_report"],
            "references": [{**literature.reference(r), "origin": origin,
                            **e["notes"].get(r.key, {})} for r, origin in e["pooled"]],
        } for e in entries],
        "seed": config["seed"],
    }, indent=2, sort_keys=True), encoding="utf-8")

    (out_dir / "channel.json").write_text(json.dumps({
        "channel": "channel_e_prior",
        "endpoint": ENDPOINT,
        "method": {
            "seed_set": "curated aneuploidy-stress compounds, src/l2_channels/"
                        "aneuploidy_prior.tsv",
            "retrieval": f"Europe PMC search over {settings['search_field']}: the compound's "
                         "names AND the configured aneuploidy vocabulary",
            "verification": "every curated citation resolved against Europe PMC; "
                            "unresolvable citations dropped and reported",
            "grading": "clinical / in_vivo / in_vitro / ungraded from NLM publication types "
                       "and MeSH descriptors",
            "direction": "fixed lexical patterns over title and abstract; a record matching "
                         "an adverse pattern is excluded from support even if it also "
                         "matches a supporting one",
            "scoring": "weighted sum of supporting records by grade",
            "ranking": "descending score, then supporting-record count, then ChEMBL id",
        },
        "parameters": {**{k: settings[k] for k in DEFAULTS}},
        "counts": {
            "compounds_curated": len(prior),
            "compounds_ranked": len(entries),
            "compounds_excluded": len(excluded),
            "records_supporting": sum(e["counts"]["n_supporting"] for e in entries),
            # Over every compound examined, not just the ranked ones. A retracted citation
            # or a genotoxicity paper matters most for a compound that did not make the
            # table, and counting only survivors would hide it.
            "records_withdrawn": sum(e["counts"]["n_withdrawn"] for e in examined),
            "records_adverse_direction": sum(e["counts"]["n_adverse_direction"]
                                             for e in examined),
            "records_without_a_resolvable_identifier": sum(
                e["counts"]["n_uncitable"] for e in examined),
            "curated_citations_unresolved": sum(
                1 for e in examined for a in e["anchor_report"]
                if a["status"] == "unresolved"),
            "curated_citations_retracted": sum(
                1 for e in examined for a in e["anchor_report"]
                if a["status"] in ("retracted", "expression_of_concern")),
            "compounds_resolved_by_an_ambiguous_synonym": sum(
                1 for e in entries if e["ambiguous_names"]),
            # The headline honesty number: for these compounds the search added nothing the
            # filters kept, so the rank rests entirely on the curated citation. A reader
            # comparing this with compounds_ranked sees how much of the table the retrieval
            # actually earned.
            "compounds_supported_only_by_curation": sum(
                1 for e in entries
                if e["counts"]["n_supporting"] == e["counts"]["n_supporting_curated"]),
            **{k: v for k, v in name_stats.items() if k != "ambiguous_synonyms"},
        },
        "per_compound": [_summary(e) for e in entries],
        "excluded": excluded,
        "queries": queries,
        "sources": {"literature": literature.EUROPEPMC_LICENCE,
                    "crossref": literature.CROSSREF_LICENCE,
                    "molecules": molecule_provenance},
        "licensing": {
            "published_fields": list(enrichment.PUBLISHABLE_FIELDS),
            "note": ("Drug identity comes from Open Targets and is ChEMBL-derived, so it "
                     "stays in the enrichment zone (decision D7) and only whitelisted "
                     "fields reach candidates.tsv. Bibliographic metadata is "
                     "redistributable; abstracts are not reproduced and stay in the local "
                     "response cache. See src/l4_validate/sources.md."),
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }, indent=2, sort_keys=True), encoding="utf-8")

    LOGGER.info("channel E written: %d compound(s) ranked, %d excluded, %d supporting "
                "record(s)", len(entries), len(excluded),
                sum(e["counts"]["n_supporting"] for e in entries))


def _summary(entry: dict) -> dict:
    """The per-compound record written into ``channel.json`` -- counts and identity only."""
    return {"compound": entry["row"]["compound"], "axis": entry["row"]["axis"],
            "chembl_id": entry["record"].chembl_id,
            "clinical_stage": entry["record"].clinical_stage,
            "score": entry["score"], "by_grade": entry["by_grade"], **entry["counts"],
            "ambiguous_names": entry["ambiguous_names"],
            "n_curated_citations": len(entry["anchor_report"]),
            "n_curated_citations_verified": entry["n_anchors_verified"]}
