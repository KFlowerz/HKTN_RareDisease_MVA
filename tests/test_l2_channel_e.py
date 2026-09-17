"""Tests for L2 channel E: the curated aneuploidy-stress prior, ranked by verified papers.

Every compound, identifier and abstract here is invented -- compounds named ``inventib``
and the like, PMIDs in the ``99xxxxxx`` range, DOIs under the reserved ``10.5555`` test
prefix. No test asserts anything about a real paper or a real drug, and none reaches the
network: the Europe PMC and Open Targets calls are stubbed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.l2_channels import channel_e_prior as channel
from src.l2_channels import enrichment, literature

PRIOR_HEADER = "\t".join(channel.PRIOR_COLUMNS)


def _prior(tmp_path: Path, *rows: str) -> Path:
    path = tmp_path / "prior.tsv"
    path.write_text("# an invented seed file\n" + PRIOR_HEADER + "\n"
                    + "".join(row + "\n" for row in rows), encoding="utf-8")
    return path


def _row(compound="inventib", aliases="", target="Invented kinase",
         axis="invented stress", direction="exploit", caution="",
         anchors="DOI:10.5555/invented.1") -> str:
    return "\t".join([compound, aliases, target, axis, direction, caution, anchors])


def _record(pmid="99000001", title="Inventib is selective for aneuploid cells.",
            abstract="Inventib selectively killed aneuploid cells.", mesh=("cell line",),
            pub_types=("journal article",), **flags) -> literature.Record:
    return literature.Record(pmid=pmid, doi=f"10.5555/invented.{pmid}", title=title,
                             abstract=abstract, mesh=mesh, pub_types=pub_types,
                             year="2020", journal="Journal of Invented Results", **flags)


# ------------------------------------------------------------------- the seed file


def test_loads_the_committed_prior():
    """The real seed file must parse; it is an input to a gate, not a sample."""
    prior = channel.load_prior()
    assert prior, "the committed prior is empty"
    for row in prior:
        assert row["compound"] and row["axis"] and row["anchors"], row
        for citation in row["anchors"]:
            assert citation.startswith(("DOI:", "PMID:")), citation
        assert row["direction"] in {"exploit", "buffer"}, row


def test_the_committed_prior_contains_no_patient_data():
    """Constraint 7: no phenotype term, coordinate or identifier in a committed file."""
    literature.refuse_private(channel.PRIOR_FILE.read_text(encoding="utf-8"))


def test_a_misspelled_column_is_refused(tmp_path):
    path = tmp_path / "prior.tsv"
    path.write_text(PRIOR_HEADER.replace("axis", "axes") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header"):
        channel.load_prior(path)


def test_a_row_with_a_missing_tab_is_refused(tmp_path):
    """Silently short rows would shift every field and mis-assign a compound's evidence."""
    with pytest.raises(ValueError, match="field"):
        channel.load_prior(_prior(tmp_path, "inventib\t\tInvented kinase"))


def test_a_duplicated_compound_is_refused(tmp_path):
    with pytest.raises(ValueError, match="more than once"):
        channel.load_prior(_prior(tmp_path, _row(), _row()))


def test_aliases_and_anchors_split_on_semicolons(tmp_path):
    row = channel.load_prior(_prior(tmp_path, _row(
        aliases="inv-1;inv 1", anchors="DOI:10.5555/a;PMID:99000002")))[0]
    assert row["aliases"] == ["inv-1", "inv 1"]
    assert row["anchors"] == ["DOI:10.5555/a", "PMID:99000002"]


# ----------------------------------------------------------------- query building


def test_the_query_is_built_from_the_seed_file_and_config_only(tmp_path):
    row = channel.load_prior(_prior(tmp_path, _row(aliases="inv-1")))[0]
    query = channel.build_query(row, channel.DEFAULTS)
    assert '"inventib"' in query and '"inv-1"' in query
    assert '"aneuploid*"' in query
    literature.refuse_private(query)


def test_an_empty_context_is_refused(tmp_path):
    """A compound name alone retrieves that drug's whole literature, aneuploid or not."""
    row = channel.load_prior(_prior(tmp_path, _row()))[0]
    with pytest.raises(ValueError, match="context_terms"):
        channel.build_query(row, {**channel.DEFAULTS, "context_terms": []})


def test_the_channel_never_reads_the_phenotype(monkeypatch, tmp_path):
    """The structural guarantee behind the query guard: no patient input is even opened.

    Channel D parses the subject's phenotype from ``data_dir``. This channel must not, so
    a query cannot carry one even by accident. Breaking the phenotype reader and the data
    directory leaves this channel working.
    """
    from src.l2_channels import phenotype

    monkeypatch.setattr(phenotype, "find_document", lambda config: (_ for _ in ()).throw(
        AssertionError("channel E opened the phenotype document")))
    config = _config(tmp_path, data_dir="/nonexistent-on-purpose")
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)
    assert (Path(config["results_dir"]) / "l2" / "channel_e_prior"
            / "candidates.tsv").is_file()


# -------------------------------------------------------------------- classifying


def test_a_retracted_record_never_supports_a_candidate():
    supporting, counts, notes = channel.classify(
        [(_record(retracted=True), "retrieved")], require_direction=False)
    assert supporting == [] and counts["n_withdrawn"] == 1
    assert notes["99000001"]["excluded"] == "withdrawn"


def test_a_retracted_curated_citation_is_dropped_too():
    """Curation exempts a record from the lexical checks, never from retraction."""
    supporting, counts, _ = channel.classify(
        [(_record(retracted=True), "curated")], require_direction=True)
    assert supporting == [] and counts["n_withdrawn"] == 1


def test_a_record_reporting_the_compound_causing_aneuploidy_is_excluded():
    """For a cancer-predisposition syndrome, the conservative reading wins."""
    adverse = _record(title="Inventib is aneugenic in human lymphocytes.",
                      abstract="Inventib induced micronuclei and chromosome aberrations.")
    supporting, counts, notes = channel.classify([(adverse, "retrieved")],
                                                 require_direction=False)
    assert supporting == [] and counts["n_adverse_direction"] == 1
    assert notes["99000001"]["adverse_patterns"]


def test_an_ambiguous_record_is_read_as_adverse():
    """Matching both a supporting and an adverse phrase must not nominate a drug."""
    both = _record(abstract="Inventib selectively killed cells but was genotoxic.")
    supporting, counts, _ = channel.classify([(both, "retrieved")], require_direction=True)
    assert supporting == [] and counts["n_adverse_direction"] == 1


def test_a_retrieved_record_without_a_direction_phrase_is_excluded():
    bare = _record(title="Inventib pharmacokinetics in aneuploid cells.",
                   abstract="We measured plasma levels of inventib in aneuploid cells.")
    supporting, counts, _ = channel.classify([(bare, "retrieved")], require_direction=True)
    assert supporting == [] and counts["n_no_direction"] == 1


def test_a_curated_record_is_kept_without_a_direction_phrase():
    """A person read the paper and asserted the link; the abstract need not restate it."""
    bare = _record(title="Identification of selective antiproliferation compounds.",
                   abstract="A screen identified several compounds.")
    supporting, counts, _ = channel.classify(
        [(bare, "curated")],
        require_direction=True,
        cooccurrence=([literature.term_pattern("inventib")],
                      [literature.term_pattern("aneuploid*")]))
    assert len(supporting) == 1 and counts["n_supporting_curated"] == 1


def test_a_retrieved_record_mentioning_both_terms_apart_is_excluded():
    apart = _record(
        title="A case report.",
        abstract="The patient has trisomy X. Inventib was given for a selective effect.")
    supporting, counts, notes = channel.classify(
        [(apart, "retrieved")], require_direction=True,
        cooccurrence=([literature.term_pattern("inventib")],
                      [literature.term_pattern("trisomy")]))
    assert supporting == [] and counts["n_no_cooccurrence"] == 1
    assert notes[apart.key]["excluded"] == "no_sentence_cooccurrence"


def test_score_weights_the_grades_it_is_given():
    records = [_record(pmid="99000001", pub_types=("randomized controlled trial",)),
               _record(pmid="99000002", mesh=("mice",)),
               _record(pmid="99000003", mesh=("cell line",)),
               _record(pmid="99000004", mesh=())]
    total, by_grade = channel.score(records, channel.DEFAULTS["grade_weights"])
    assert by_grade == {"clinical": 1, "in_vivo": 1, "in_vitro": 1, "ungraded": 1}
    assert total == pytest.approx(4.0 + 2.0 + 1.0 + 0.25)


def test_an_unknown_grade_weight_is_refused(monkeypatch, tmp_path):
    """A typo in config would silently drop a grade's contribution to every rank."""
    config = _config(tmp_path)
    config["l2"]["channel_e"]["grade_weights"] = {"clinical": 1.0, "invitro": 1.0}
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    with pytest.raises(ValueError, match="unknown grade"):
        channel.generate(config)


# --------------------------------------------------------------- the channel itself


def _config(tmp_path: Path, **overrides) -> dict:
    config = {
        "seed": 42,
        "data_dir": str(tmp_path / "data"),
        "results_dir": str(tmp_path / "results"),
        "reference_dir": str(tmp_path / "reference"),
        "enrichment_dir": str(tmp_path / "enrichment"),
        "l2": {"channel_e": {"max_records_per_compound": 10}},
    }
    config.update(overrides)
    return config


def _stub(monkeypatch, tmp_path, *, retrieved, anchors=None, molecules=None,
          crossref=None, prior_rows=(_row(),)) -> None:
    """Replace the two outside dependencies: Europe PMC and the Open Targets molecules."""
    monkeypatch.setattr(channel, "PRIOR_FILE", _prior(tmp_path, *prior_rows))

    identities = molecules if molecules is not None else {
        "CHEMBL9000001": {"name": "INVENTIB", "drug_type": "Small molecule",
                          "clinical_stage": enrichment.APPROVED}}
    names = {enrichment.normalise_name(v["name"]): k for k, v in identities.items()}

    monkeypatch.setattr(channel.enrichment, "ensure_datasets",
                        lambda config, datasets, channel=None: (
                            {"drug_molecule": [tmp_path / "molecules.parquet"]},
                            [{"source": "invented", "release": "0.0"}]))
    monkeypatch.setattr(channel.enrichment, "load_molecules", lambda paths: identities)
    monkeypatch.setattr(channel.enrichment, "name_index", lambda paths: (
        names, {"molecule_names": len(names), "molecule_synonyms": 0,
                "ambiguous_synonyms": frozenset()}))

    monkeypatch.setattr(channel.literature, "search", lambda config, query, **kw: (
        list(retrieved), {"service": "invented", "query": query, "hit_count": len(retrieved),
                          "retrieved": "2026-09-17", "response_sha256": "0" * 64,
                          "records_returned": len(retrieved), "licence": "invented"}))
    monkeypatch.setattr(channel.literature, "fetch_identifier", lambda config, citation: (
        (anchors or {}).get(citation), {"service": "invented", "query": citation,
                                        "hit_count": 0, "retrieved": "2026-09-17",
                                        "response_sha256": "0" * 64,
                                        "records_returned": 0, "licence": "invented"}))
    # Stubbed too, or a resolved anchor would reach the real Crossref from a test.
    monkeypatch.setattr(channel.literature, "crossref", lambda config, doi: (
        crossref or {}).get(doi, {"resolved": True, "update_to": [], "doi": doi}))


def _outputs(config: dict) -> tuple:
    out = Path(config["results_dir"]) / "l2" / "channel_e_prior"
    return (out, json.loads((out / "channel.json").read_text(encoding="utf-8")),
            json.loads((out / "references.json").read_text(encoding="utf-8")))


def test_writes_a_ranked_table_with_its_references(monkeypatch, tmp_path):
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)

    out, meta, references = _outputs(config)
    rows = (out / "candidates.tsv").read_text(encoding="utf-8").splitlines()
    assert rows[0].split("\t")[:2] == ["rank", "chembl_id"]
    assert len(rows) == 2
    assert "INVENTIB" in rows[1]
    assert meta["counts"]["compounds_ranked"] == 1
    assert references["compounds"][0]["references"][0]["pmid"] == "99000001"


def test_every_candidate_carries_the_chemoprevention_endpoint(monkeypatch, tmp_path):
    """G2 decision D4: L3 and L5 must not be able to read this as symptomatic relief."""
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)
    out, _, _ = _outputs(config)
    header, row = (out / "candidates.tsv").read_text(encoding="utf-8").splitlines()[:2]
    assert row.split("\t")[header.split("\t").index("endpoint")] == "chemoprevention"


def test_no_abstract_text_reaches_any_output(monkeypatch, tmp_path):
    """Third-party content stays in the local cache."""
    secret = "A proprietary abstract that must not be republished."
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record(
        abstract=f"Inventib selectively killed aneuploid cells. {secret}")])
    channel.generate(config)
    out, _, _ = _outputs(config)
    for path in out.rglob("*"):
        if path.is_file():
            assert secret not in path.read_text(encoding="utf-8"), path


def test_an_unresolvable_curated_citation_is_dropped_and_reported(monkeypatch, tmp_path):
    """A fabricated reference is worse than a missing candidate."""
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()], anchors={})
    channel.generate(config)
    _, meta, references = _outputs(config)
    assert meta["counts"]["curated_citations_unresolved"] == 1
    assert references["compounds"][0]["curated_citations"][0]["status"] == "unresolved"
    assert meta["per_compound"][0]["n_curated_citations_verified"] == 0


def test_a_retraction_crossref_reports_drops_the_citation(monkeypatch, tmp_path):
    """Europe PMC is the primary retraction signal; Crossref is the second opinion.

    Checked 2026-09-17: Crossref's ``update-to`` was empty for a record Europe PMC reports
    as retracted, so Crossref adds only what it knows and never subtracts. When it *does*
    report one, the citation goes.
    """
    config = _config(tmp_path)
    anchor = _record(pmid="99000800")
    _stub(monkeypatch, tmp_path, retrieved=[],
          anchors={"DOI:10.5555/invented.1": anchor},
          crossref={anchor.doi: {"resolved": True, "update_to": [
              {"type": "retraction", "doi": "10.5555/invented.retraction"}]}})
    with pytest.raises(ValueError, match="unmet threshold"):
        channel.generate(config)


def test_a_citation_crossref_cannot_resolve_is_still_reported(monkeypatch, tmp_path):
    """A Crossref outage must not drop a record Europe PMC already established exists."""
    config = _config(tmp_path)
    anchor = _record(pmid="99000801")
    _stub(monkeypatch, tmp_path, retrieved=[],
          anchors={"DOI:10.5555/invented.1": anchor},
          crossref={anchor.doi: {"resolved": False, "update_to": []}})
    channel.generate(config)
    _, meta, references = _outputs(config)
    reported = references["compounds"][0]["curated_citations"][0]
    assert reported["status"] == "verified" and reported["crossref_resolved"] is False
    assert meta["counts"]["compounds_ranked"] == 1


def test_an_unapproved_compound_is_excluded_with_its_reason(monkeypatch, tmp_path):
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()], molecules={
        "CHEMBL9000001": {"name": "INVENTIB", "drug_type": "Small molecule",
                          "clinical_stage": "PHASE_2"}})
    with pytest.raises(ValueError, match="no compound reached"):
        channel.generate(config)


def test_an_unresolvable_compound_name_is_excluded_not_guessed(monkeypatch, tmp_path):
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()], molecules={})
    with pytest.raises(ValueError, match="no compound reached"):
        channel.generate(config)


def test_an_empty_result_is_an_error_not_an_empty_table(monkeypatch, tmp_path):
    """A stub that returns [] looks like a finding of 'no candidates' (CLAUDE.md)."""
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[], anchors={})
    with pytest.raises(ValueError, match="unmet threshold"):
        channel.generate(config)


def test_the_channel_record_says_what_rests_on_curation_alone(monkeypatch, tmp_path):
    config = _config(tmp_path)
    anchor = _record(pmid="99000900", title="A curated paper.", abstract="Nothing lexical.")
    _stub(monkeypatch, tmp_path, retrieved=[],
          anchors={"DOI:10.5555/invented.1": anchor})
    channel.generate(config)
    _, meta, _ = _outputs(config)
    assert meta["counts"]["compounds_supported_only_by_curation"] == 1
    assert meta["counts"]["compounds_ranked"] == 1


def test_the_channel_record_carries_the_queries_and_their_response_hashes(monkeypatch,
                                                                          tmp_path):
    """The index moves under a fixed query, so the response is pinned by hash."""
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)
    _, meta, _ = _outputs(config)
    assert meta["queries"] and all("query" in q and "response_sha256" in q
                                   for q in meta["queries"])


def test_the_caveats_and_parameters_travel_with_the_result(monkeypatch, tmp_path):
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)
    _, meta, _ = _outputs(config)
    assert set(meta["parameters"]) == set(channel.DEFAULTS)
    assert meta["parameters"]["grade_weights"] == channel.DEFAULTS["grade_weights"]
    assert any("curated, not discovered" in c for c in meta["caveats"])


def test_no_restricted_field_reaches_the_candidate_table(monkeypatch, tmp_path):
    """Open Targets drug content is ChEMBL-derived and stays in the enrichment zone (D7).

    What may leave is the whitelist in ``enrichment.PUBLISHABLE_FIELDS``. The curated
    ``target`` and ``axis`` columns come from the committed seed file, not from the
    restricted source, so they are this project's own text and may be published.
    """
    config = _config(tmp_path)
    _stub(monkeypatch, tmp_path, retrieved=[_record()])
    channel.generate(config)
    out, _, _ = _outputs(config)
    table = (out / "candidates.tsv").read_text(encoding="utf-8")
    header = table.splitlines()[0].split("\t")

    curated = set(channel.PRIOR_COLUMNS) | {"rank", "score", "endpoint", "drug_name"}
    counts = {c for c in header if c.startswith("n_")}
    allowed = curated | counts | set(enrichment.PUBLISHABLE_FIELDS)
    assert set(header) <= allowed, sorted(set(header) - allowed)
    # The restricted content itself: target sets and mechanism-of-action text.
    assert "mechanism_of_action" not in table and "ENSP" not in table


def test_the_module_does_not_import_the_phenotype_reader():
    """The structural half of the guarantee: there is no code path to patient input.

    Read from the parsed module rather than its text, so a mention in a docstring or a
    comment -- of which this module has several -- does not pass for an import.
    """
    import ast

    tree = ast.parse(Path(channel.__file__).read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[-1] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.update(alias.name.split(".")[-1] for alias in node.names)
            if node.module:
                imported.add(node.module.split(".")[-1])
    assert "phenotype" not in imported, "channel E must not reach the phenotype reader"
