"""Tests for the Europe PMC / Crossref lookup behind L2 channel E.

Nothing here touches the network: every test either stubs the HTTP call or works on a
parsed record. The identifiers are invented -- PMIDs in the ``99xxxxxx`` range and DOIs
under the reserved ``10.5555`` test prefix -- so no test asserts anything about a real
paper, and a test can never be mistaken for a verified citation.
"""

from __future__ import annotations

import json

import pytest

from src.l2_channels import literature


def _raw(**overrides) -> dict:
    """One Europe PMC ``resultList`` entry, in the service's own shape."""
    raw = {
        "id": "99000001", "source": "MED", "pmid": "99000001",
        "doi": "10.5555/invented.1", "title": "An invented study of something.",
        "pubYear": "2020", "abstractText": "An invented abstract.",
        "journalInfo": {"journal": {"title": "Journal of Invented Results"}},
        "pubTypeList": {"pubType": ["Journal Article"]},
        "meshHeadingList": {"meshHeading": [{"descriptorName": "Cell Line"}]},
    }
    raw.update(overrides)
    return raw


# --------------------------------------------------------------------------- guardrail


@pytest.mark.parametrize("query, what", [
    ('TITLE_ABS:"HP:0001250"', "an HPO term"),
    ('TITLE_ABS:"chr15:40000000"', "a coordinate"),
    ('"15-40161020"', "a bare coordinate"),
    ('TITLE_ABS:"c.2211_2212insA"', "an HGVS change"),
    ('"rs121913529"', "a dbSNP id"),
    ('"A>G" AND "BUB1B"', "an allele change"),
    ('"OMIM:257300"', "an OMIM id"),
])
def test_refuses_a_query_carrying_patient_derived_content(query, what):
    """Nothing about the subject may leave the machine, so the check is in code."""
    with pytest.raises(ValueError, match="refusing to send"):
        literature.refuse_private(query)


def test_the_refusal_message_does_not_repeat_the_offending_text():
    """An error string is the one place patient data would escape the guard."""
    with pytest.raises(ValueError) as caught:
        literature.refuse_private('TITLE_ABS:"HP:0004322"')
    assert "HP:0004322" not in str(caught.value)


@pytest.mark.parametrize("query", [
    '(TITLE_ABS:"chloroquine") AND (TITLE_ABS:"aneuploid*")',
    '(TITLE_ABS:"BUB1B" OR TITLE_ABS:"BUBR1") AND (TITLE_ABS:"trisomy")',
    'EXT_ID:21315436 AND SRC:"MED"',
])
def test_allows_a_query_of_public_vocabulary(query):
    literature.refuse_private(query)


def test_search_checks_the_query_before_any_request(monkeypatch, tmp_path):
    """The guard must run before the socket opens, not after the response comes back."""
    def explode(*args, **kwargs):
        raise AssertionError("a request was made despite a refused query")

    monkeypatch.setattr(literature, "_get_json", explode)
    with pytest.raises(ValueError, match="refusing to send"):
        literature.search({"reference_dir": str(tmp_path)}, 'TITLE_ABS:"HP:0001250"')


# ------------------------------------------------------------------------- parsing


def test_parses_the_fields_an_output_will_carry():
    record = literature.parse_record(_raw())
    assert (record.pmid, record.doi, record.year) == (
        "99000001", "10.5555/invented.1", "2020")
    assert record.journal == "Journal of Invented Results"
    assert record.pub_types == ("journal article",)
    assert record.mesh == ("cell line",)
    assert not record.withdrawn


def test_a_retracted_publication_type_marks_the_record():
    record = literature.parse_record(_raw(
        pubTypeList={"pubType": ["Journal Article", "Retracted Publication"]}))
    assert record.retracted and record.withdrawn


def test_a_retraction_notice_marks_the_record_even_without_the_publication_type():
    """Europe PMC reports the notice as a link from the article to the retraction."""
    record = literature.parse_record(_raw(commentCorrectionList={
        "commentCorrection": [{"type": "Retraction in", "id": "99000002"}]}))
    assert record.retracted and record.withdrawn


def test_an_expression_of_concern_withdraws_support_without_being_a_retraction():
    record = literature.parse_record(_raw(commentCorrectionList={
        "commentCorrection": [{"type": "Expression of concern in", "id": "99000003"}]}))
    assert record.withdrawn and record.expression_of_concern and not record.retracted


def test_an_erratum_is_a_correction_not_a_withdrawal():
    record = literature.parse_record(_raw(commentCorrectionList={
        "commentCorrection": [{"type": "Erratum in", "id": "99000004"}]}))
    assert record.corrected and not record.withdrawn


def test_a_record_missing_every_optional_block_still_parses():
    """Preprints and unindexed records arrive without MeSH, types or a journal."""
    record = literature.parse_record({"id": "99000009", "source": "PPR"})
    assert record.mesh == () and record.pub_types == () and record.journal == ""
    assert record.key == "99000009" or record.key == ""


# -------------------------------------------------------------------------- grading


@pytest.mark.parametrize("overrides, expected", [
    ({"pubTypeList": {"pubType": ["Randomized Controlled Trial"]}}, "clinical"),
    ({"meshHeadingList": {"meshHeading": [{"descriptorName": "Mice"}]}}, "in_vivo"),
    ({"meshHeadingList": {"meshHeading": [{"descriptorName": "Cells, Cultured"}]}},
     "in_vitro"),
    ({"meshHeadingList": {"meshHeading": []}}, "ungraded"),
])
def test_grades_from_indexing(overrides, expected):
    assert literature.grade(literature.parse_record(_raw(**overrides))) == expected


def test_clinical_outranks_an_animal_descriptor():
    """A trial indexed with Animals is still clinical evidence; the ladder is ordered."""
    record = literature.parse_record(_raw(
        pubTypeList={"pubType": ["Clinical Trial"]},
        meshHeadingList={"meshHeading": [{"descriptorName": "Mice"}]}))
    assert literature.grade(record) == "clinical"


def test_every_grade_is_one_of_the_declared_vocabulary():
    """The grade travels into L3 and the report, so the vocabulary is closed."""
    for overrides in ({}, {"pubTypeList": {"pubType": ["Review"]}},
                      {"meshHeadingList": {"meshHeading": [{"descriptorName": "Humans"}]}}):
        assert literature.grade(literature.parse_record(_raw(**overrides))) \
            in literature.GRADES


# ------------------------------------------------------------- proximity and matching


def test_cooccurrence_requires_one_sentence_not_one_record():
    """The whole point: both terms present somewhere is not evidence about the two."""
    apart = literature.parse_record(_raw(
        title="A study of invented things.",
        abstractText="Inventib was given to patients. Separately, aneuploidy was seen."))
    together = literature.parse_record(_raw(
        title="A study of invented things.",
        abstractText="Inventib selectively killed aneuploid cells."))
    left = [literature.term_pattern("inventib")]
    right = [literature.term_pattern("aneuploid*")]
    assert not literature.cooccur(apart, left, right)
    assert literature.cooccur(together, left, right)


def test_the_title_counts_as_a_sentence():
    record = literature.parse_record(_raw(
        title="Inventib is selective for aneuploid cells", abstractText=""))
    assert literature.cooccur(record, [literature.term_pattern("inventib")],
                              [literature.term_pattern("aneuploid*")])


def test_a_wildcard_term_means_the_same_in_the_query_and_in_the_text():
    pattern = literature.term_pattern("aneuploid*")
    assert pattern.search("aneuploidy") and pattern.search("aneuploid cells")
    assert not pattern.search("euploid")


def test_a_phrase_term_tolerates_a_line_break():
    assert literature.term_pattern("chromosomal instability").search(
        "chromosomal\n  instability")


def test_matches_returns_pattern_names_never_the_matched_text():
    """Abstracts are third-party content; outputs carry the name of what matched."""
    record = literature.parse_record(_raw(
        abstractText="Inventib showed selective antiproliferative activity."))
    found = literature.matches(record, (("selectivity", literature.term_pattern("selectiv")),
                                        ("absent", literature.term_pattern("nothing here"))))
    assert found == ("selectivity",)


# --------------------------------------------------------------------- search and cache


def _page(results, next_cursor=None) -> dict:
    page = {"hitCount": len(results), "resultList": {"result": results}}
    if next_cursor:
        page["nextCursorMark"] = next_cursor
    return page


def test_search_caches_the_response_and_does_not_refetch(monkeypatch, tmp_path):
    """A re-run must reproduce the run that produced a result, not re-query the index."""
    calls = []

    def fake(url, timeout=60):
        calls.append(url)
        return _page([_raw()])

    monkeypatch.setattr(literature, "_get_json", fake)
    config = {"reference_dir": str(tmp_path)}
    first, provenance = literature.search(config, 'TITLE_ABS:"inventib"')
    second, again = literature.search(config, 'TITLE_ABS:"inventib"')

    assert len(calls) == 1
    assert [r.pmid for r in first] == [r.pmid for r in second] == ["99000001"]
    assert provenance["response_sha256"] == again["response_sha256"]
    assert provenance["query"] == 'TITLE_ABS:"inventib"'
    assert provenance["hit_count"] == 1


def test_search_deduplicates_and_respects_max_records(monkeypatch, tmp_path):
    duplicate = [_raw(), _raw(), _raw(id="99000002", pmid="99000002")]
    monkeypatch.setattr(literature, "_get_json",
                        lambda url, timeout=60: _page(duplicate))
    records, _ = literature.search({"reference_dir": str(tmp_path)},
                                   'TITLE_ABS:"inventib"', max_records=5)
    assert [r.pmid for r in records] == ["99000001", "99000002"]


def test_search_stops_paging_when_the_cursor_repeats(monkeypatch, tmp_path):
    """A service that keeps returning the same cursor must not loop forever."""
    monkeypatch.setattr(literature, "_get_json",
                        lambda url, timeout=60: _page([_raw()], next_cursor="*"))
    records, _ = literature.search({"reference_dir": str(tmp_path)},
                                   'TITLE_ABS:"inventib"', max_records=50)
    assert len(records) == 1


def test_fetch_identifier_requires_a_typed_citation(tmp_path):
    """'21315436' alone could be a PMID or a truncated anything; never guess."""
    with pytest.raises(ValueError, match="PMID:"):
        literature.fetch_identifier({"reference_dir": str(tmp_path)}, "21315436")


def test_fetch_identifier_returns_nothing_for_an_unresolvable_citation(monkeypatch,
                                                                       tmp_path):
    monkeypatch.setattr(literature, "_get_json", lambda url, timeout=60: _page([]))
    record, provenance = literature.fetch_identifier(
        {"reference_dir": str(tmp_path)}, "DOI:10.5555/does.not.exist")
    assert record is None and provenance["hit_count"] == 0


def test_crossref_reports_an_unresolvable_doi_without_raising(monkeypatch, tmp_path):
    """Europe PMC has already established the record; Crossref is the second opinion."""
    def refuse(url, timeout=60):
        raise ValueError("no such work")

    monkeypatch.setattr(literature, "_get_json", refuse)
    assert literature.crossref({"reference_dir": str(tmp_path)},
                               "10.5555/missing")["resolved"] is False


def test_the_response_cache_lives_outside_the_repository(tmp_path):
    from src import refcache

    path = literature.cache_dir({"reference_dir": str(tmp_path)})
    assert refcache.REPO_ROOT not in path.parents and path != refcache.REPO_ROOT


# ------------------------------------------------------------------------ publishing


def test_reference_carries_bibliographic_metadata_and_no_abstract():
    """Abstracts are third-party content under mixed licences; they stay in the cache."""
    record = literature.parse_record(_raw(abstractText="Proprietary abstract text."))
    published = literature.reference(record)
    assert "Proprietary abstract text." not in json.dumps(published)
    assert "abstract" not in published
    assert published["pmid"] == "99000001" and published["grade"] == "in_vitro"
