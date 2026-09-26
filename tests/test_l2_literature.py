"""Tests for the Europe PMC / Crossref lookup behind L2 channel E.

Nothing here touches the network: every test either stubs the HTTP call or works on a
parsed record. The identifiers are invented -- PMIDs in the ``99xxxxxx`` range and DOIs
under the reserved ``10.5555`` test prefix -- so no test asserts anything about a real
paper, and a test can never be mistaken for a verified citation.
"""

from __future__ import annotations

import http.client
import json
import urllib.error

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


# Every identifier below is invented. The guard matches on shape, not on meaning, so a
# made-up term exercises it exactly as a real one would -- and a real HPO term or
# coordinate has no business in a committed file (CLAUDE.md, "Clinical phenotype is
# patient data"). The OMIM id is the disease's own public identifier, printed in
# CLAUDE.md, and carries nothing about this subject.
@pytest.mark.parametrize("query, what", [
    ('TITLE_ABS:"HP:9000001"', "an HPO term"),
    ('TITLE_ABS:"chr9:99990000"', "a coordinate"),
    ('"9-99999999"', "a bare coordinate"),
    ('TITLE_ABS:"c.9999_9999insA"', "an HGVS change"),
    ('"rs99999999"', "a dbSNP id"),
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
        literature.refuse_private('TITLE_ABS:"HP:9000002"')
    assert "HP:9000002" not in str(caught.value)


@pytest.mark.parametrize("query", [
    '(TITLE_ABS:chloroquine) AND (TITLE_ABS:aneuploid*)',
    '(TITLE_ABS:"BUB1B" OR TITLE_ABS:"BUBR1") AND (TITLE_ABS:"trisomy")',
    'EXT_ID:21315436 AND SRC:"MED"',
])
def test_allows_a_query_of_public_vocabulary(query):
    literature.refuse_private(query)


@pytest.mark.parametrize("doi", [
    # An AACR DOI ends in something shaped exactly like chromosome 13 at position 1174.
    "10.1158/0008-5472.CAN-13-1174",
    "10.1158/1078-0432.CCR-20-1234",
    # A suffix containing '.g.' reads as an HGVS genomic expression.
    "10.1002/j.g.12345",
])
def test_a_public_doi_is_not_mistaken_for_patient_data(doi):
    """A guard that cries wolf on a citation gets widened until it stops guarding."""
    literature.refuse_private(f"DOI:{doi}")


def test_a_coordinate_beside_a_doi_is_still_refused():
    """Masking DOIs must not blind the check to a coordinate elsewhere in the query."""
    with pytest.raises(ValueError, match="genomic coordinate"):
        literature.refuse_private('DOI:"10.1158/0008-5472.CAN-13-1174" OR "chr9:99999999"')


def test_search_checks_the_query_before_any_request(monkeypatch, tmp_path):
    """The guard must run before the socket opens, not after the response comes back."""
    def explode(*args, **kwargs):
        raise AssertionError("a request was made despite a refused query")

    monkeypatch.setattr(literature, "_get_json", explode)
    with pytest.raises(ValueError, match="refusing to send"):
        literature.search({"reference_dir": str(tmp_path)}, 'TITLE_ABS:"HP:9000001"')


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
    assert record.key == "PPR:99000009"
    assert record.citable, "Europe PMC's own id is a resolvable handle"


def test_a_record_with_no_identifier_at_all_is_not_citable():
    """A claim needs a citation a reader can resolve; a title is not one."""
    record = literature.parse_record({"title": "An invented study with no identifiers."})
    assert not record.citable and record.key == ""


def test_the_deduplication_key_never_falls_back_to_the_title():
    """Two unrelated records sharing a title must not collapse into one."""
    a = literature.parse_record({"id": "99000010", "source": "PPR", "title": "Same title."})
    b = literature.parse_record({"id": "99000011", "source": "PPR", "title": "Same title."})
    assert a.key != b.key


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


@pytest.mark.parametrize("error", [
    ValueError("no such work"),
    urllib.error.URLError("unreachable"),
    TimeoutError("read timed out"),          # socket.timeout: an OSError, NOT a URLError
    http.client.HTTPException("truncated"),
])
def test_crossref_survives_any_transport_failure(monkeypatch, tmp_path, error):
    """Europe PMC has already established the record; Crossref is the second opinion.

    A read timeout is the case that matters: ``socket.timeout`` is a ``TimeoutError``,
    which is an ``OSError`` but not a ``URLError``, so the narrower spelling let a slow
    Crossref kill the whole channel.
    """
    def refuse(url, timeout=60):
        raise error

    monkeypatch.setattr(literature, "_get_json", refuse)
    assert literature.crossref({"reference_dir": str(tmp_path)},
                               "10.5555/missing")["resolved"] is False


def test_a_crossref_outage_is_not_cached(monkeypatch, tmp_path):
    """Otherwise one bad minute makes a good DOI unresolvable for the life of the cache.

    Worse, a retraction the publisher deposits later would never be seen, because the
    lookup that would have found it never runs again.
    """
    config = {"reference_dir": str(tmp_path)}
    calls = []

    def flaky(url, timeout=60):
        calls.append(url)
        if len(calls) == 1:
            raise urllib.error.URLError("transient")
        return {"message": {"DOI": "10.5555/invented.1", "title": ["Recovered"],
                            "update-to": [{"type": "retraction", "DOI": "10.5555/r"}]}}

    monkeypatch.setattr(literature, "_get_json", flaky)
    assert literature.crossref(config, "10.5555/invented.1")["resolved"] is False
    second = literature.crossref(config, "10.5555/invented.1")
    assert len(calls) == 2, "the failure was cached instead of retried"
    assert second["resolved"] is True
    assert second["update_to"][0]["type"] == "retraction"


def test_a_successful_crossref_lookup_is_cached(monkeypatch, tmp_path):
    calls = []

    def once(url, timeout=60):
        calls.append(url)
        return {"message": {"DOI": "10.5555/invented.1", "title": ["A work"]}}

    monkeypatch.setattr(literature, "_get_json", once)
    config = {"reference_dir": str(tmp_path)}
    literature.crossref(config, "10.5555/invented.1")
    literature.crossref(config, "10.5555/invented.1")
    assert len(calls) == 1


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
