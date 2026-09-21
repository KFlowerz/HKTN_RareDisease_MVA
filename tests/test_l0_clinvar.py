"""Tests for the L0 ClinVar cross-reference.

Every input is invented: gene symbols (GENE1, GENE2), coordinates, alleles, ClinVar
variation ids, conditions and review statuses. No patient variant appears here, and no
real ClinVar record does either -- a real pathogenic allele in a committed test could be
read as a statement about the subject.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src import refcache
from src.l0_genomics import annotate, causal, clinvar, transcripts as tx

REPO_ROOT = Path(__file__).resolve().parent.parent

REGIONS = {"GENE1": ("1", 1000, 2000), "GENE2": ("2", 5000, 6000)}
CFG = {"annotator": {"tool": "snpeff", "database": "DBX", "java_heap": "8g",
                     "transcript_policy": "mane_select_tiered", "mane_release": "MANE.test"}}
PAD = 100

HEADER = [
    "##fileformat=VCFv4.1\n",
    "##fileDate=2026-09-05\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
]


def _row(contig="1", pos=1500, vid="111", ref="C", alt="T",
         clnsig="Pathogenic", review="criteria_provided,_single_submitter", extra="") -> str:
    info = (f"ALLELEID=1;CLNDN=Invented_condition;CLNREVSTAT={review};CLNSIG={clnsig};"
            f"CLNVC=single_nucleotide_variant;GENEINFO=GENE1:1;MC=SO:0001587|nonsense{extra}")
    return f"{contig}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n"


def _allele(vid="v1", gene="GENE1", contig="1", pos=1500, ref="C", alt="T") -> causal.Allele:
    return causal.Allele(variant_id=vid, gene=gene, contig=contig, pos=pos, ref=ref, alt=alt,
                         zygosity="het", filter="PASS")


def _call(hgvs_p="p.Arg550Ter", transcript="ENST00000000001.1") -> tx.VariantCall:
    return tx.VariantCall(gene="GENE1", transcript=transcript, annotation_pass="primary",
                          lof_tier="mane_select", mane_release="MANE.test",
                          transcript_policy="mane_select_tiered", consequence="stop_gained",
                          impact="HIGH", hgvs_c="c.1648C>T", hgvs_p=hgvs_p, variant_id="v1",
                          zygosity="het", effect_class="lof")


# ------------------------------------------------------- minimal representation


def test_minimal_leaves_an_snv_alone() -> None:
    assert clinvar.minimal(1500, "C", "T") == (1500, "C", "T")


def test_minimal_trims_a_padded_deletion() -> None:
    # The same deletion carrying one extra shared base must compare equal to the parsimonious
    # spelling. This is the difference between two VCFs' padding conventions.
    assert clinvar.minimal(100, "AGG", "AG") == clinvar.minimal(100, "AG", "A")


def test_minimal_does_not_realign_a_repeat() -> None:
    """The stated limit, asserted so it cannot quietly become a claim of more.

    ``100 AGG>AG`` and ``101 GG>G`` delete the same G, but only the first is left-aligned.
    Minimal representation removes padding; it does not move an event through a repeat, so
    a source that left-aligned differently is reported as absent, not as a match. That is
    the caveat carried in the artifact.
    """
    assert clinvar.minimal(100, "AGG", "AG") != clinvar.minimal(101, "GG", "G")
    assert any("left-align" in caveat for caveat in clinvar.CAVEATS)


def test_minimal_trims_trailing_before_leading() -> None:
    assert clinvar.minimal(10, "CTT", "CAT") == (11, "T", "A")


def test_minimal_never_empties_an_allele() -> None:
    pos, ref, alt = clinvar.minimal(10, "AAAA", "A")
    assert ref and alt and (pos, ref, alt) == (10, "AAAA", "A")


def test_minimal_is_case_insensitive() -> None:
    assert clinvar.minimal(10, "c", "t") == clinvar.minimal(10, "C", "T")


# ------------------------------------------------------------- CLNSIG buckets


@pytest.mark.parametrize("clnsig,expected", [
    ("Pathogenic", "pathogenic"),
    ("Likely_pathogenic", "likely_pathogenic"),
    ("Pathogenic/Likely_pathogenic", "pathogenic"),
    ("Pathogenic|risk_factor", "pathogenic"),
    ("Pathogenic_low_penetrance", "pathogenic"),
    ("Benign", "benign"),
    ("Benign/Likely_benign", "benign"),
    ("Likely_benign", "likely_benign"),
    ("Uncertain_significance", "uncertain"),
    ("Conflicting_classifications_of_pathogenicity", "conflicting"),
    ("Conflicting_interpretations_of_pathogenicity", "conflicting"),
    ("drug_response", "other"),
    ("", "other"),
])
def test_significance_class(clnsig: str, expected: str) -> None:
    assert clinvar.significance_class(clnsig) == expected


def test_conflicting_is_not_reported_as_pathogenic() -> None:
    """A conflicting record must never fall into the pathogenic bucket.

    It is the one significance that names pathogenicity while meaning "submitters
    disagree", and counting it as support would manufacture evidence for a G1 call.
    """
    assert clinvar.significance_class(
        "Conflicting_classifications_of_pathogenicity") not in clinvar.PATHOGENIC_CLASSES


# --------------------------------------------------------------- review status


def test_review_stars_cover_the_documented_statuses() -> None:
    assert clinvar.REVIEW_STARS["practice_guideline"] == 4
    assert clinvar.REVIEW_STARS["reviewed_by_expert_panel"] == 3
    assert clinvar.REVIEW_STARS["criteria_provided,_multiple_submitters,_no_conflicts"] == 2
    assert clinvar.REVIEW_STARS["criteria_provided,_single_submitter"] == 1
    assert clinvar.REVIEW_STARS["no_assertion_criteria_provided"] == 0


def test_unknown_review_status_is_none_not_zero() -> None:
    """An unrecognised status must not be scored. Zero stars is a claim of its own."""
    _, records = clinvar.parse_records(HEADER + [_row(review="some_future_status")],
                                       REGIONS, pad=PAD)
    assert records[0].stars is None
    assert records[0].review == "some_future_status"


# ------------------------------------------------------------------- parsing


def test_parse_records_keeps_the_file_date_and_region_records() -> None:
    date, records = clinvar.parse_records(
        HEADER + [_row(pos=1500), _row(contig="2", pos=5500, vid="222")], REGIONS, pad=PAD)
    assert date == "2026-09-05"
    assert {r.variation_id for r in records} == {"111", "222"}


def test_parse_records_drops_records_outside_every_region() -> None:
    _, records = clinvar.parse_records(
        HEADER + [_row(pos=9_000_000), _row(contig="7", pos=1500, vid="333")],
        REGIONS, pad=PAD)
    assert records == []


def test_parse_records_keeps_the_padding_window() -> None:
    _, records = clinvar.parse_records(HEADER + [_row(pos=2000 + PAD)], REGIONS, pad=PAD)
    assert len(records) == 1


def test_parse_records_splits_multiallelic_records() -> None:
    _, records = clinvar.parse_records(HEADER + [_row(alt="T,G")], REGIONS, pad=PAD)
    assert sorted(r.alt for r in records) == ["G", "T"]


def test_parse_records_does_not_merge_two_regions_on_one_contig() -> None:
    """Each gene gets its own interval, not a min/max span across the contig.

    ``l0_panel_extension`` can put two panel genes on the same chromosome. Merging their
    spans would sweep in every record between them -- on a big chromosome, tens of
    thousands of records, annotated and matched for nothing.
    """
    regions = {"GENE1": ("1", 1000, 2000), "GENE3": ("1", 9_000_000, 9_001_000)}
    _, records = clinvar.parse_records(
        HEADER + [_row(pos=1500), _row(pos=5_000_000, vid="444"), _row(pos=9_000_500, vid="555")],
        regions, pad=PAD)
    assert {r.variation_id for r in records} == {"111", "555"}


def test_multiallelic_alleles_keep_separate_protein_annotations() -> None:
    """ClinVar gives every ALT of a multi-allelic record the same variation id.

    Keyed by that id, one allele's amino-acid change would be attached to the other. The
    panel's current records are biallelic throughout, so this would go unnoticed.
    """
    _, records = clinvar.parse_records(HEADER + [_row(alt="T,A")], REGIONS, pad=PAD)
    assert len({r.variation_id for r in records}) == 1
    assert len({r.key for r in records}) == 2

    by_alt = {r.alt: r for r in records}
    protein = {by_alt["T"].key: (("ENST00000000001.1", "p.Arg550Ter"),),
               by_alt["A"].key: (("ENST00000000001.1", "p.Arg550Gln"),)}
    match = clinvar.match_allele(_allele(alt="T"), _call("p.Arg550Ter"), records, protein)
    residue = match["same_protein_position_pathogenic"]
    assert [r["hgvs_p"] for r in residue] == ["p.Arg550Gln"]


def test_parse_records_reads_the_legacy_frequency_fields() -> None:
    _, records = clinvar.parse_records(HEADER + [_row(extra=";AF_TGP=0.02;AF_ESP=0.01")],
                                       REGIONS, pad=PAD)
    assert records[0].frequencies == {"af_tgp": 0.02, "af_esp": 0.01}


def test_parse_records_keeps_the_molecular_consequence_term() -> None:
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    assert records[0].molecular_consequence == "nonsense"


# ------------------------------------------------------------------- matching


def test_exact_match_reports_the_record() -> None:
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    match = clinvar.match_allele(_allele(), None, records, {})
    assert match["match"] == "exact"
    assert match["records"][0]["significance"] == "pathogenic"
    assert match["records"][0]["stars"] == 1


def test_a_different_allele_at_the_same_position_is_not_a_match() -> None:
    """Context, not identity: a different base at the same site is a different variant."""
    _, records = clinvar.parse_records(HEADER + [_row(alt="G")], REGIONS, pad=PAD)
    match = clinvar.match_allele(_allele(alt="T"), None, records, {})
    assert match["match"] == "same_position"
    assert match["records"] == []
    assert len(match["same_position_records"]) == 1


def test_absence_is_reported_as_none() -> None:
    _, records = clinvar.parse_records(HEADER + [_row(pos=1600)], REGIONS, pad=PAD)
    match = clinvar.match_allele(_allele(pos=1500), None, records, {})
    assert match["match"] == "none"
    assert match["records"] == [] and match["same_position_records"] == []


def test_padded_indel_representations_still_match() -> None:
    row = _row(pos=1500, ref="AG", alt="A")
    _, records = clinvar.parse_records(HEADER + [row], REGIONS, pad=PAD)
    match = clinvar.match_allele(_allele(pos=1500, ref="AGG", alt="AG"), None, records, {})
    assert match["match"] == "exact"


def test_same_protein_position_is_reported_with_same_change_flagged() -> None:
    rows = [_row(vid="222", alt="A"), _row(vid="333", alt="G")]
    _, records = clinvar.parse_records(HEADER + rows, REGIONS, pad=PAD)
    protein = {records[0].key: (("ENST00000000001.1", "p.Arg550Ter"),),
               records[1].key: (("ENST00000000001.1", "p.Arg550Gln"),)}
    match = clinvar.match_allele(_allele(alt="T"), _call("p.Arg550Ter"), records, protein)
    residue = match["same_protein_position_pathogenic"]
    assert [r["variation_id"] for r in residue] == ["222", "333"]
    assert residue[0]["same_change"] is True and residue[1]["same_change"] is False


def test_the_exact_match_is_not_repeated_as_same_residue_evidence() -> None:
    """Otherwise one record reads as two -- the match, plus support for its own residue."""
    _, records = clinvar.parse_records(HEADER + [_row(vid="111", alt="T")], REGIONS, pad=PAD)
    protein = {records[0].key: (("ENST00000000001.1", "p.Arg550Ter"),)}
    match = clinvar.match_allele(_allele(alt="T"), _call("p.Arg550Ter"), records, protein)
    assert match["match"] == "exact"
    assert match["same_protein_position_pathogenic"] == []


def test_same_protein_position_ignores_another_transcript() -> None:
    """Residue numbers are only comparable on the transcript the call was made on."""
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    protein = {records[0].key: (("ENST00000000002.1", "p.Arg550Ter"),)}
    match = clinvar.match_allele(_allele(), _call(), records, protein)
    assert match["same_protein_position_pathogenic"] == []


def test_same_protein_position_needs_a_residue_number() -> None:
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    protein = {records[0].key: (("ENST00000000001.1", "p.Arg550Ter"),)}
    match = clinvar.match_allele(_allele(), _call(hgvs_p=""), records, protein)
    assert match["same_protein_position_pathogenic"] == []


@pytest.mark.parametrize("hgvs_p,expected", [
    ("p.Arg550Ter", 550), ("p.Leu1012Pro", 1012), ("p.Ter1051ext*?", 1051),
    ("", None), ("p.=", None), ("c.1648C>T", None),
])
def test_protein_position(hgvs_p: str, expected) -> None:
    assert clinvar.protein_position(hgvs_p) == expected


# ------------------------------------------------- attachment to the G1 artifact


def test_attach_clinvar_puts_the_summary_on_the_configuration() -> None:
    genes = {"GENE1": {"configurations": [{"variant_ids": ["v1", "v2"]}]}}
    crossref = {"status": "ok", "alleles": [
        {"variant_id": "v1", "match": "exact", "records": [
            {"significance": "pathogenic", "stars": 2}],
         "same_position_records": [], "same_protein_position_pathogenic": []},
        {"variant_id": "v2", "match": "none", "records": [],
         "same_position_records": [], "same_protein_position_pathogenic": [{"x": 1}]},
    ]}
    causal.attach_clinvar(genes, crossref)
    attached = genes["GENE1"]["configurations"][0]["clinvar"]
    assert [c["variant_id"] for c in attached] == ["v1", "v2"]
    assert attached[0]["significance"] == "pathogenic" and attached[0]["stars"] == 2
    assert attached[1]["significance"] is None
    assert attached[1]["n_same_protein_position_pathogenic"] == 1


def test_attach_clinvar_is_a_no_op_when_the_crossref_did_not_run() -> None:
    genes = {"GENE1": {"configurations": [{"variant_ids": ["v1"]}]}}
    causal.attach_clinvar(genes, {"status": "disabled"})
    assert "clinvar" not in genes["GENE1"]["configurations"][0]


# ----------------------------------------------------------------- guardrails


def test_enabled_defaults_on_and_can_be_turned_off() -> None:
    assert clinvar.enabled({}) is True
    assert clinvar.enabled({"clinvar": {"enabled": False}}) is False


def test_the_reference_cache_may_not_live_inside_the_repo(tmp_path, monkeypatch) -> None:
    """The cache holds a ClinVar ``.vcf.gz``, so it must sit outside the working tree.

    ``test_no_patient_data_files_present`` forbids a genomic file anywhere in the repo.
    That rule stays absolute -- a public release is not a reason to carve an exception
    into it, because the next genomic file to land there may not be public. So the cache
    refuses a repo-internal location instead.
    """
    monkeypatch.delenv("MVA_REF_ROOT", raising=False)
    with pytest.raises(ValueError, match="inside the repository"):
        refcache.reference_dir({"reference_dir": str(REPO_ROOT / "reference")})

    monkeypatch.setenv("MVA_REF_ROOT", str(tmp_path / "cache"))
    assert refcache.reference_dir({}).is_dir()


def test_the_default_cache_is_outside_the_repo(monkeypatch) -> None:
    monkeypatch.delenv("MVA_REF_ROOT", raising=False)
    resolved = Path(refcache.DEFAULT_DIR).expanduser()
    assert resolved.is_absolute()
    assert REPO_ROOT not in resolved.parents and resolved != REPO_ROOT


def test_the_configured_cache_is_outside_the_repo() -> None:
    """The committed config must not point the cache back into the working tree."""
    import yaml

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    resolved = Path(config["reference_dir"]).expanduser()
    assert REPO_ROOT not in resolved.parents and resolved != REPO_ROOT


def test_the_source_never_builds_a_per_variant_query() -> None:
    """The whole release is downloaded and matched locally.

    A per-variant lookup would send this child's coordinates to NCBI, which COMPLIANCE.md
    prohibits. The rule is easy to break later by "just checking one variant", so it is
    asserted here rather than left to the docstring.
    """
    source = Path(clinvar.__file__).read_text(encoding="utf-8")
    for forbidden in ("eutils", "esearch", "efetch", "/api/", "?term=", "requests."):
        assert forbidden not in source
    # A whole VCF release with no query attached. Asserted on the shape rather than on a
    # literal filename: the URL is now pinned to a DATED weekly release
    # (clinvar_YYYYMMDD.vcf.gz) instead of the rolling clinvar.vcf.gz, and the property
    # that matters here is "a whole file, not a lookup", which both forms satisfy.
    for url in (clinvar.CLINVAR_URL, clinvar.CLINVAR_ARCHIVE_URL):
        assert url.endswith(".vcf.gz")
        assert "?" not in url and "=" not in url


# ------------------------------------------------- reading a cached release


def _gzipped(tmp_path, lines, name="release.txt.gz"):
    """Write invented VCF text to a gzip file.

    Deliberately not named ``.vcf.gz``: no genomic file may exist in the repo or its test
    run, and :func:`clinvar.load_records` cares about gzip, not about the extension.
    """
    import gzip

    path = tmp_path / name
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.writelines(lines)
    return path


def test_load_records_reads_a_gzipped_release(tmp_path) -> None:
    path = _gzipped(tmp_path, HEADER + [_row()])
    date, records = clinvar.load_records(path, REGIONS, pad=PAD)
    assert date == "2026-09-05" and [r.variation_id for r in records] == ["111"]


def test_load_records_refuses_an_empty_result(tmp_path) -> None:
    """Zero records in six well-studied genes means the contigs do not line up.

    Returning an empty set would read as "none of these variants are known", which is the
    opposite of what a `chr`-prefix mismatch actually means -- the same failure
    ``read_alleles`` guards against on the subject's own VCF.
    """
    path = _gzipped(tmp_path, HEADER + [_row(contig="chr1")])
    with pytest.raises(ValueError, match="contig-naming"):
        clinvar.load_records(path, REGIONS, pad=PAD)


# ---------------------------------------------------------- region assignment


def test_gene_of_assigns_a_record_to_its_region() -> None:
    _, records = clinvar.parse_records(HEADER + [_row(contig="2", pos=5500)], REGIONS, pad=PAD)
    assert clinvar.gene_of(records[0], REGIONS, pad=PAD) == "GENE2"


def test_gene_of_returns_empty_outside_every_region() -> None:
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    assert clinvar.gene_of(records[0], {"GENE2": ("2", 5000, 6000)}, pad=PAD) == ""


def test_gene_of_is_deterministic_where_regions_overlap() -> None:
    """Two panel genes can overlap. Symbol order decides, not dict order."""
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    overlapping = {"GENE9": ("1", 1000, 2000), "GENE1": ("1", 1400, 1600)}
    assert clinvar.gene_of(records[0], overlapping, pad=PAD) == "GENE1"
    assert clinvar.gene_of(records[0], dict(reversed(list(overlapping.items()))),
                           pad=PAD) == "GENE1"


# ------------------------------------------ annotating the pathogenic records


def _fake_annotate(monkeypatch, by_id: dict, seen: list = None):
    """Replace the snpEff call with a lookup, recording what was submitted."""
    def fake(config, records, flags):
        if seen is not None:
            seen.extend(records)
        return {vid: (tuple(by_id.get(vid, ())), frozenset()) for vid, *_ in records}

    monkeypatch.setattr(clinvar.annotate, "annotate", fake)


def _tx_ann(transcript="ENST00000000001.1", hgvs_p="p.Arg550Ter", feature="transcript"):
    return annotate.Annotation(allele="T", effects=("missense_variant",), impact="MODERATE",
                               gene="GENE1", feature_type=feature, transcript=transcript,
                               biotype="protein_coding", rank="3/10", hgvs_c="c.10C>T",
                               hgvs_p=hgvs_p)


def test_annotate_pathogenic_submits_only_pathogenic_records(monkeypatch) -> None:
    """A same-residue benign or uncertain record says nothing a person would act on."""
    rows = [_row(vid="111", clnsig="Pathogenic"),
            _row(vid="222", alt="A", clnsig="Likely_pathogenic"),
            _row(vid="333", alt="G", clnsig="Benign"),
            _row(vid="444", alt="C", clnsig="Uncertain_significance")]
    _, records = clinvar.parse_records(HEADER + rows, REGIONS, pad=PAD)
    seen = []
    _fake_annotate(monkeypatch, {}, seen)
    clinvar.annotate_pathogenic(CFG, records)
    assert len(seen) == 2


def test_annotate_pathogenic_keys_by_allele_and_keeps_every_transcript(monkeypatch) -> None:
    _, records = clinvar.parse_records(HEADER + [_row(alt="T,A")], REGIONS, pad=PAD)
    by_alt = {r.alt: r for r in records}
    _fake_annotate(monkeypatch, {
        "cv0": (_tx_ann(hgvs_p="p.Arg550Ter"), _tx_ann("ENST00000000002.1", "p.Arg9Gln")),
        "cv1": (_tx_ann(hgvs_p="p.Arg550Gln"),),
    })
    protein = clinvar.annotate_pathogenic(CFG, records)
    assert set(protein) == {by_alt["T"].key, by_alt["A"].key}
    assert len(protein[by_alt["T"].key]) == 2


def test_annotate_pathogenic_skips_records_with_no_residue(monkeypatch) -> None:
    """A splice or UTR record has no amino-acid position to compare against."""
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    _fake_annotate(monkeypatch, {"cv0": (_tx_ann(hgvs_p=""),)})
    assert clinvar.annotate_pathogenic(CFG, records) == {}


def test_annotate_pathogenic_ignores_non_transcript_features(monkeypatch) -> None:
    _, records = clinvar.parse_records(HEADER + [_row()], REGIONS, pad=PAD)
    _fake_annotate(monkeypatch, {"cv0": (_tx_ann(feature="gene_variant"),)})
    assert clinvar.annotate_pathogenic(CFG, records) == {}


def test_annotate_pathogenic_runs_nothing_without_pathogenic_records(monkeypatch) -> None:
    """No snpEff invocation at all -- it costs a database load."""
    _, records = clinvar.parse_records(HEADER + [_row(clnsig="Benign")], REGIONS, pad=PAD)

    def explode(*args, **kwargs):
        raise AssertionError("snpEff must not run when there is nothing to annotate")

    monkeypatch.setattr(clinvar.annotate, "annotate", explode)
    assert clinvar.annotate_pathogenic(CFG, records) == {}


# --------------------------------------------------------- the download cache


def _fake_urlopen(monkeypatch, payload: bytes, calls: list):
    class Response:
        def __init__(self):
            self._data = payload

        def read(self, size):
            chunk, self._data = self._data[:size], self._data[size:]
            return chunk

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def fake(request, timeout=None):
        calls.append(request.full_url)
        return Response()

    monkeypatch.setattr(refcache.urllib.request, "urlopen", fake)


def test_fetch_downloads_once_and_reuses_the_cache(tmp_path, monkeypatch) -> None:
    calls = []
    _fake_urlopen(monkeypatch, b"payload", calls)
    dest = tmp_path / "release.txt"

    first = refcache.fetch("https://example.org/release.txt", dest)
    second = refcache.fetch("https://example.org/release.txt", dest)

    assert calls == ["https://example.org/release.txt"]
    assert dest.read_bytes() == b"payload"
    assert first["bytes"] == second["bytes"] == 7
    assert first["sha256"] == second["sha256"] == refcache.sha256(dest)


def test_fetch_leaves_no_partial_file_behind(tmp_path, monkeypatch) -> None:
    """An interrupted download must not be mistaken for a complete one next run."""
    calls = []
    _fake_urlopen(monkeypatch, b"payload", calls)
    dest = tmp_path / "release.txt"
    refcache.fetch("https://example.org/release.txt", dest)
    assert list(tmp_path.iterdir()) == [dest]


def test_fetch_identifies_itself(tmp_path, monkeypatch) -> None:
    """A public archive is entitled to know who is pulling a whole release."""
    assert "User-Agent" in refcache.USER_AGENT
    assert refcache.USER_AGENT["User-Agent"].startswith("mva-track2")


# ------------------------------------- what survives a cross-reference failure


def _stub_causal(monkeypatch, crossref):
    """Run ``causal.call`` without a VCF, snpEff or a network."""
    alleles = [
        causal.Allele(variant_id="v1", gene="GENE1", contig="1", pos=1500, ref="C", alt="T",
                      zygosity="het", filter="PASS"),
        causal.Allele(variant_id="v2", gene="GENE1", contig="1", pos=1700, ref="G", alt="A",
                      zygosity="het", filter="PASS"),
    ]
    lof = annotate.Annotation(allele="T", effects=("stop_gained",), impact="HIGH", gene="GENE1",
                              feature_type="transcript", transcript="ENST00000000001.1",
                              biotype="protein_coding", rank="3/10", hgvs_c="c.10C>T",
                              hgvs_p="p.Arg550Ter")
    missense = annotate.Annotation(allele="A", effects=("missense_variant",), impact="MODERATE",
                                   gene="GENE1", feature_type="transcript",
                                   transcript="ENST00000000001.1", biotype="protein_coding",
                                   rank="4/10", hgvs_c="c.20G>A", hgvs_p="p.Gly7Asp")

    monkeypatch.setattr(causal.annotate, "gene_regions",
                        lambda config, genes: {"GENE1": ("1", 1000, 2000)})
    monkeypatch.setattr(causal, "read_alleles", lambda path, regions, pad=0: alleles)
    monkeypatch.setattr(causal.annotate, "annotate", lambda config, records, flags: {
        "v1": ((lof,), frozenset({"GENE1"})), "v2": ((missense,), frozenset({"GENE1"}))})
    monkeypatch.setattr(causal.annotate, "snpeff_version", lambda: "SnpEff\ttest")
    monkeypatch.setattr(causal._clinvar, "cross_reference", crossref)
    # The gnomAD lookup has its own integration tests (tests/test_l0_gnomad.py).
    monkeypatch.setattr(causal._gnomad, "enabled", lambda config: False)


def test_annotation_survives_a_failed_cross_reference(tmp_path, monkeypatch) -> None:
    """The cross-reference downloads a release; the annotation before it must not be lost.

    Same rule as run.py writing the burden before this half starts: a network failure may
    not discard work that already succeeded.
    """
    def boom(*args, **kwargs):
        raise OSError("network unreachable")

    _stub_causal(monkeypatch, boom)
    with pytest.raises(OSError):
        causal.call("unused.vcf.gz", dict(CFG, seed=42, clinvar={"enabled": True}),
                    tmp_path, panel=("GENE1",))

    assert (tmp_path / "variant_calls.json").exists()
    assert not (tmp_path / "causal_gene_call.json").exists(), (
        "the G1 artifact must not be written without the evidence it claims to weigh"
    )


def test_the_g1_artifact_carries_the_clinvar_status(tmp_path, monkeypatch) -> None:
    import json

    crossref = {
        "status": "ok", "release": {"file_date": "2026-09-13"},
        "alleles": [
            {"variant_id": "v1", "match": "exact",
             "records": [{"significance": "pathogenic", "stars": 2}],
             "same_position_records": [], "same_protein_position_pathogenic": []},
            {"variant_id": "v2", "match": "none", "records": [],
             "same_position_records": [], "same_protein_position_pathogenic": []},
        ],
    }
    _stub_causal(monkeypatch, lambda *a, **k: crossref)
    summary = causal.call("unused.vcf.gz", dict(CFG, seed=42, clinvar={"enabled": True}),
                          tmp_path, panel=("GENE1",))

    assert summary["clinvar"] == "ok"
    written = json.loads((tmp_path / "causal_gene_call.json").read_text(encoding="utf-8"))
    assert written["clinvar"]["release"] == "2026-09-13"
    attached = written["per_gene"]["GENE1"]["configurations"][0]["clinvar"]
    assert [c["match"] for c in attached] == ["exact", "none"]
    assert (tmp_path / "clinvar_crossref.json").exists()


def test_disabling_the_cross_reference_says_so_in_the_artifact(tmp_path, monkeypatch) -> None:
    """"Disabled" and "nothing known" must never look alike in the G1 artifact."""
    import json

    def explode(*args, **kwargs):
        raise AssertionError("cross_reference must not run when disabled")

    _stub_causal(monkeypatch, explode)
    causal.call("unused.vcf.gz", dict(CFG, seed=42, clinvar={"enabled": False}),
                tmp_path, panel=("GENE1",))

    written = json.loads((tmp_path / "causal_gene_call.json").read_text(encoding="utf-8"))
    assert written["clinvar"]["status"] == "disabled"
    assert "clinvar" not in written["per_gene"]["GENE1"]["configurations"][0]
    assert not (tmp_path / "clinvar_crossref.json").exists()


def test_the_configured_release_url_is_a_whole_file(tmp_path) -> None:
    """A query string in the configured URL would be a per-variant lookup by another name."""
    import yaml

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    url = config["clinvar"]["url"]
    assert "?" not in url and url.endswith(".vcf.gz")
    assert config["clinvar"]["enabled"] is True


# ---------------------------------------------- the assembled artifact itself


def test_summarize_survives_an_allele_with_no_match() -> None:
    """The compact form must not assume a record exists -- most alleles have none."""
    summary = clinvar.summarize({"variant_id": "v9", "match": "none", "records": [],
                                 "same_position_records": [{"x": 1}],
                                 "same_protein_position_pathogenic": []})
    assert summary["significance"] is None and summary["stars"] is None
    assert summary["n_same_position_records"] == 1


def test_cross_reference_assembles_the_artifact(tmp_path, monkeypatch) -> None:
    """Pins the artifact's contract: provenance, per-gene denominators, caveats, seed.

    L5 and the report read this file, and a silently renamed key would surface as a
    missing figure rather than an error. The network fetch and the snpEff pass are stubbed
    -- both are covered separately, and the real ones run in L0.
    """
    _, records = clinvar.parse_records(
        HEADER + [_row(vid="111"), _row(vid="222", alt="A", clnsig="Benign")],
        REGIONS, pad=PAD)

    monkeypatch.setattr(refcache, "reference_dir", lambda config: tmp_path)
    monkeypatch.setattr(refcache, "fetch", lambda url, dest: {
        "url": url, "file": dest.name, "bytes": 7, "sha256": "abc", "retrieved": "2026-09-16"})
    monkeypatch.setattr(clinvar, "load_records",
                        lambda path, regions, pad: ("2026-09-13", records))
    monkeypatch.setattr(clinvar, "annotate_pathogenic", lambda config, recs: {})
    monkeypatch.setattr(clinvar.annotate, "snpeff_version", lambda: "SnpEff\ttest")

    config = dict(CFG, seed=42, clinvar={"enabled": True, "url": "https://x/clinvar.vcf.gz"})
    result = clinvar.cross_reference(config, [_allele()], {"v1": _call()},
                                     {"GENE1": ("1", 1000, 2000)}, pad=PAD)

    assert result["status"] == "ok" and result["seed"] == 42
    assert result["release"]["file_date"] == "2026-09-13"
    assert result["release"]["licence"] == clinvar.LICENCE
    assert result["gene_context"]["GENE1"] == {
        "n_records": 2, "n_pathogenic": 1, "n_annotated_on_call_transcript": 0}
    assert [m["variant_id"] for m in result["alleles"]] == ["v1"]
    assert result["alleles"][0]["match"] == "exact"
    assert result["caveats"] == list(clinvar.CAVEATS)
    assert "landrum2018 doi:10.1093/nar/gkx1153" in result["provenance"]["citations"]


def test_cross_reference_counts_what_the_vocabulary_did_not_map(tmp_path, monkeypatch) -> None:
    """Format drift is the one failure the suite cannot catch, so the artifact reports it.

    Every row in these tests is invented, so a renamed ClinVar INFO key would pass the
    whole suite while the real run produced records with no significance and no stars.
    Counting the unmapped values makes that visible to whoever reads the artifact.
    """
    rows = [_row(vid="111"),
            _row(vid="222", alt="A", review="renamed_in_a_future_release"),
            _row(vid="333", alt="G", review="renamed_in_a_future_release")]
    _, records = clinvar.parse_records(HEADER + rows, REGIONS, pad=PAD)

    monkeypatch.setattr(refcache, "reference_dir", lambda config: tmp_path)
    monkeypatch.setattr(refcache, "fetch", lambda url, dest: {
        "url": url, "file": dest.name, "bytes": 0, "sha256": "", "retrieved": "2026-09-16"})
    monkeypatch.setattr(clinvar, "load_records",
                        lambda path, regions, pad: ("2026-09-13", records))
    monkeypatch.setattr(clinvar, "annotate_pathogenic", lambda config, recs: {})
    monkeypatch.setattr(clinvar.annotate, "snpeff_version", lambda: "SnpEff\ttest")

    result = clinvar.cross_reference(dict(CFG, seed=42, clinvar={}), [], {},
                                     {"GENE1": ("1", 1000, 2000)}, pad=PAD)
    check = result["format_check"]
    assert check["n_records_parsed"] == 3
    assert check["unmapped_review_statuses"] == {"renamed_in_a_future_release": 2}


def test_cross_reference_reads_the_release_it_downloaded(tmp_path, monkeypatch) -> None:
    """One path, built once: the file fetched is the file parsed."""
    fetched, parsed = [], []
    monkeypatch.setattr(refcache, "reference_dir", lambda config: tmp_path)
    monkeypatch.setattr(refcache, "fetch", lambda url, dest: fetched.append(dest) or {
        "url": url, "file": dest.name, "bytes": 0, "sha256": "", "retrieved": "2026-09-16"})
    monkeypatch.setattr(clinvar, "load_records", lambda path, regions, pad: (
        parsed.append(path) or ("2026-09-13", list(clinvar.parse_records(
            HEADER + [_row()], REGIONS, pad=PAD)[1]))))
    monkeypatch.setattr(clinvar, "annotate_pathogenic", lambda config, recs: {})
    monkeypatch.setattr(clinvar.annotate, "snpeff_version", lambda: "SnpEff\ttest")

    clinvar.cross_reference(dict(CFG, seed=42, clinvar={"url": "https://x/clinvar.vcf.gz"}),
                            [], {}, {"GENE1": ("1", 1000, 2000)}, pad=PAD)
    assert fetched == parsed == [tmp_path / "clinvar.vcf.gz"]
