"""Tests for L2 channel D: phenotype similarity, inputs, and the channel.

Every term, disease and drug here is invented. HPO ids use the ``HP:9xxxxxx`` range, which
the real ontology does not use, so no test can encode a real phenotype -- let alone the
subject's (CLAUDE.md constraint 7). Mondo and ChEMBL ids are likewise made up.
"""

from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path

import pytest

from src.l2_channels import channel_d_phenotype as channel
from src.l2_channels import enrichment, phenosim, phenotype

ROOT, ABN, MOI = "HP:0000001", "HP:0000118", "HP:0000005"
A, A1, A2, B, B1, C = ("HP:9000001", "HP:9000011", "HP:9000012", "HP:9000002",
                       "HP:9000021", "HP:9000003")

PARENTS = {ROOT: [], ABN: [ROOT], MOI: [ROOT], A: [ABN], A1: [A], A2: [A], B: [ABN],
           B1: [B], C: [ABN]}

OBO = f"""format-version: 1.2
data-version: hp/releases/2099-01-01

[Term]
id: {ROOT}
name: All

[Term]
id: {ABN}
name: Invented abnormality root
is_a: {ROOT} ! All

[Term]
id: {MOI}
name: Invented mode of inheritance
is_a: {ROOT} ! All

[Term]
id: {A}
name: Invented A
alt_id: HP:9900001
is_a: {ABN}

[Term]
id: {A1}
name: Invented A1
is_a: {A}

[Term]
id: {A2}
name: Invented A2
is_a: {A}

[Term]
id: {B}
name: Invented B
is_a: {ABN}

[Term]
id: {B1}
name: Invented B1
is_a: {B}

[Term]
id: {C}
name: Invented C
is_a: {ABN}

[Term]
id: HP:9000099
name: Invented obsolete
is_obsolete: true
replaced_by: {B1}

[Typedef]
id: part_of
""".splitlines(keepends=True)


def _onto():
    return phenosim.Ontology.build(PARENTS)


def _corpus(onto, annotations=None):
    annotations = annotations or {
        "MONDO:9000001": {A1, B1},
        "MONDO:9000002": {A2},
        "MONDO:9000003": {B1, C},
        "MONDO:9000004": {C},
    }
    return phenosim.Corpus.build(onto, {d: {onto.index[t] for t in ts}
                                        for d, ts in annotations.items()})


# ------------------------------------------------------------------ ontology


def test_closures_include_self_and_follow_is_a() -> None:
    onto = _onto()
    anc = {onto.ids[i] for i in onto.ancestors[onto.index[A1]]}
    assert anc == {A1, A, ABN, ROOT}
    desc = {onto.ids[i] for i in onto.descendants[onto.index[A]]}
    assert desc == {A, A1, A2}


def test_is_under() -> None:
    onto = _onto()
    assert onto.is_under(onto.index[B1], onto.index[ABN])
    assert not onto.is_under(onto.index[MOI], onto.index[ABN])


def test_information_content_rises_with_rarity() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    ic = {t: corpus.ic[onto.index[t]] for t in PARENTS}
    assert ic[ABN] == pytest.approx(0.0)          # every disease is under it
    assert ic[A1] > ic[A] > ic[ABN]
    assert ic[A2] == pytest.approx(ic[A1])        # one disease each


def test_unannotated_terms_get_a_finite_information_content() -> None:
    onto = _onto()
    corpus = _corpus(onto, {"MONDO:9000001": {C}})
    assert corpus.ic[onto.index[A1]] == pytest.approx(0.0)  # floored at one of one disease


def test_term_similarity_is_the_most_informative_common_ancestor() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    sim = phenosim.term_similarity(onto, corpus, onto.index[A1])
    assert sim[onto.index[A2]] == pytest.approx(corpus.ic[onto.index[A]])
    assert sim[onto.index[A1]] == pytest.approx(corpus.ic[onto.index[A1]])
    assert sim[onto.index[B1]] == pytest.approx(0.0)


def test_an_exact_match_scores_one() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    query = [onto.index[A1], onto.index[B1]]
    scores, _ = phenosim.score_diseases(onto, corpus, query, seed=1, permutations=20)
    best = {s.disease: s for s in scores}["MONDO:9000001"]
    assert best.score == pytest.approx(1.0)
    assert all(s.score <= best.score + 1e-9 for s in scores)


def test_scores_do_not_depend_on_query_order_and_are_seeded() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    q = [onto.index[C], onto.index[A2]]
    first, _ = phenosim.score_diseases(onto, corpus, q, seed=7, permutations=50)
    second, _ = phenosim.score_diseases(onto, corpus, list(reversed(q)), seed=7,
                                        permutations=50)
    assert first == second


def test_p_values_are_never_zero() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    scores, info = phenosim.score_diseases(onto, corpus, [onto.index[A1]], seed=3,
                                           permutations=10)
    assert min(s.p_value for s in scores) >= info["min_attainable_p"] == pytest.approx(1 / 11)


def test_an_empty_query_or_too_small_a_pool_is_refused() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    with pytest.raises(ValueError, match="no usable term"):
        phenosim.score_diseases(onto, corpus, [], seed=1)
    with pytest.raises(ValueError, match="smaller than the query"):
        phenosim.score_diseases(onto, corpus, [onto.index[A1], onto.index[B1]], seed=1,
                                pool=[onto.index[C]])


def test_term_ratio_counts_only_the_query_terms_lineage() -> None:
    onto = _onto()
    corpus = _corpus(onto)
    q = [onto.index[A1]]
    ic = corpus.ic
    assert phenosim.term_ratio(onto, corpus, q, onto.index[A1]) == pytest.approx(1.0)
    assert phenosim.term_ratio(onto, corpus, q, onto.index[A]) == pytest.approx(
        ic[onto.index[A]] / ic[onto.index[A1]])
    # A sibling shares an ancestor but is a different feature.
    assert phenosim.term_ratio(onto, corpus, q, onto.index[A2]) == 0.0
    assert phenosim.term_ratio(onto, corpus, q, onto.index[B1]) == 0.0


# ------------------------------------------------------------------- the document


def _docx(path: Path, rows, outside="") -> Path:
    """A minimal .docx: one table of rows (lists of cell strings), optional body text."""
    def cell(text):
        runs = "".join(f"<w:r><w:t>{part}</w:t></w:r>" for part in text.split("|"))
        return f"<w:tc><w:p>{runs}</w:p></w:tc>"

    table = "<w:tbl>" + "".join("<w:tr>" + "".join(cell(c) for c in r) + "</w:tr>"
                                for r in rows) + "</w:tbl>"
    body = f"<w:p><w:r><w:t>{outside}</w:t></w:r></w:p>" if outside else ""
    xml = ('<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://schemas.'
           'openxmlformats.org/wordprocessingml/2006/main"><w:body>' + body + table +
           "</w:body></w:document>")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("word/document.xml", xml)
    return path


HEADER = ["Clinical Feature", "HPO Term", "HPO ID", "Presentation / Notes"]


def test_read_terms_takes_ids_from_table_rows_in_order(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", B1, "invented note"],
                                      ["f", "t", A1, "invented"], ["f", "t", B1, "dup"]])
    terms, counts = phenotype.read_terms(doc)
    assert terms == [B1, A1]
    assert counts == {"narrative_negations": 0, "document_reviewed": False}


def test_an_id_split_across_runs_is_joined(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", "HP:90|00011", "n"]])
    assert phenotype.read_terms(doc)[0] == [A1]


def test_a_negated_row_is_refused_without_quoting_it(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", A1, "ok"],
                                      ["SECRETFEATURE", "t", B1, "absent on exam"]])
    with pytest.raises(ValueError, match="row 3") as err:
        phenotype.read_terms(doc)
    assert "SECRETFEATURE" not in str(err.value) and B1 not in str(err.value)


def test_ids_outside_a_table_are_refused(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", A1, "n"]], outside=f"see {C}")
    with pytest.raises(ValueError, match="outside a table"):
        phenotype.read_terms(doc)


def test_a_document_without_ids_is_refused(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", "", "n"]])
    with pytest.raises(ValueError, match="no HPO id"):
        phenotype.read_terms(doc)


def test_an_hpo_term_name_may_contain_a_negation_word(tmp_path) -> None:
    """"Absent speech" is a feature that is PRESENT; HPO's wording may not veto its row."""
    doc = _docx(tmp_path / "x.docx", [HEADER, ["invented feature", "Absent invented thing",
                                               A1, "documented"]])
    assert phenotype.read_terms(doc)[0] == [A1]


def test_a_negation_in_a_clinical_cell_still_refuses_the_row(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["invented", "Invented A1", A1,
                                               "excluded on examination"]])
    with pytest.raises(ValueError, match="row 2"):
        phenotype.read_terms(doc)


def test_narrative_negation_needs_a_persons_review(tmp_path) -> None:
    """Prose can negate a block of features, and no word list can tell that from
    "no single feature is diagnostic" -- so it defers to a person."""
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", A1, "n"]],
                outside="The following were excluded from consideration.")
    with pytest.raises(ValueError, match="MVA_PHENOTYPE_REVIEWED") as err:
        phenotype.read_terms(doc)
    assert phenotype.fingerprint(doc) in str(err.value)

    terms, counts = phenotype.read_terms(doc, reviewed_fingerprint=phenotype.fingerprint(doc))
    assert terms == [A1]
    assert counts == {"narrative_negations": 1, "document_reviewed": True}


def test_a_review_does_not_survive_an_edited_document(tmp_path) -> None:
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", A1, "n"]], outside="none excluded")
    stale = phenotype.fingerprint(doc)
    _docx(doc, [HEADER, ["f", "t", A1, "n"], ["f", "t", B1, "n"]], outside="none excluded")
    with pytest.raises(ValueError, match="negation word"):
        phenotype.read_terms(doc, reviewed_fingerprint=stale)


def test_find_document_wants_exactly_one(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        phenotype.find_document({"data_dir": tmp_path})
    _docx(tmp_path / "A_Phenotype_1.docx", [HEADER])
    assert phenotype.find_document({"data_dir": tmp_path}).name == "A_Phenotype_1.docx"
    _docx(tmp_path / "B_Phenotype_2.docx", [HEADER])
    with pytest.raises(ValueError, match="2 documents"):
        phenotype.find_document({"data_dir": tmp_path})


# ------------------------------------------------------------------ OBO and Monarch


def test_parse_obo_reads_terms_obsoletes_and_version() -> None:
    obo = phenotype.parse_obo(OBO)
    assert obo["version"] == "hp/releases/2099-01-01"
    assert obo["parents"][A1] == [A]
    assert "part_of" not in obo["parents"]
    assert obo["replaced_by"] == {"HP:9000099": B1}
    assert obo["alt_ids"] == {"HP:9900001": A}


def test_resolve_follows_alt_ids_and_replacements() -> None:
    obo = phenotype.parse_obo(OBO)
    assert phenotype.resolve(A1, obo) == A1
    assert phenotype.resolve("HP:9900001", obo) == A
    assert phenotype.resolve("HP:9000099", obo) == B1
    assert phenotype.resolve("HP:9999999", obo) is None


def test_resolve_survives_a_replacement_cycle() -> None:
    obo = {"parents": {}, "alt_ids": {"HP:1": "HP:2"}, "replaced_by": {"HP:2": "HP:1"}}
    assert phenotype.resolve("HP:1", obo) is None


MONARCH_HEADER = "subject\tsubject_label\tnegated\tpredicate\tobject\n"


def test_load_associations_drops_negated_and_non_mondo_rows() -> None:
    lines = [MONARCH_HEADER,
             f"MONDO:9000001\tinvented disease one\tFalse\thas_phenotype\t{A1}\n",
             f"MONDO:9000001\tinvented disease one\tTrue\thas_phenotype\t{B1}\n",
             f"HGNC:1\tinvented gene\tFalse\thas_phenotype\t{C}\n",
             "MONDO:9000002\tinvented disease two\tFalse\thas_phenotype\tUPHENO:1\n"]
    annotations, labels, stats = phenotype.load_associations(lines)
    assert annotations == {"MONDO:9000001": {A1}}
    assert labels == {"MONDO:9000001": "invented disease one"}
    assert (stats["negated"], stats["non_mondo_subject"], stats["non_hpo_object"]) == (1, 1, 1)


def test_a_changed_monarch_format_is_refused() -> None:
    with pytest.raises(ValueError, match="'negated'"):
        phenotype.load_associations(["subject\tsubject_label\tobject\n"])


def test_the_cached_names_carry_the_release(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("MVA_REF_ROOT", str(tmp_path / "ref"))
    fetched = []
    monkeypatch.setattr(phenotype.refcache, "fetch",
                        lambda url, dest: fetched.append(dest.name) or {"url": url})
    paths, provenance = phenotype.ensure_sources({})
    assert paths["hpo"].name == "2026-09-01.hp.obo"
    assert paths["monarch"].name == "2026-09-02.disease_phenotype.all.tsv.gz"
    assert provenance["hpo"]["release"] == "2026-09-01"


def test_the_configured_releases_are_pinned() -> None:
    """"latest" would change what a re-run matches against without anyone deciding it."""
    import yaml

    config = yaml.safe_load((Path(__file__).resolve().parent.parent / "config" /
                             "pipeline.yaml").read_text(encoding="utf-8"))
    settings = config["l2"]["channel_d"]
    for key in ("hpo_url", "monarch_url"):
        assert "latest" not in settings[key]
        assert any(part[:2] == "20" and part.count("-") == 2 for part in settings[key].split("/"))


# ------------------------------------------------------------------ enrichment


def test_load_indications_keeps_only_the_requested_stage(monkeypatch) -> None:
    monkeypatch.setattr(enrichment, "_rows", lambda path: [
        {"drugId": "CHEMBL9001", "diseaseId": "MONDO_9000001", "maxClinicalStage": "APPROVAL"},
        {"drugId": "CHEMBL9002", "diseaseId": "MONDO_9000001", "maxClinicalStage": "PHASE_2"},
        {"drugId": None, "diseaseId": "MONDO_9000002", "maxClinicalStage": "APPROVAL"},
    ])
    out, stats = enrichment.load_indications(["p"])
    assert out == {"MONDO_9000001": {"CHEMBL9001"}}
    assert (stats["indication_rows"], stats["indications_at_stage"]) == (3, 1)
    widened, stats = enrichment.load_indications(["p"], stage="PHASE_2")
    assert widened == {"MONDO_9000001": {"CHEMBL9002"}} and stats["stage"] == "PHASE_2"


def test_each_channel_uses_its_own_release_setting(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("MVA_ENRICHMENT_ROOT", str(tmp_path / "enrich"))
    bases = []
    monkeypatch.setattr(enrichment, "part_urls", lambda base: bases.append(base) or [base + "x.parquet"])
    monkeypatch.setattr(enrichment.refcache, "fetch", lambda url, dest: {"url": url})
    enrichment.ensure_datasets({"reference_dir": str(tmp_path / "ref"),
                                "l2": {"channel_d": {"open_targets_release": "99.01"}}},
                               ("clinical_indication",), channel="channel_d")
    assert bases == [f"{enrichment.OPEN_TARGETS_FTP}/99.01/output/clinical_indication/"]


# ------------------------------------------------------------------ ranking


def _hit(disease, score, p=0.001):
    return phenosim.DiseaseScore(disease=disease, score=score, raw=score, p_value=p)


MOLECULES = {c: {"name": f"INVENTED {c}", "drug_type": "Small molecule",
                 "clinical_stage": "APPROVAL"} for c in ("CHEMBL9001", "CHEMBL9002", "CHEMBL9003")}


def test_rank_drugs_takes_the_better_route_and_breaks_ties_deterministically() -> None:
    indications = {"MONDO_9000001": {"CHEMBL9001", "CHEMBL9002"},
                   "HP_9000011": {"CHEMBL9002", "CHEMBL9003"}}
    entries = channel.rank_drugs([_hit("MONDO:9000001", 0.4)], {A1: 0.9}, indications,
                                 MOLECULES, approved_only=True)
    assert [(e["record"].chembl_id, e["route"], e["score"]) for e in entries] == [
        ("CHEMBL9002", "phenotype", 0.9),   # both routes -> more support wins the tie
        ("CHEMBL9003", "phenotype", 0.9),
        ("CHEMBL9001", "disease", 0.4),
    ]


def test_rank_drugs_drops_unapproved_drugs_when_asked() -> None:
    molecules = {**MOLECULES, "CHEMBL9001": {**MOLECULES["CHEMBL9001"],
                                             "clinical_stage": "PHASE_3"}}
    entries = channel.rank_drugs([_hit("MONDO:9000001", 0.4)], {},
                                 {"MONDO_9000001": {"CHEMBL9001"}}, molecules,
                                 approved_only=True)
    assert entries == []


def test_query_terms_reports_counts_not_terms() -> None:
    obo = phenotype.parse_obo(OBO)
    onto = phenosim.Ontology.build(obo["parents"])
    indices, counts = channel.query_terms([A1, MOI, "HP:9999999", "HP:9000099"], obo, onto)
    assert sorted(onto.ids[i] for i in indices) == sorted([A1, B1])
    assert counts == {"terms_in_document": 4, "terms_used": 2, "terms_unresolved": 1,
                      "terms_outside_phenotypic_abnormality": 1}
    assert not any(str(v).startswith("HP:") for v in counts.values())


# ------------------------------------------------------------------ end to end


def _stub_sources(tmp_path, monkeypatch, *, indications):
    ref = tmp_path / "ref"
    ref.mkdir()
    (ref / "hp.obo").write_text("".join(OBO), encoding="utf-8")
    monarch = ref / "assoc.tsv.gz"
    rows = [MONARCH_HEADER] + [
        f"{d}\tinvented {d}\tFalse\thas_phenotype\t{t}\n"
        for d, ts in {"MONDO:9000001": [A1, B1], "MONDO:9000002": [A2],
                      "MONDO:9000003": [B1, C], "MONDO:9000004": [C],
                      "MONDO:9000005": [A2, C], "MONDO:9000006": [B]}.items()
        for t in ts]
    with gzip.open(monarch, "wt", encoding="utf-8") as handle:
        handle.writelines(rows)
    monkeypatch.setattr(channel.phenotype, "ensure_sources",
                        lambda config: ({"hpo": ref / "hp.obo", "monarch": monarch},
                                        {"hpo": {"release": "2099-01-01"}}))
    monkeypatch.setattr(channel.enrichment, "ensure_datasets",
                        lambda config, datasets, channel: ({"clinical_indication": ["i"],
                                                            "drug_molecule": ["m"]}, []))
    monkeypatch.setattr(channel.enrichment, "load_indications",
                        lambda paths, stage=enrichment.APPROVED: (
                            indications, {"indication_rows": 1, "stage": stage}))
    monkeypatch.setattr(channel.enrichment, "load_molecules", lambda paths: MOLECULES)

    data = tmp_path / "data"
    data.mkdir()
    _docx(data / "Invented_Phenotype_1.docx", [HEADER, ["f", "t", A1, "documented"],
                                              ["f", "t", B1, "documented"],
                                              ["f", "t", MOI, "documented"]])
    return {"data_dir": data, "results_dir": tmp_path / "results", "seed": 42,
            "l2": {"channel_d": {"permutations": 30, "p_max": 1.0, "top_diseases": 3}}}


def test_generate_writes_a_symptomatic_table_with_no_term_or_disease(tmp_path, monkeypatch,
                                                                     caplog) -> None:
    config = _stub_sources(tmp_path, monkeypatch, indications={
        "MONDO_9000001": {"CHEMBL9001"}, "HP_9000011": {"CHEMBL9002"}})
    with caplog.at_level("DEBUG"):
        channel.generate(config)

    out = tmp_path / "results" / "l2" / "channel_d_phenotype"
    table = (out / "candidates.tsv").read_text(encoding="utf-8")
    lines = table.strip().split("\n")
    assert lines[0].split("\t")[-1] == "endpoint"
    assert all(line.endswith("\tsymptomatic") for line in lines[1:])
    assert {line.split("\t")[1] for line in lines[1:]} == {"CHEMBL9001", "CHEMBL9002"}

    # Patient-derived identifiers stay out of the table, the channel record and the log.
    for text in (table, (out / "channel.json").read_text(encoding="utf-8"), caplog.text):
        assert "HP:9" not in text and "HP_9" not in text and "MONDO:9" not in text

    evidence = json.loads((out / "evidence.json").read_text(encoding="utf-8"))
    assert evidence["patient_derived"] is True and evidence["redistributable"] is False
    record = json.loads((out / "channel.json").read_text(encoding="utf-8"))
    # The mean information content of the subject's own terms is a number someone holding
    # the public files could recompute for a guessed feature set and match; it belongs
    # with the patient-derived evidence, not in the shareable record.
    assert "self_similarity" not in json.dumps(record)
    assert "query_self_similarity" in evidence
    assert record["parameters"]["indication_stage"] == "APPROVAL"
    assert record["counts"]["terms_used"] == 2
    assert record["counts"]["terms_outside_phenotypic_abnormality"] == 1
    assert "2099-01-01" in record["attribution"]     # HPO licence: show the version


def test_an_unattainable_p_max_is_a_config_error(tmp_path, monkeypatch) -> None:
    """With too few permutations no disease can pass; say so before doing any work."""
    config = _stub_sources(tmp_path, monkeypatch, indications={})
    config["l2"]["channel_d"].update(permutations=99, p_max=0.001)
    with pytest.raises(ValueError, match="smallest attainable"):
        channel.generate(config)


def test_generate_refuses_to_write_an_empty_channel(tmp_path, monkeypatch) -> None:
    config = _stub_sources(tmp_path, monkeypatch, indications={"MONDO_9999999": {"CHEMBL9001"}})
    with pytest.raises(ValueError, match="unmet threshold"):
        channel.generate(config)
    assert not (tmp_path / "results" / "l2" / "channel_d_phenotype" / "candidates.tsv").exists()


def test_a_word_lock_file_does_not_count_as_a_second_document(tmp_path) -> None:
    """Word writes ~$name.docx while a document is open -- including while it is reviewed."""
    _docx(tmp_path / "Invented_Phenotype_1.docx", [HEADER, ["f", "t", A1, "n"]])
    (tmp_path / "~$vented_Phenotype_1.docx").write_bytes(b"word owner file")
    assert phenotype.find_document({"data_dir": tmp_path}).name == "Invented_Phenotype_1.docx"


def test_narrative_text_is_not_confused_with_cell_text(tmp_path) -> None:
    """Prose is identified structurally, not by subtracting cell strings from the text.

    A note repeating a cell's wording must not cancel out a real negation elsewhere.
    """
    doc = _docx(tmp_path / "x.docx", [HEADER, ["f", "t", A1, "documented"]],
                outside="documented documented none excluded")
    with pytest.raises(ValueError, match="negation word"):
        phenotype.read_terms(doc)
    _, counts = phenotype.read_terms(doc, reviewed_fingerprint=phenotype.fingerprint(doc))
    assert counts["narrative_negations"] == 2
