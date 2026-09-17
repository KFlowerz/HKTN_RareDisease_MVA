"""Tests for Channel B's wiring and the licence boundary it depends on.

Every drug, gene and protein identifier here is invented. The proximity measure itself is
tested in ``test_l2_proximity.py``; this file is about what crosses the boundary between
the restricted enrichment zone and a published output.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src import refcache
from src.l2_channels import channel_b_proximity as channel, enrichment

REPO_ROOT = Path(__file__).resolve().parent.parent


def _record(chembl_id="CHEMBL1", **kwargs) -> enrichment.DrugRecord:
    defaults = dict(name="INVENTEDAZOLE", drug_type="Small molecule",
                    clinical_stage=enrichment.APPROVED,
                    targets=frozenset({"9606.ENSP00000000001"}),
                    unmapped_targets=frozenset({"ENSG00000000009"}))
    defaults.update(kwargs)
    return enrichment.DrugRecord(chembl_id=chembl_id, **defaults)


# ------------------------------------------------------------ the boundary


def test_publishable_returns_exactly_the_whitelist() -> None:
    assert set(enrichment.publishable(_record())) == set(enrichment.PUBLISHABLE_FIELDS)


def test_publishable_never_returns_targets() -> None:
    """The one rule this zone exists for."""
    public = enrichment.publishable(_record())
    assert "targets" not in public and "unmapped_targets" not in public
    flat = " ".join(str(v) for v in public.values())
    assert "ENSP" not in flat and "ENSG" not in flat


def test_publishable_is_a_whitelist_not_a_blacklist() -> None:
    """A field added to DrugRecord later must stay inside the zone until someone decides.

    A blacklist would leak every new field by default; that is the failure mode this is
    shaped to avoid, so it is asserted rather than left to the docstring.
    """
    source = (REPO_ROOT / "src" / "l2_channels" / "enrichment.py").read_text(encoding="utf-8")
    assert "PUBLISHABLE_FIELDS" in source
    for field in enrichment.PUBLISHABLE_FIELDS:
        assert field in ("chembl_id", "name", "clinical_stage", "drug_type", "n_targets")


def test_the_published_count_is_a_count_not_the_targets() -> None:
    assert enrichment.publishable(_record())["n_targets"] == 1


def test_the_restricted_source_is_marked_non_redistributable() -> None:
    assert enrichment.REDISTRIBUTABLE is False
    assert "ChEMBL" in enrichment.LICENCE and "BY-SA" in enrichment.LICENCE


def test_the_enrichment_cache_is_separate_from_the_public_one() -> None:
    """Two zones, two directories: the separation is physical, not a convention."""
    assert enrichment.DEFAULT_DIR != refcache.DEFAULT_DIR


def test_the_enrichment_cache_may_not_live_inside_the_repo(monkeypatch) -> None:
    monkeypatch.delenv("MVA_ENRICHMENT_ROOT", raising=False)
    with pytest.raises(ValueError, match="inside the repository"):
        enrichment.enrichment_dir({"enrichment_dir": str(REPO_ROOT / "enrichment")})


# ------------------------------------------------------------ parsing


@pytest.mark.parametrize("value,expected", [
    (["a", "b"], ["a", "b"]),
    ("['a', 'b']", ["a", "b"]),
    ("plain", ["plain"]),
    (None, []),
    ([], []),
])
def test_as_list_survives_a_stringified_column(value, expected) -> None:
    assert enrichment._as_list(value) == expected


def _fake_rows(monkeypatch, by_path: dict):
    monkeypatch.setattr(enrichment, "_rows", lambda path: by_path[path])


def test_load_mechanisms_accumulates_targets_per_drug(monkeypatch) -> None:
    """A mechanism row names one target and the drugs acting on it."""
    _fake_rows(monkeypatch, {"p": [
        {"targets": ["ENSG1"], "chemblIds": ["CHEMBL1", "CHEMBL2"]},
        {"targets": ["ENSG2"], "chemblIds": ["CHEMBL1"]},
    ]})
    mapping, rows = enrichment.load_mechanisms(["p"], {"ENSG1": "P1", "ENSG2": "P2"})
    assert rows == 2
    assert mapping["CHEMBL1"][0] == frozenset({"P1", "P2"})
    assert mapping["CHEMBL2"][0] == frozenset({"P1"})


def test_load_mechanisms_records_targets_it_could_not_map(monkeypatch) -> None:
    """A gene absent from the interactome is reported, not silently dropped."""
    _fake_rows(monkeypatch, {"p": [{"targets": ["ENSG1", "ENSGX"], "chemblIds": ["CHEMBL1"]}]})
    mapping, _ = enrichment.load_mechanisms(["p"], {"ENSG1": "P1"})
    targets, unmapped = mapping["CHEMBL1"]
    assert targets == frozenset({"P1"}) and unmapped == frozenset({"ENSGX"})


def test_load_molecules_keys_on_the_chembl_id(monkeypatch) -> None:
    _fake_rows(monkeypatch, {"m": [
        {"id": "CHEMBL1", "name": "INVENTEDAZOLE", "drugType": "Small molecule",
         "maximumClinicalStage": "APPROVAL"},
        {"id": None, "name": "ignored"},
    ]})
    molecules = enrichment.load_molecules(["m"])
    assert list(molecules) == ["CHEMBL1"]
    assert molecules["CHEMBL1"]["clinical_stage"] == "APPROVAL"


def test_name_index_prefers_a_name_over_another_molecule_s_synonym(monkeypatch) -> None:
    _fake_rows(monkeypatch, {"m": [
        {"id": "CHEMBL1", "name": "INVENTEDAZOLE", "synonyms": ["inv-1"]},
        {"id": "CHEMBL2", "name": "OTHERAZOLE", "synonyms": ["inventedazole"]},
    ]})
    index, stats = enrichment.name_index(["m"])
    assert index["inventedazole"] == "CHEMBL1"
    assert index["inv 1"] == "CHEMBL1"
    # CHEMBL2's synonym lost to CHEMBL1's preferred name, which settles the key outright.
    assert "inventedazole" not in stats["ambiguous_synonyms"]


def test_name_index_reports_two_molecules_sharing_a_preferred_name(monkeypatch) -> None:
    """The worst collision, not the benign one: picking either is a coin toss on identity.

    An earlier version subtracted preferred names from the ambiguous set, which reported
    exactly this case as unambiguous.
    """
    _fake_rows(monkeypatch, {"m": [
        {"id": "CHEMBL2", "name": "INVENTEDAZOLE"},
        {"id": "CHEMBL1", "name": "inventedazole"},
    ]})
    index, stats = enrichment.name_index(["m"])
    assert "inventedazole" in stats["ambiguous_synonyms"]
    assert index["inventedazole"] == "CHEMBL1", "the winner must not depend on row order"


def test_name_index_is_independent_of_row_order(monkeypatch) -> None:
    rows = [{"id": "CHEMBL9", "synonyms": ["shared"]},
            {"id": "CHEMBL3", "synonyms": ["shared"]}]
    _fake_rows(monkeypatch, {"m": rows})
    forwards, _ = enrichment.name_index(["m"])
    _fake_rows(monkeypatch, {"m": list(reversed(rows))})
    backwards, stats = enrichment.name_index(["m"])
    assert forwards["shared"] == backwards["shared"] == "CHEMBL3"
    assert "shared" in stats["ambiguous_synonyms"]


def test_ensembl_index_filters_on_the_alias_source(tmp_path) -> None:
    path = tmp_path / "aliases.txt"
    path.write_text(
        "#string_protein_id\talias\tsource\n"
        "P1\tENSG1\tEnsembl_gene\n"
        "P2\tENSG2\tSomething_else\n", encoding="utf-8")
    assert enrichment.ensembl_index(path) == {"ENSG1": "P1"}


def test_drug_records_keeps_only_approved_when_asked(monkeypatch) -> None:
    monkeypatch.setattr(enrichment, "ensure_datasets",
                        lambda config: ({"drug_mechanism_of_action": ["x"],
                                         "drug_molecule": ["m"]}, []))
    monkeypatch.setattr(enrichment, "load_mechanisms", lambda paths, index: (
        {"CHEMBL1": (frozenset({"P1"}), frozenset()),
         "CHEMBL2": (frozenset({"P2"}), frozenset())}, 2))
    monkeypatch.setattr(enrichment, "load_molecules", lambda paths: {
        "CHEMBL1": {"name": "A", "drug_type": "Small molecule", "clinical_stage": "APPROVAL"},
        "CHEMBL2": {"name": "B", "drug_type": "Small molecule", "clinical_stage": "PHASE_2"}})

    approved, _, stats = enrichment.drug_records(
        {"l2": {"channel_b": {"approved_only": True}}}, {})
    assert [r.chembl_id for r in approved] == ["CHEMBL1"] and stats["drugs_kept"] == 1

    everything, _, _ = enrichment.drug_records(
        {"l2": {"channel_b": {"approved_only": False}}}, {})
    assert [r.chembl_id for r in everything] == ["CHEMBL1", "CHEMBL2"]


def test_a_drug_with_no_mappable_target_is_not_scored(monkeypatch) -> None:
    """Scoring it would put a drug in the ranking on the strength of nothing."""
    monkeypatch.setattr(enrichment, "ensure_datasets",
                        lambda config: ({"drug_mechanism_of_action": ["x"],
                                         "drug_molecule": ["m"]}, []))
    monkeypatch.setattr(enrichment, "load_mechanisms", lambda paths, index: (
        {"CHEMBL1": (frozenset(), frozenset({"ENSGX"}))}, 1))
    monkeypatch.setattr(enrichment, "load_molecules", lambda paths: {
        "CHEMBL1": {"name": "A", "drug_type": "Small molecule", "clinical_stage": "APPROVAL"}})
    records, _, stats = enrichment.drug_records({}, {})
    assert records == [] and stats["drugs_with_unmappable_targets_only"] == 1


def test_approved_is_read_from_the_clinical_stage() -> None:
    assert _record().approved is True
    assert _record(clinical_stage="PHASE_2").approved is False


# ------------------------------------------------------------ the channel


def test_the_channel_needs_l1_to_have_run(tmp_path) -> None:
    """Inventing a module would still produce a ranking -- a meaningless one."""
    with pytest.raises(FileNotFoundError, match="L1 must run first"):
        channel._module_nodes(tmp_path)


def test_module_nodes_are_read_by_column_name(tmp_path) -> None:
    l1 = tmp_path / "l1_target"
    l1.mkdir()
    (l1 / "module_nodes.tsv").write_text(
        "protein\tsymbol\trank\n9606.ENSP1\tGENE1\t1\n9606.ENSP2\tGENE2\t2\n", encoding="utf-8")
    assert channel._module_nodes(tmp_path) == ["9606.ENSP1", "9606.ENSP2"]


def test_the_written_table_carries_no_target_identifier(tmp_path, monkeypatch) -> None:
    """The end-to-end form of the boundary rule, asserted on the bytes that get written."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from((f"P{i}", f"P{i + 1}") for i in range(30))

    l1 = tmp_path / "l1_target"
    l1.mkdir()
    (l1 / "module_nodes.tsv").write_text("protein\n" + "\n".join(f"P{i}" for i in range(3)),
                                         encoding="utf-8")

    monkeypatch.setattr(channel.sources, "ensure_sources",
                        lambda config: ({"links": "x", "aliases": "y"}, []))
    monkeypatch.setattr(channel.l1_module, "build_graph",
                        lambda edges, score_min=0: graph)
    monkeypatch.setattr(channel.sources, "load_edges", lambda path, score_min=0: [])
    monkeypatch.setattr(channel.enrichment, "ensembl_index", lambda path: {})
    monkeypatch.setattr(channel.enrichment, "drug_records", lambda config, index: (
        [_record("CHEMBL1", targets=frozenset({"P5"})),
         _record("CHEMBL2", targets=frozenset({"P25"}))],
        [{"source": "restricted", "redistributable": False}],
        {"approved_only": True, "mechanism_rows": 2, "drugs_with_a_mechanism": 2,
         "molecules": 2, "drugs_kept": 2, "drugs_with_unmappable_targets_only": 0}))

    channel.generate({"results_dir": str(tmp_path), "seed": 42,
                      "l2": {"channel_b": {"permutations": 50, "min_bin_size": 5}}})

    written = (tmp_path / "l2" / "channel_b_proximity" / "candidates.tsv").read_text(
        encoding="utf-8")
    assert "ENSP" not in written and "ENSG" not in written
    assert "P5" not in written and "P25" not in written
    assert "CHEMBL1" in written and "INVENTEDAZOLE" in written


def test_the_closer_drug_ranks_first(tmp_path, monkeypatch) -> None:
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from((f"P{i}", f"P{i + 1}") for i in range(30))
    l1 = tmp_path / "l1_target"
    l1.mkdir()
    (l1 / "module_nodes.tsv").write_text("protein\n" + "\n".join(f"P{i}" for i in range(3)),
                                         encoding="utf-8")

    monkeypatch.setattr(channel.sources, "ensure_sources",
                        lambda config: ({"links": "x", "aliases": "y"}, []))
    monkeypatch.setattr(channel.l1_module, "build_graph", lambda edges, score_min=0: graph)
    monkeypatch.setattr(channel.sources, "load_edges", lambda path, score_min=0: [])
    monkeypatch.setattr(channel.enrichment, "ensembl_index", lambda path: {})
    monkeypatch.setattr(channel.enrichment, "drug_records", lambda config, index: (
        [_record("FAR", targets=frozenset({"P25"})),
         _record("NEAR", targets=frozenset({"P4"}))],
        [], {"approved_only": True, "mechanism_rows": 2, "drugs_with_a_mechanism": 2,
             "molecules": 2, "drugs_kept": 2, "drugs_with_unmappable_targets_only": 0}))

    channel.generate({"results_dir": str(tmp_path), "seed": 42,
                      "l2": {"channel_b": {"permutations": 200, "min_bin_size": 5}}})
    lines = (tmp_path / "l2" / "channel_b_proximity" / "candidates.tsv").read_text(
        encoding="utf-8").splitlines()
    assert lines[1].split("\t")[1] == "NEAR"


def test_the_channel_names_no_gene_drug_or_disease() -> None:
    """Channel B must work for a disease in no knowledge graph, so it names none."""
    forbidden = ("BUB1B", "CEP57", "TRIP13", "BUB3", "CEP192", "aneuploidy",
                 "mosaic variegated", "spindle")
    for name in ("channel_b_proximity.py", "enrichment.py", "proximity.py"):
        source = (REPO_ROOT / "src" / "l2_channels" / name).read_text(encoding="utf-8")
        for term in forbidden:
            assert term.lower() not in source.lower(), f"{term!r} in {name}"


# ------------------------------------------- regressions found in code review


def test_a_drug_that_lost_every_target_is_counted(monkeypatch) -> None:
    """The statistic that would expose a broken crosswalk must be able to fire.

    load_mechanisms previously returned only drugs with at least one mapped target, so
    drugs_with_unmappable_targets_only was structurally always 0 -- the one number that
    reports "the crosswalk matched nothing" could never be anything but zero.
    """
    monkeypatch.setattr(enrichment, "_rows", lambda path: [
        {"targets": ["ENSG1"], "chemblIds": ["CHEMBL1"]},
        {"targets": ["ENSGX"], "chemblIds": ["CHEMBL2"]},
    ])
    mapping, _ = enrichment.load_mechanisms(["p"], {"ENSG1": "P1"})
    assert mapping["CHEMBL2"] == (frozenset(), frozenset({"ENSGX"}))

    monkeypatch.setattr(enrichment, "ensure_datasets",
                        lambda config: ({"drug_mechanism_of_action": ["p"],
                                         "drug_molecule": ["m"]}, []))
    monkeypatch.setattr(enrichment, "load_molecules", lambda paths: {
        c: {"name": c, "drug_type": "Small molecule", "clinical_stage": "APPROVAL"}
        for c in ("CHEMBL1", "CHEMBL2")})
    records, _, stats = enrichment.drug_records({}, {"ENSG1": "P1"})
    assert [r.chembl_id for r in records] == ["CHEMBL1"]
    assert stats["drugs_with_unmappable_targets_only"] == 1


def test_an_empty_release_listing_is_refused(monkeypatch, tmp_path) -> None:
    """Zero parts downstream reads as "no drug has a target", which is a finding."""
    monkeypatch.setenv("MVA_ENRICHMENT_ROOT", str(tmp_path / "enrich"))
    monkeypatch.setattr(enrichment, "part_urls", lambda base: [])
    with pytest.raises(ValueError, match="no Parquet part"):
        enrichment.ensure_datasets({"reference_dir": str(tmp_path / "ref")})


def test_the_two_zones_may_not_be_the_same_directory(monkeypatch, tmp_path) -> None:
    """Collapsing them by configuration would leave decision D7 with no boundary."""
    monkeypatch.delenv("MVA_ENRICHMENT_ROOT", raising=False)
    monkeypatch.delenv("MVA_REF_ROOT", raising=False)
    shared = str(tmp_path / "one-cache")
    with pytest.raises(ValueError, match="must not share a directory"):
        enrichment.enrichment_dir({"enrichment_dir": shared, "reference_dir": shared})


def _l1_artifacts(tmp_path, nodes=("P0", "P1", "P2"), score_min=700):
    import json

    l1 = tmp_path / "l1_target"
    l1.mkdir(exist_ok=True)
    (l1 / "module_nodes.tsv").write_text("protein\n" + "\n".join(nodes) + "\n",
                                         encoding="utf-8")
    (l1 / "module.json").write_text(
        json.dumps({"parameters": {"string_score_min": score_min}}), encoding="utf-8")
    return l1


def test_the_cutoff_comes_from_the_module_not_the_config(tmp_path) -> None:
    """Otherwise the module and the graph it is scored against can silently disagree."""
    _l1_artifacts(tmp_path, score_min=900)
    assert channel._module_score_min(tmp_path, {"l1": {"string_score_min": 400}}) == 900


def test_the_cutoff_falls_back_to_config_when_unrecorded(tmp_path) -> None:
    l1 = tmp_path / "l1_target"
    l1.mkdir()
    (l1 / "module_nodes.tsv").write_text("protein\nP0\n", encoding="utf-8")
    assert channel._module_score_min(tmp_path, {"l1": {"string_score_min": 400}}) == 400


def test_an_empty_module_is_refused(tmp_path) -> None:
    l1 = tmp_path / "l1_target"
    l1.mkdir()
    (l1 / "module_nodes.tsv").write_text("protein\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no module protein"):
        channel._module_nodes(tmp_path)


def _stub_channel(monkeypatch, graph, records):
    monkeypatch.setattr(channel.sources, "ensure_sources",
                        lambda config: ({"links": "x", "aliases": "y"}, []))
    monkeypatch.setattr(channel.l1_module, "build_graph", lambda edges, score_min=0: graph)
    monkeypatch.setattr(channel.sources, "load_edges", lambda path, score_min=0: [])
    monkeypatch.setattr(channel.enrichment, "ensembl_index", lambda path: {})
    monkeypatch.setattr(channel.enrichment, "drug_records", lambda config, index: (
        records, [], {"approved_only": True, "mechanism_rows": 2,
                      "drugs_with_a_mechanism": 2, "molecules": 2,
                      "drugs_kept": len(records),
                      "drugs_with_unmappable_targets_only": 0}))


def test_no_scorable_drug_is_an_error_not_an_empty_table(tmp_path, monkeypatch) -> None:
    """An empty channel is a claim -- "no drug is near this module" -- and this is not it."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from((f"P{i}", f"P{i + 1}") for i in range(30))
    _l1_artifacts(tmp_path)
    _stub_channel(monkeypatch, graph, [])
    with pytest.raises(ValueError, match="broken crosswalk"):
        channel.generate({"results_dir": str(tmp_path), "seed": 42,
                          "l2": {"channel_b": {"permutations": 20, "min_bin_size": 5}}})


def test_a_module_absent_from_the_graph_is_an_error(tmp_path, monkeypatch) -> None:
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from((f"Q{i}", f"Q{i + 1}") for i in range(30))
    _l1_artifacts(tmp_path)
    _stub_channel(monkeypatch, graph, [_record("CHEMBL1", targets=frozenset({"Q5"}))])
    with pytest.raises(ValueError, match="largest component"):
        channel.generate({"results_dir": str(tmp_path), "seed": 42,
                          "l2": {"channel_b": {"permutations": 20, "min_bin_size": 5}}})
