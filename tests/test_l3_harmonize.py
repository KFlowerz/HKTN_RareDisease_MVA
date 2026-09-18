"""Tests for L3 drug-identity harmonisation.

Every molecule, ChEMBL id, UNII and RxCUI here is invented. ChEMBL ids use the
``CHEMBL9xxxxxx`` range, UNIIs and RxCUIs are made up, and no test asserts anything about
a real drug.
"""

from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path

import pytest

from src.l2_channels import enrichment
from src.l3_integrate import harmonize


def _molecules(monkeypatch, rows):
    monkeypatch.setattr(enrichment, "_rows", lambda path: rows)


def _unichem(tmp_path: Path, mapping) -> Path:
    path = tmp_path / "src1src14.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("From src:'1'\tTo src:'14'\n")
        for chembl_id, unii in mapping:
            handle.write(f"{chembl_id}\t{unii}\n")
    return path


def _openfda(tmp_path: Path, results, name: str = "ndc") -> Path:
    path = tmp_path / f"{name}.json.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{name}.json", json.dumps({"results": results}))
    return path


# --------------------------------------------------------------- parent collapsing


def test_a_formulation_resolves_to_its_parent():
    parent = {"CHEMBL9000002": "CHEMBL9000001"}
    assert harmonize.resolve_parent("CHEMBL9000002", parent) == "CHEMBL9000001"
    assert harmonize.resolve_parent("CHEMBL9000001", parent) == "CHEMBL9000001"


def test_a_chain_of_parents_resolves_to_the_root():
    parent = {"C3": "C2", "C2": "C1"}
    assert harmonize.resolve_parent("C3", parent) == "C1"


def test_a_cycle_in_the_source_data_terminates():
    """A hang mid-run is worse than a wrong-but-visible answer."""
    assert harmonize.resolve_parent("C1", {"C1": "C2", "C2": "C1"}) in {"C1", "C2"}


def test_parent_index_reads_the_curated_relation(monkeypatch):
    _molecules(monkeypatch, [
        {"id": "CHEMBL9000001", "name": "INVENTIB", "drugType": "Small molecule",
         "maximumClinicalStage": "APPROVAL"},
        {"id": "CHEMBL9000002", "name": "INVENTIB CITRATE", "parentId": "CHEMBL9000001"},
        {"id": None, "name": "ignored"},
    ])
    parent, identity, stats = harmonize.parent_index(["p"])
    assert parent == {"CHEMBL9000002": "CHEMBL9000001"}
    assert identity["CHEMBL9000001"]["clinical_stage"] == "APPROVAL"
    assert stats == {"molecules": 2, "with_a_parent": 1}


# ------------------------------------------------------------------- the crosswalk


def test_reads_the_unichem_mapping(tmp_path):
    path = _unichem(tmp_path, [("CHEMBL9000001", "UNII0001"),
                               ("CHEMBL9000001", "UNII0002"),
                               ("CHEMBL9000003", "UNII0003")])
    mapping = harmonize.load_chembl_to_unii(path)
    assert mapping["CHEMBL9000001"] == {"UNII0001", "UNII0002"}


def test_refuses_a_file_that_is_not_a_unichem_mapping(tmp_path):
    """A silently wrong crosswalk would mis-identify every drug in the pipeline."""
    path = tmp_path / "wrong.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("chembl_id,something\n")
    with pytest.raises(ValueError, match="UniChem"):
        harmonize.load_chembl_to_unii(path)


def test_openfda_indexes_by_unii_and_by_name(tmp_path):
    path = _openfda(tmp_path, [
        {"generic_name": "Inventib", "brand_name": "Inventor",
         "openfda": {"rxcui": ["999001"], "unii": ["UNII0001"],
                     "pharm_class_moa": ["Invented Inhibitor [MoA]"]}},
        {"generic_name": "NoCode", "openfda": {"unii": ["UNII0009"]}},
    ])
    by_unii, by_name, stats = harmonize.load_openfda([path])
    assert by_unii["UNII0001"]["rxcui"] == {"999001"}
    assert by_name["inventib"]["rxcui"] == {"999001"}
    assert by_name["inventor"]["rxcui"] == {"999001"}
    assert "nocode" not in by_name, "a record with no RxCUI carries no crosswalk"
    assert stats["openfda_records"] == 2 and stats["records_with_rxcui"] == 1


# -------------------------------------------------------------- the rxcui merge


def _identity(key, name, rxcui, chembl_ids=None):
    return harmonize.Identity(key=key, name=name, rxcui=tuple(rxcui),
                              chembl_ids=frozenset(chembl_ids or {key}))


def test_identities_sharing_a_complete_rxcui_set_are_merged():
    """On the real data ixazomib and ixazomib citrate took ranks 2 and 3 as two rows."""
    identities = {
        "C1": _identity("C1", "INVENTIB", ("999001", "999002")),
        "C2": _identity("C2", "INVENTIB CITRATE", ("999001", "999002")),
    }
    merged, count = harmonize.merge_on_rxcui(identities)
    assert count == 1
    assert merged["C1"].key == merged["C2"].key
    assert merged["C1"].name == "INVENTIB", "the shorter name is the parent's"
    assert merged["C1"].chembl_ids == frozenset({"C1", "C2"})


def test_a_partial_rxcui_overlap_is_not_merged():
    """A combination product legitimately contains its components' codes."""
    identities = {
        "C1": _identity("C1", "INVENTIB", ("999001",)),
        "C2": _identity("C2", "INVENTIB AND OTHERAZOLE", ("999001", "999002")),
    }
    merged, count = harmonize.merge_on_rxcui(identities)
    assert count == 0
    assert merged["C1"].key != merged["C2"].key


def test_identities_with_no_rxcui_are_never_merged():
    """Empty sets are all equal; merging on them would collapse unrelated drugs."""
    identities = {"C1": _identity("C1", "AAA", ()), "C2": _identity("C2", "BBB", ())}
    merged, count = harmonize.merge_on_rxcui(identities)
    assert count == 0 and merged["C1"].key != merged["C2"].key


def test_the_merge_winner_does_not_depend_on_iteration_order():
    a = {"C2": _identity("C2", "SAME", ("1",)), "C1": _identity("C1", "SAME", ("1",))}
    b = {"C1": _identity("C1", "SAME", ("1",)), "C2": _identity("C2", "SAME", ("1",))}
    assert harmonize.merge_on_rxcui(a)[0]["C1"].key == \
        harmonize.merge_on_rxcui(b)[0]["C1"].key == "C1"


# --------------------------------------------------------------------- publishing


def test_publishable_carries_no_restricted_field():
    """ChEMBL-derived content stays in the enrichment zone (D7); openFDA fields may go."""
    identity = harmonize.Identity(key="CHEMBL9000001", name="INVENTIB",
                                  rxcui=("999001",), unii=("UNII0001",),
                                  pharm_class_moa=("Invented Inhibitor [MoA]",),
                                  drug_type="Small molecule", clinical_stage="APPROVAL")
    published = identity.publishable()
    assert set(published) == {"chembl_id", "drug_name", "drug_type", "clinical_stage",
                              "rxcui", "unii", "pharm_class_epc", "pharm_class_moa",
                              "rxcui_resolution", "n_chembl_ids_collapsed"}
    assert published["rxcui"] == "999001"


def test_an_unresolved_rxcui_is_reported_not_dropped():
    """~45% of candidates have no RxCUI; dropping them would hide non-US drugs."""
    identity = harmonize.Identity(key="CHEMBL9000001", name="INVENTIB")
    published = identity.publishable()
    assert published["rxcui"] == "" and published["rxcui_resolution"] == "unresolved"


def test_collapsed_counts_the_ids_folded_in():
    identity = harmonize.Identity(key="C1", chembl_ids=frozenset({"C1", "C2", "C3"}))
    assert identity.collapsed == 2


# ------------------------------------------------------------------- end to end


def test_build_resolves_by_unii_and_falls_back_to_name(monkeypatch, tmp_path):
    _molecules(monkeypatch, [
        {"id": "CHEMBL9000001", "name": "INVENTIB", "drugType": "Small molecule",
         "maximumClinicalStage": "APPROVAL"},
        {"id": "CHEMBL9000002", "name": "INVENTIB CITRATE", "parentId": "CHEMBL9000001"},
        {"id": "CHEMBL9000003", "name": "OTHERAZOLE"},
        {"id": "CHEMBL9000004", "name": "UNKNOWNAZOLE"},
    ])
    unichem = _unichem(tmp_path, [("CHEMBL9000001", "UNII0001")])
    openfda = _openfda(tmp_path, [
        {"generic_name": "Inventib", "openfda": {"rxcui": ["999001"], "unii": ["UNII0001"]}},
        {"generic_name": "Otherazole", "openfda": {"rxcui": ["999003"]}},
    ])
    monkeypatch.setattr(harmonize, "ensure_sources",
                        lambda config: ({"unichem": unichem, "openfda": [openfda]}, []))
    monkeypatch.setattr(harmonize.enrichment, "ensure_datasets",
                        lambda config, datasets, channel=None: ({"drug_molecule": ["p"]}, []))

    identities, stats, _ = harmonize.build(
        {}, ["CHEMBL9000001", "CHEMBL9000002", "CHEMBL9000003", "CHEMBL9000004"])

    # The formulation collapsed onto its parent and shares one identity object.
    assert identities["CHEMBL9000002"].key == "CHEMBL9000001"
    assert identities["CHEMBL9000001"].resolution == "unii"
    assert identities["CHEMBL9000003"].resolution == "name"
    assert identities["CHEMBL9000004"].resolution == "unresolved"
    assert stats["unii"] == 1 and stats["name"] == 1 and stats["unresolved"] == 1
    assert stats["identities_out"] == 3 and stats["collapsed_ids"] == 1


def test_build_returns_an_identity_for_every_input_id(monkeypatch, tmp_path):
    """Nothing is dropped for being unresolvable."""
    _molecules(monkeypatch, [{"id": f"CHEMBL900000{i}", "name": f"INV{i}"} for i in range(3)])
    unichem = _unichem(tmp_path, [])
    openfda = _openfda(tmp_path, [])
    monkeypatch.setattr(harmonize, "ensure_sources",
                        lambda config: ({"unichem": unichem, "openfda": [openfda]}, []))
    monkeypatch.setattr(harmonize.enrichment, "ensure_datasets",
                        lambda config, datasets, channel=None: ({"drug_molecule": ["p"]}, []))
    ids = [f"CHEMBL900000{i}" for i in range(3)]
    identities, stats, _ = harmonize.build({}, ids)
    assert set(identities) == set(ids)
    assert stats["unresolved"] == 3
