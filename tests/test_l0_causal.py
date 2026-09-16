"""Tests for the L0 causal-gene call.

Every input is invented: gene symbols (GENE1, GENE2), coordinates, alleles, phasing and
annotations. No patient variant, coordinate, genotype or phenotype appears here -- and no
real pathogenic variant does either, because a known allele in a committed test could be
read as a statement about the subject.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.l0_genomics import annotate, causal

CFG = {"annotator": {"tool": "snpeff", "database": "DBX", "java_heap": "8g",
                     "transcript_policy": "mane_select_tiered", "mane_release": "MANE.test"}}


def _ann(effects="stop_gained", impact="HIGH", gene="GENE1", transcript="ENST00000000001.1",
         biotype="protein_coding", rank="3/10") -> annotate.Annotation:
    return annotate.Annotation(allele="T", effects=tuple(effects.split("&")), impact=impact,
                               gene=gene, feature_type="transcript", transcript=transcript,
                               biotype=biotype, rank=rank, hgvs_c="c.10C>T", hgvs_p="")


def _allele(vid="v1", zygosity="het", filter="PASS", phase_gt="", phase_id="",
            trans_partner="", gene="GENE1") -> causal.Allele:
    return causal.Allele(variant_id=vid, gene=gene, contig="1", pos=100, ref="C", alt="T",
                         zygosity=zygosity, filter=filter, phase_gt=phase_gt,
                         phase_id=phase_id, trans_partner=trans_partner)


def _classify(allele, primary=(), extra=(), nmd=("GENE1",)):
    """Primary pass sees only ``primary``; the all-transcript pass sees it plus ``extra``."""
    genes = frozenset(nmd)
    return causal.classify(
        allele,
        {allele.variant_id: (tuple(primary), genes)},
        {allele.variant_id: (tuple(primary) + tuple(extra), genes)},
        policy="mane_select_tiered", mane_release="MANE.test",
    )


def _pair(allele, primary=(), extra=()):
    return allele, _classify(allele, primary, extra)


# ------------------------------------------------------------------------- snpEff


def test_snpeff_argv_carries_every_privacy_flag() -> None:
    argv = annotate.snpeff_argv(CFG, "ann", "-lof", "DBX", "-")
    for flag in ("-noLog", "-nodownload", "-noStats", "-Xmx8g"):
        assert flag in argv
    bed = annotate.snpeff_argv(CFG, "genes2bed", "DBX", "GENE1")
    assert "-noLog" in bed and "-nodownload" in bed


def test_parse_ann_and_nmd_tag() -> None:
    info = (
        "AC=1;ANN=T|stop_gained&splice_region_variant|HIGH|GENE1|ENSG00000000001|transcript|"
        "ENST00000000001.1|protein_coding|3/10|c.10C>T|p.Gln4*|10/900|10/800|4/260||,"
        "T|intron_variant|MODIFIER|GENE1|ENSG00000000001|transcript|ENST00000000002.1|"
        "protein_coding|2/9|c.5+20C>T||||||;NMD=(GENE1|ENSG00000000001|1|1.00)"
    )
    anns = annotate.parse_ann(info)
    assert len(anns) == 2
    assert anns[0].effects == ("stop_gained", "splice_region_variant")
    assert anns[0].transcript == "ENST00000000001.1" and anns[0].rank == "3/10"
    assert annotate.parse_tag_genes(info, "NMD") == frozenset({"GENE1"})
    assert annotate.parse_tag_genes("AC=1", "NMD") == frozenset()


# --------------------------------------------------------------------- consequence


@pytest.mark.parametrize(
    ("ann", "expected"),
    [
        (_ann("stop_gained"), "lof"),
        (_ann("splice_donor_variant&intron_variant"), "lof"),
        (_ann("splice_donor_variant", biotype="retained_intron"), "other"),
        (_ann("missense_variant", "MODERATE"), "protein_altering"),
        (_ann("synonymous_variant", "LOW"), "other"),
    ],
)
def test_effect_class(ann, expected: str) -> None:
    assert causal.effect_class(ann) == expected


def test_lof_on_mane_select_is_the_call() -> None:
    call = _classify(_allele(), [_ann("stop_gained")])
    assert (call.lof_tier, call.annotation_pass, call.effect_class) == ("mane_select", "primary", "lof")


def test_lof_only_on_another_transcript_is_non_mane() -> None:
    call = _classify(_allele(), [_ann("missense_variant", "MODERATE")],
                     [_ann("stop_gained", transcript="ENST00000000002.1")])
    assert call.lof_tier == "non_mane"
    assert call.annotation_pass == "secondary"
    assert call.transcript == "ENST00000000002.1"


def test_missense_is_not_lof_but_is_kept() -> None:
    call = _classify(_allele(), [_ann("missense_variant", "MODERATE")])
    assert (call.lof_tier, call.effect_class) == ("not_lof", "protein_altering")


def test_allele_outside_the_genes_transcripts_is_uncovered() -> None:
    assert _classify(_allele(), [_ann("stop_gained", gene="GENE2")]) is None


def test_last_exon_truncation_is_flagged_not_demoted() -> None:
    call = _classify(_allele(), [_ann("stop_gained", rank="10/10")], nmd=())
    assert call.lof_tier == "mane_select"
    assert call.last_exon is True
    assert "nmd_escape_last_exon" in call.flags
    assert "snpeff_predicts_no_nmd" in call.flags


def test_intronic_rank_is_not_read_as_an_exon() -> None:
    call = _classify(_allele(), [_ann("splice_donor_variant&intron_variant", rank="3/9")])
    assert call.exon_rank == "" and call.last_exon is None
    assert not any(f.startswith("nmd") or "nmd" in f for f in call.flags)


def test_filtered_allele_is_flagged() -> None:
    call = _classify(_allele(filter="FS60"), [_ann("stop_gained")])
    assert "filtered" in call.flags


# ---------------------------------------------------------------------------- phase


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        (_allele("v1", trans_partner="v2"), _allele("v2", trans_partner="v1"), "trans"),
        (_allele("v1", phase_gt="0|1", phase_id="p1"), _allele("v2", phase_gt="1|0", phase_id="p1"), "trans"),
        (_allele("v1", phase_gt="0|1", phase_id="p1"), _allele("v2", phase_gt="0|1", phase_id="p1"), "cis"),
        (_allele("v1", phase_gt="0|1", phase_id="p1"), _allele("v2", phase_gt="1|0", phase_id="p2"), "unphased"),
        (_allele("v1"), _allele("v2"), "unphased"),
    ],
)
def test_phase_relation(a, b, expected: str) -> None:
    assert causal.phase_relation(a, b) == expected


# -------------------------------------------------------------------- configuration

LOF = [_ann("stop_gained")]
MISSENSE = [_ann("missense_variant", "MODERATE")]


def test_homozygous_lof_is_biallelic() -> None:
    result = causal.resolve_gene([_pair(_allele(zygosity="hom_alt"), LOF)])
    assert result["verdict"] == "biallelic_lof"
    assert result["configurations"][0]["phase"] == "homozygous"


def test_two_het_lof_in_trans_is_biallelic() -> None:
    a = _allele("v1", phase_gt="0|1", phase_id="p1")
    b = _allele("v2", phase_gt="1|0", phase_id="p1")
    result = causal.resolve_gene([_pair(a, LOF), _pair(b, LOF)])
    assert result["verdict"] == "biallelic_lof"
    assert result["configurations"][0]["phase"] == "trans"


def test_two_het_lof_in_cis_is_not_biallelic() -> None:
    a = _allele("v1", phase_gt="0|1", phase_id="p1")
    b = _allele("v2", phase_gt="0|1", phase_id="p1")
    result = causal.resolve_gene([_pair(a, LOF), _pair(b, LOF)])
    assert result["verdict"] == "monoallelic_lof"
    assert result["configurations"] == []


def test_unphased_pair_is_a_candidate_with_a_flag() -> None:
    result = causal.resolve_gene([_pair(_allele("v1"), LOF), _pair(_allele("v2"), LOF)])
    assert result["verdict"] == "biallelic_lof"
    assert "unphased_pair" in result["configurations"][0]["flags"]


def test_pair_with_a_non_mane_lof_is_demoted() -> None:
    other = [_ann("stop_gained", transcript="ENST00000000002.1")]
    result = causal.resolve_gene([_pair(_allele("v1"), LOF), _pair(_allele("v2"), MISSENSE, other)])
    assert result["verdict"] == "biallelic_lof_non_mane"


def test_lof_plus_missense_is_its_own_class() -> None:
    result = causal.resolve_gene([_pair(_allele("v1"), LOF), _pair(_allele("v2"), MISSENSE)])
    assert result["verdict"] == "lof_plus_protein_altering"


def test_filtered_lof_never_enters_a_configuration_but_is_listed() -> None:
    result = causal.resolve_gene([_pair(_allele("v1", zygosity="hom_alt", filter="QD2"), LOF)])
    assert result["verdict"] == "no_lof"
    assert result["filtered_lof_variant_ids"] == ["v1"]


def test_single_het_lof_is_monoallelic() -> None:
    assert causal.resolve_gene([_pair(_allele(), LOF)])["verdict"] == "monoallelic_lof"


def test_no_alleles_is_no_lof() -> None:
    assert causal.resolve_gene([])["verdict"] == "no_lof"


# -------------------------------------------------------------------------- verdict


def _gene(verdict: str) -> dict:
    return {"verdict": verdict}


def test_overall_single_and_multiple_candidates() -> None:
    assert causal.overall_verdict({"GENE1": _gene("biallelic_lof"), "GENE2": _gene("no_lof")}) == (
        "single_candidate", ["GENE1"])
    assert causal.overall_verdict({"GENE1": _gene("biallelic_lof"), "GENE2": _gene("biallelic_lof")})[0] == (
        "multiple_candidates")


def test_non_mane_only_is_not_a_candidate_call() -> None:
    verdict, genes = causal.overall_verdict({"GENE1": _gene("biallelic_lof_non_mane")})
    assert verdict == "non_mane_candidates_only" and genes == ["GENE1"]


def test_nothing_found_is_never_worded_as_no_causal_variant() -> None:
    verdict, genes = causal.overall_verdict({"GENE1": _gene("monoallelic_lof")})
    assert (verdict, genes) == ("no_biallelic_lof", [])
    note = causal.VERDICT_NOTES[verdict]
    assert "not evidence that no causal variant exists" in note


# ------------------------------------------------------------------- reading the VCF


def _write_vcf(tmp_path: Path, rows: list) -> Path:
    """A synthetic, invented single-sample VCF in pytest's temp dir -- never in the repo."""
    pysam = pytest.importorskip("pysam")
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=100000>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=FS60,Description="invented">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype">',
        '##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    raw = tmp_path / "synthetic.vcf"
    raw.write_text("\n".join(header + rows) + "\n")
    gz = str(tmp_path / "synthetic.vcf.gz")
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)
    return Path(gz)


def test_read_alleles_zygosity_phase_filter_and_region(tmp_path: Path) -> None:
    rows = [
        "1\t1100\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ:PGT:PID\t0/1:30:99:0|1:1100_C_T",
        "1\t1200\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t1/1:30:99",
        "1\t1300\t.\tA\tT,G\t50\tPASS\t.\tGT:DP:GQ\t1/2:30:99",
        "1\t1400\t.\tA\tC\t50\tPASS\t.\tGT:DP:GQ\t0/0:30:99",
        "1\t1500\t.\tA\tC\t50\tPASS\t.\tGT:DP:GQ\t./.:0:0",
        "1\t1600\t.\tT\tC\t50\tFS60\t.\tGT:DP:GQ\t0/1:30:99",
        "1\t1700\t.\tT\tC,*\t50\tPASS\t.\tGT:DP:GQ\t1/2:30:99",
        "1\t9000\t.\tT\tC\t50\tPASS\t.\tGT:DP:GQ\t1/1:30:99",
    ]
    alleles = causal.read_alleles(_write_vcf(tmp_path, rows), {"GENE1": ("1", 1000, 2000)}, pad=100)
    by_pos = {}
    for a in alleles:
        by_pos.setdefault(a.pos, []).append(a)

    assert sorted(by_pos) == [1100, 1200, 1300, 1600, 1700]      # 0/0, ./. and 9000 dropped
    assert by_pos[1100][0].zygosity == "het"
    assert (by_pos[1100][0].phase_gt, by_pos[1100][0].phase_id) == ("0|1", "1100_C_T")
    assert by_pos[1200][0].zygosity == "hom_alt"
    t, g = by_pos[1300]
    assert (t.zygosity, g.zygosity) == ("het", "het")
    assert (t.trans_partner, g.trans_partner) == (g.variant_id, t.variant_id)
    assert by_pos[1600][0].filter == "FS60"
    assert [a.alt for a in by_pos[1700]] == ["C"]                 # spanning deletion skipped
    assert len({a.variant_id for a in alleles}) == len(alleles)


def test_read_alleles_rejects_a_contig_mismatch(tmp_path: Path) -> None:
    vcf = _write_vcf(tmp_path, ["1\t1100\t.\tC\tT\t50\tPASS\t.\tGT\t0/1"])
    with pytest.raises(ValueError, match="prefix"):
        causal.read_alleles(vcf, {"GENE1": ("chr1", 1000, 2000)})
