"""Tests for the L0 gnomAD population-frequency lookup.

Every input is invented: gene symbols (GENE1, GENE2), coordinates, alleles and counts. No
patient variant appears here, and no real gnomAD site does either.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.l0_genomics import annotate, causal, gnomad

REPO_ROOT = Path(__file__).resolve().parent.parent

REGIONS = {"GENE1": ("1", 1000, 2000), "GENE2": ("2", 5000, 6000)}
PAD = 100


def _row(contig="chr1", pos=1500, ref="C", alt="T", filt="PASS", ac="3", an="1000",
         af="0.003", nhomalt="0", grpmax="0.005") -> str:
    return "\t".join((contig, str(pos), ref, alt, filt, ac, an, af, nhomalt, grpmax)) + "\n"


def _allele(vid="v1", gene="GENE1", contig="1", pos=1500, ref="C", alt="T") -> causal.Allele:
    return causal.Allele(variant_id=vid, gene=gene, contig=contig, pos=pos, ref=ref, alt=alt,
                         zygosity="het", filter="PASS")


# ------------------------------------------------------------------ parsing


def test_contig_names_normalise_both_ways() -> None:
    assert gnomad.bare_contig("chr15") == "15" == gnomad.bare_contig("15")
    assert gnomad.remote_contig("15") == "chr15" == gnomad.remote_contig("chr15")


def test_parse_sites_reads_every_field() -> None:
    (site,) = gnomad.parse_sites([_row()])
    assert (site.contig, site.pos, site.ref, site.alt, site.filter) == ("chr1", 1500, "C", "T", "PASS")
    assert (site.ac, site.an, site.af, site.nhomalt, site.af_grpmax) == (3, 1000, 0.003, 0, 0.005)


def test_parse_sites_keeps_missing_values_missing() -> None:
    """bcftools prints '.' for an absent field; a guessed 0 would read as "never seen"."""
    (site,) = gnomad.parse_sites([_row(nhomalt=".", grpmax=".")])
    assert site.nhomalt is None and site.af_grpmax is None


def test_parse_sites_splits_a_multiallelic_row() -> None:
    sites = gnomad.parse_sites([_row(alt="T,G", ac="3,1", af="0.003,0.001", nhomalt="1,0",
                                     grpmax="0.004,0.002")])
    assert [(s.alt, s.ac, s.af, s.nhomalt, s.an) for s in sites] == [
        ("T", 3, 0.003, 1, 1000), ("G", 1, 0.001, 0, 1000)]


def test_parse_sites_skips_blank_lines() -> None:
    assert len(gnomad.parse_sites(["\n", _row(), ""])) == 1


def test_site_key_is_minimal_and_contig_normalised() -> None:
    (site,) = gnomad.parse_sites([_row(ref="CAG", alt="CG")])
    allele = _allele(ref="CA", alt="C")
    assert site.key == ("1", 1500, "CA", "C") == ("1", *gnomad.minimal(allele.pos, allele.ref, allele.alt))


# ------------------------------------------------------------------- lookup


def test_lookup_reports_an_observed_allele() -> None:
    sites = gnomad.parse_sites([_row(), _row(pos=1510, ac="1", an="990")])
    entry = gnomad.lookup(_allele(), sites)
    assert entry["observed"] is True
    assert (entry["ac"], entry["an"], entry["nhomalt"], entry["filter"]) == (3, 1000, 0, "PASS")


def test_lookup_does_not_match_another_allele_at_the_same_position() -> None:
    entry = gnomad.lookup(_allele(alt="G"), gnomad.parse_sites([_row(alt="T")]))
    assert entry["observed"] is False
    assert "ac" not in entry


def test_absence_carries_how_well_the_position_was_sampled() -> None:
    """An allele absent where the position was barely sampled is unmeasured, not rare."""
    sites = gnomad.parse_sites([
        _row(pos=1490, an="400"), _row(pos=1520, an="600"),   # nearby, poorly sampled
        _row(pos=1900, an="1000"),                            # the gene's best-sampled site
    ])
    entry = gnomad.lookup(_allele(alt="A"), sites)
    assert entry["observed"] is False
    assert entry["n_sites_nearby"] == 2
    assert entry["an_nearby_median"] == 500
    assert entry["an_max_in_gene"] == 1000
    assert entry["an_nearby_fraction_of_max"] == 0.5


def test_no_nearby_site_gives_no_proxy_rather_than_a_guess() -> None:
    entry = gnomad.lookup(_allele(), gnomad.parse_sites([_row(pos=1900, alt="G")]))
    assert entry["n_sites_nearby"] == 0
    assert entry["an_nearby_median"] is None
    assert entry["an_nearby_fraction_of_max"] is None


def test_load_sites_refuses_an_empty_extract(tmp_path) -> None:
    """Every panel gene has gnomAD sites; zero means a broken read, not a clean gene."""
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="holds no site"):
        gnomad.load_sites(path, "GENE1", "exomes")


# ------------------------------------------------ the guardrail (decision D8)


def _record_extracts(monkeypatch):
    calls = []

    def fake_extract(config, dataset, gene, span):
        calls.append((dataset, gene, span))
        return Path("unused"), {"url": "u", "region": "r"}

    monkeypatch.setattr(gnomad, "extract", fake_extract)
    monkeypatch.setattr(gnomad, "load_sites", lambda path, gene, dataset: gnomad.parse_sites([_row()]))
    return calls


def test_the_queried_regions_never_depend_on_the_alleles(monkeypatch) -> None:
    """The remote read must be the same whatever the subject carries.

    If an allele's coordinates ever shaped the request, a server log would reveal them --
    the thing decision D8 exists to prevent.
    """
    calls = _record_extracts(monkeypatch)
    gnomad.frequencies({"seed": 42}, [_allele()], REGIONS, pad=PAD)
    first = list(calls)
    calls.clear()
    gnomad.frequencies({"seed": 42}, [_allele(pos=1234, alt="G"), _allele("v2", "GENE2", "2", 5500)],
                       REGIONS, pad=PAD)
    assert calls == first
    assert sorted({span for _, _, span in first}) == [("1", 900, 2100), ("2", 4900, 6100)]


def test_extract_cannot_see_an_allele() -> None:
    import inspect

    assert list(inspect.signature(gnomad.extract).parameters) == ["config", "dataset", "gene", "span"]


def test_only_core_fields_are_read() -> None:
    """SpliceAI annotations in the same files are CC BY-NC; nothing but counts is parsed."""
    assert gnomad.FIELDS == ("AC", "AN", "AF", "nhomalt", "AF_grpmax")
    source = Path(gnomad.__file__).read_text(encoding="utf-8")
    for forbidden in ("%INFO/spliceai", "%INFO/cadd", "%INFO/revel", "%INFO/pangolin",
                      "%INFO/phylop", "%INFO/vep", "%INFO/vrs"):
        assert forbidden not in source.lower()


def _stub_bcftools(monkeypatch, tmp_path, *, returncode=0, rows=(_row(),)):
    monkeypatch.setenv("MVA_REF_ROOT", str(tmp_path / "ref"))
    monkeypatch.setattr(gnomad, "_remote_headers",
                        lambda url: {"last_modified": "then", "etag": "e", "bytes": 1})
    seen = []

    def fake_run(argv, stdout, stderr, text, cwd, timeout, check):
        seen.append({"argv": argv, "cwd": Path(cwd)})
        stdout.writelines(rows)
        (Path(cwd) / "remote.vcf.bgz.tbi").write_text("index", encoding="utf-8")

        class Result:
            pass

        result = Result()
        result.returncode, result.stderr = returncode, "network unreachable"
        return result

    monkeypatch.setattr(gnomad.subprocess, "run", fake_run)
    return seen


def test_extract_queries_the_gene_span_and_caches_it(tmp_path, monkeypatch) -> None:
    seen = _stub_bcftools(monkeypatch, tmp_path)
    path, provenance = gnomad.extract({}, "exomes", "GENE1", ("1", 900, 2100))

    (call,) = seen
    assert call["argv"][:4] == ["bcftools", "query", "-r", "chr1:900-2100"]
    assert call["argv"][-1] == gnomad.URL_TEMPLATE.format(
        release=gnomad.RELEASE, dataset="exomes", contig="chr1")
    assert call["cwd"] == path.parent, "htslib's downloaded index must land in the cache"
    assert not list(path.parent.glob("*.tbi"))
    assert provenance["region"] == "chr1:900-2100"
    assert provenance["remote"]["etag"] == "e"
    assert len(provenance["sha256"]) == 64

    gnomad.extract({}, "exomes", "GENE1", ("1", 900, 2100))
    assert len(seen) == 1, "a cached extract must not be read again"


def test_a_failed_read_is_an_error_and_leaves_nothing_behind(tmp_path, monkeypatch) -> None:
    """A network failure must never read as "gnomAD has no variant in this gene"."""
    _stub_bcftools(monkeypatch, tmp_path, returncode=1, rows=())
    with pytest.raises(RuntimeError, match="network unreachable"):
        gnomad.extract({}, "exomes", "GENE1", ("1", 900, 2100))
    cache = gnomad.gnomad_dir({})
    assert not list(cache.iterdir())


def test_a_timed_out_read_leaves_no_partial_table(tmp_path, monkeypatch) -> None:
    """A partial table left behind would be read as a complete gene next time."""
    _stub_bcftools(monkeypatch, tmp_path)

    def slow(argv, stdout, stderr, text, cwd, timeout, check):
        stdout.write(_row())
        (Path(cwd) / "remote.vcf.bgz.tbi").write_text("index", encoding="utf-8")
        raise gnomad.subprocess.TimeoutExpired(argv, timeout)

    monkeypatch.setattr(gnomad.subprocess, "run", slow)
    with pytest.raises(gnomad.subprocess.TimeoutExpired):
        gnomad.extract({}, "exomes", "GENE1", ("1", 900, 2100))
    assert not list(gnomad.gnomad_dir({}).iterdir())


def test_the_cache_is_outside_the_repo(monkeypatch) -> None:
    monkeypatch.delenv("MVA_REF_ROOT", raising=False)
    assert REPO_ROOT not in gnomad.gnomad_dir({}).parents


# ------------------------------------------------------------ the artifact


def test_frequencies_assembles_the_artifact(monkeypatch) -> None:
    _record_extracts(monkeypatch)
    result = gnomad.frequencies({"seed": 42, "gnomad": {"datasets": ["genomes"]}},
                                [_allele(), _allele("v2", alt="G")], REGIONS, pad=PAD)
    assert result["status"] == "ok"
    assert result["release"]["datasets"] == ["genomes"]
    assert result["release"]["extracts"]["genomes"]["GENE1"]["n_sites"] == 1
    assert [a["genomes"]["observed"] for a in result["alleles"]] == [True, False]
    assert "exomes" not in result["alleles"][0]
    assert result["caveats"] == list(gnomad.CAVEATS)
    assert result["seed"] == 42


def test_summarize_keeps_the_numbers_a_reviewer_needs() -> None:
    entry = {"variant_id": "v1",
             "exomes": {"observed": True, "site_in_release": True, "ac": 1, "an": 1000,
                        "nhomalt": 0, "an_nearby_fraction_of_max": 1.0, "af": 0.001,
                        "filter": "PASS"},
             "genomes": {"observed": False, "site_in_release": False,
                         "an_nearby_fraction_of_max": 0.9}}
    out = gnomad.summarize(entry)
    assert out["exomes"] == {"observed": True, "site_in_release": True, "ac": 1, "an": 1000,
                             "af": 0.001, "nhomalt": 0, "filter": "PASS",
                             "an_nearby_fraction_of_max": 1.0}
    assert out["genomes"]["observed"] is False and out["genomes"]["ac"] is None


def test_enabled_defaults_on_and_can_be_turned_off() -> None:
    assert gnomad.enabled({}) is True
    assert gnomad.enabled({"gnomad": {"enabled": False}}) is False


def test_the_configured_release_is_the_documented_one() -> None:
    import yaml

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    assert config["gnomad"]["release"] == gnomad.RELEASE
    assert set(config["gnomad"]["datasets"]) <= set(gnomad.DATASETS)


# -------------------------------------------------- integration with causal.call


CFG = {"annotator": {"tool": "snpeff", "database": "DBX", "java_heap": "8g",
                     "transcript_policy": "mane_select_tiered", "mane_release": "MANE.test"},
       "seed": 42, "clinvar": {"enabled": False}}


def _stub_causal(monkeypatch, frequencies):
    alleles = [_allele("v1", pos=1500, alt="T"), _allele("v2", pos=1700, ref="G", alt="A")]
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
    monkeypatch.setattr(causal._gnomad, "frequencies", frequencies)


def _fake_frequencies(config, alleles, regions, pad):
    return {"status": "ok", "release": {"version": "4.1.1", "datasets": ["exomes"]},
            "alleles": [{"variant_id": a.variant_id,
                         "exomes": {"observed": a.variant_id == "v1", "ac": 5, "an": 1000,
                                    "nhomalt": 0, "an_nearby_fraction_of_max": 1.0}}
                        for a in alleles]}


def test_the_g1_artifact_carries_the_frequencies(tmp_path, monkeypatch) -> None:
    _stub_causal(monkeypatch, _fake_frequencies)
    summary = causal.call("unused.vcf.gz", CFG, tmp_path, panel=("GENE1",))

    assert summary["gnomad"] == "ok"
    written = json.loads((tmp_path / "causal_gene_call.json").read_text(encoding="utf-8"))
    assert written["gnomad"] == {"status": "ok", "release": "4.1.1",
                                 "artifact": "gnomad_frequencies.json"}
    attached = written["per_gene"]["GENE1"]["configurations"][0]["gnomad"]
    assert [a["exomes"]["observed"] for a in attached] == [True, False]
    assert causal.FREQUENCY_CAVEATS["ok"] in written["caveats"]
    assert causal.FREQUENCY_CAVEATS["disabled"] not in written["caveats"]
    assert (tmp_path / "gnomad_frequencies.json").exists()


def test_disabling_the_lookup_says_frequency_was_not_assessed(tmp_path, monkeypatch) -> None:
    """"Not assessed" and "rare" must never look alike in the G1 artifact."""
    def explode(*args, **kwargs):
        raise AssertionError("frequencies must not run when disabled")

    _stub_causal(monkeypatch, explode)
    causal.call("unused.vcf.gz", dict(CFG, gnomad={"enabled": False}), tmp_path, panel=("GENE1",))

    written = json.loads((tmp_path / "causal_gene_call.json").read_text(encoding="utf-8"))
    assert written["gnomad"]["status"] == "disabled"
    assert causal.FREQUENCY_CAVEATS["disabled"] in written["caveats"]
    assert "gnomad" not in written["per_gene"]["GENE1"]["configurations"][0]
    assert not (tmp_path / "gnomad_frequencies.json").exists()


def test_a_failed_lookup_keeps_the_annotation_but_writes_no_g1_artifact(tmp_path, monkeypatch) -> None:
    def boom(*args, **kwargs):
        raise RuntimeError("bcftools could not read gnomAD")

    _stub_causal(monkeypatch, boom)
    with pytest.raises(RuntimeError):
        causal.call("unused.vcf.gz", CFG, tmp_path, panel=("GENE1",))
    assert (tmp_path / "variant_calls.json").exists()
    assert not (tmp_path / "causal_gene_call.json").exists()


def test_every_caveat_is_kept_whichever_way_the_lookup_went(tmp_path, monkeypatch) -> None:
    _stub_causal(monkeypatch, _fake_frequencies)
    causal.call("unused.vcf.gz", CFG, tmp_path, panel=("GENE1",))
    written = json.loads((tmp_path / "causal_gene_call.json").read_text(encoding="utf-8"))
    assert len(written["caveats"]) == len(causal.CAVEATS) + 1
    for caveat in causal.CAVEATS:
        assert caveat in written["caveats"]


def test_a_site_with_no_passing_carrier_is_not_observed() -> None:
    """gnomAD keeps sites whose carriers all failed its filters (AC=0).

    Reporting those as observed would read, at gate G1, as "the population carries this".
    """
    sites = gnomad.parse_sites([_row(ac="0", af="0", filt="AC0"), _row(pos=1510, an="900")])
    entry = gnomad.lookup(_allele(), sites)
    assert entry["observed"] is False
    assert entry["site_in_release"] is True
    assert (entry["ac"], entry["filter"]) == (0, "AC0")


def test_an_allele_absent_from_the_release_says_so_too() -> None:
    entry = gnomad.lookup(_allele(alt="A"), gnomad.parse_sites([_row()]))
    assert entry["observed"] is False and entry["site_in_release"] is False


def test_the_summary_carries_the_filter_and_frequency() -> None:
    """A flagged site and a clean one are not the same evidence for a G1 reviewer."""
    entry = {"variant_id": "v1",
             "exomes": {"observed": True, "site_in_release": True, "ac": 3, "an": 1000,
                        "af": 0.003, "nhomalt": 0, "filter": "PASS",
                        "an_nearby_fraction_of_max": 1.0}}
    out = gnomad.summarize(entry, ("exomes",))
    assert out["exomes"]["filter"] == "PASS" and out["exomes"]["af"] == 0.003
    assert out["exomes"]["site_in_release"] is True
