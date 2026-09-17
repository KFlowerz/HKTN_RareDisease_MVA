"""L0 -- causal-gene call: panel variants, consequence, zygosity, configuration.

Purpose
    Find which panel gene, if any, carries a biallelic loss-of-function (LoF)
    configuration in the subject's VCF, and ship that candidate call with every variant
    that was considered and every one that was rejected -- so a reader can disagree with
    it rather than take it on trust.

Inputs
    The single-sample VCF (bgzipped, tabix-indexed); ``config["annotator"]``; the gene
    panel; optional ``config["l0_region_pad"]`` (bp around each gene span).

Outputs
    Under ``results_dir/l0_genomics/``:
      - ``variant_calls.json`` -- each panel allele's genotype record and its
        :class:`~src.l0_genomics.transcripts.VariantCall`; alleles that no panel
        transcript covers are listed separately
      - ``causal_gene_call.json`` -- per-gene configurations and verdicts, the overall
        verdict, provenance, and caveats; each configuration carries what ClinVar holds
        about its alleles
      - ``clinvar_crossref.json`` -- the full cross-reference behind that summary
        (:mod:`src.l0_genomics.clinvar`)
      - ``gnomad_frequencies.json`` -- each allele's population frequency
        (:mod:`src.l0_genomics.gnomad`)

Method
    1. Regions: each panel gene's span from the local snpEff database, padded.
    2. Alleles: every non-reference allele the sample carries there, with zygosity,
       FILTER status, and GATK physical phasing (``PGT``/``PID``) where the caller
       emitted it.
    3. Consequence: the two decision-D5 passes -- MANE Select, then all transcripts --
       give each allele a LoF tier. LoF means a nonsense, frameshift, canonical
       splice-site or exon-loss consequence [macarthur2012] doi:10.1126/science.1215040,
       on a protein-coding transcript.
    4. Configuration: MVA is caused by biallelic mutation [hanks2004] doi:10.1038/ng1449,
       so one heterozygous LoF is not a call. Two heterozygous alleles are cis, trans or
       unphased; a cis pair leaves the other copy intact and is never biallelic.
    5. Cross-reference: each allele is matched against the cached public ClinVar release
       [landrum2018] doi:10.1093/nar/gkx1153, so a configuration is weighed against what
       has already been submitted about its alleles and not against prediction alone.
       See :mod:`src.l0_genomics.clinvar`; it interprets nothing.
    6. Population frequency: each allele is looked up in gnomAD [chen2024]
       doi:10.1038/s41586-023-06045-0, read by whole panel-gene span so no subject
       coordinate leaves the machine (decision D8). See :mod:`src.l0_genomics.gnomad`; it
       interprets nothing either.

Guardrail
    - The output is a *candidate* for gate G1. Nothing here writes
      ``config["causal_gene"]``; a person reads the evidence and sets it.
    - "No biallelic LoF found" is reported as exactly that, never as "no causal variant".
      Missense, structural, deep-intronic and mosaic second hits are not assessed, and 6
      of the 118 known pathogenic panel variants are not LoF under any transcript policy
      (``results/transcript_policy/clinvar_summary.tsv``, seed=42).
    - Only FILTER=PASS alleles enter a configuration. Filtered alleles stay in the output
      with the failing filter named: in an n=1 search they deserve a look, not a silent
      drop.
    - Per-variant records identify the child. They are written under ``results_dir`` only
      and never logged; the return value and logs carry counts and categories.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

from . import annotate
from . import clinvar as _clinvar
from . import gnomad as _gnomad
from . import transcripts as tx

#: Padding around each gene span. Covers either BED coordinate convention and variants
#: just outside the outermost transcript; consequences are still read only from the
#: panel gene's own transcripts.
DEFAULT_PAD = 1000

#: LoF consequences [macarthur2012] doi:10.1126/science.1215040 -- nonsense, frameshift,
#: canonical splice site, and deletions removing exons.
LOF_EFFECTS = frozenset({
    "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
    "exon_loss_variant", "transcript_ablation",
})
#: Truncating LoF, whose effect depends on nonsense-mediated decay: a premature stop in
#: the last exon is not expected to trigger NMD [aboutayoun2018] doi:10.1002/humu.23626.
TRUNCATING = frozenset({"stop_gained", "frameshift_variant"})
#: Consequences whose snpEff rank counts introns, not exons.
INTRONIC = frozenset({
    "intron_variant", "splice_donor_variant", "splice_acceptor_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
})
_CLASS_RANK = {"lof": 2, "protein_altering": 1, "other": 0}

#: Per-gene verdicts, strongest first.
GENE_VERDICTS = (
    "biallelic_lof", "biallelic_lof_non_mane", "lof_plus_protein_altering",
    "monoallelic_lof", "no_lof",
)

VERDICT_NOTES = {
    "single_candidate": (
        "One panel gene carries a biallelic LoF configuration on its MANE Select "
        "transcript. It is a candidate for gate G1, to be reviewed with its evidence and "
        "caveats -- not a setting."
    ),
    "multiple_candidates": (
        "More than one panel gene carries a biallelic LoF configuration on MANE Select. "
        "L0 does not choose between them; review each."
    ),
    "non_mane_candidates_only": (
        "Biallelic LoF configurations exist only when LoF on non-MANE transcripts is "
        "admitted (decision D5 tier non_mane). None can be a call without evidence that "
        "the transcript is biologically relevant."
    ),
    "lof_plus_protein_altering_only": (
        "No biallelic LoF configuration. At least one gene pairs a LoF allele with a "
        "protein-altering second allele, which is not interpreted here. Any ClinVar and "
        "population-frequency evidence attached to the configuration is reported, not "
        "weighed."
    ),
    "no_biallelic_lof": (
        "No biallelic LoF configuration was found in the panel. This is not evidence "
        "that no causal variant exists: missense, structural, deep-intronic and mosaic "
        "second hits are not assessed."
    ),
}

#: Swapped into the caveats according to whether the gnomAD lookup ran.
FREQUENCY_CAVEATS = {
    "ok": "Population allele frequency is reported from gnomAD (gnomad_frequencies.json) "
          "but not interpreted: rarity is necessary, not sufficient, for pathogenicity -- "
          "see the caveats in that file.",
    "disabled": "Population allele frequency is not assessed -- the gnomAD lookup is "
                "disabled in config -- so a protein-altering second hit may be a common "
                "benign variant.",
}

CAVEATS = (
    "Copy number is not assessed: an apparently homozygous allele is not distinguished "
    "from a heterozygous allele over a deletion of the other copy.",
    "Structural variants, deep-intronic and regulatory variants, and mosaic second hits "
    "are not assessed from a called VCF.",
    "Unphased pairs are reported as candidate biallelic configurations. Without parental "
    "samples or physical phasing, a cis arrangement cannot be excluded.",
    "Missense and other protein-altering variants are listed but not interpreted; no "
    "pathogenicity predictor is applied.",
    "A ClinVar match reports what submitters have concluded, weighted by review status. "
    "It is not a classification made here, and absence from ClinVar is not evidence of "
    "benignity -- see the caveats in clinvar_crossref.json.",
)


@dataclass(frozen=True)
class Allele:
    """One non-reference allele the sample carries in a panel region."""

    variant_id: str          # internal id ("v1", "v2", ...) -- never the dataset's own ID
    gene: str                # the panel gene whose region yielded it
    contig: str
    pos: int
    ref: str
    alt: str
    zygosity: str            # "het" | "hom_alt" | "hemizygous"
    filter: str              # "PASS", or the failing filter names
    phase_gt: str = ""       # FORMAT/PGT, e.g. "0|1"; "" when not physically phased
    phase_id: str = ""       # FORMAT/PID; physically phased records share one
    trans_partner: str = ""  # the other ALT at a 1/2 site -- trans by genotype
    dp: int | None = None
    gq: int | None = None


def _text(value) -> str:
    return "" if value is None else str(value)


def _int(value):
    try:
        return None if value is None else int(value)
    except (TypeError, ValueError):
        return None


def read_alleles(vcf_path, regions: dict, *, pad: int = DEFAULT_PAD) -> list:
    """Every non-reference allele the sample carries in the padded panel regions.

    Spanning-deletion (``*``) and symbolic (``<...>``) ALTs are skipped: they describe an
    overlapping event recorded elsewhere, not a sequence change at this position.

    Raises:
        ValueError: If the VCF is not single-sample, or a region's contig is absent from
            its header -- usually a ``chr``-prefix mismatch, which would otherwise return
            an empty panel that reads as "no variants".
    """
    import pysam

    # htslib writes warnings to stderr that quote the VCF's path, and this dataset's
    # filename embeds a lab accession and a sequencer flowcell id, which must never
    # propagate into a log or an artifact (see run.py's guardrail). Silence them.
    pysam.set_verbosity(0)

    vcf = pysam.VariantFile(str(vcf_path))
    if len(vcf.header.samples) != 1:
        raise ValueError(f"expected exactly 1 sample, found {len(vcf.header.samples)}")
    present = set(vcf.header.contigs)

    alleles, seen = [], set()
    for gene, (contig, start, end) in sorted(regions.items()):
        if contig not in present:
            raise ValueError(
                f"{gene}: contig {contig!r} is not in the VCF header -- a `chr`-prefix mismatch?"
            )
        for rec in vcf.fetch(contig, max(0, start - pad), end + pad):
            if not rec.alts:
                continue
            sample = rec.samples[0]
            gt = tuple(sample.get("GT") or ())
            called = [g for g in gt if g is not None]
            keys = list(rec.filter.keys())
            status = "PASS" if not keys or keys == ["PASS"] else ",".join(keys)

            new = []
            for idx in sorted({g for g in called if g > 0}):
                alt = rec.alts[idx - 1]
                if alt == "*" or alt.startswith("<"):
                    continue
                key = (rec.contig, rec.pos, rec.ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                if len(gt) == 1:
                    zygosity = "hemizygous"
                elif called.count(idx) == len(gt):
                    zygosity = "hom_alt"
                else:
                    zygosity = "het"
                new.append(dict(
                    variant_id=f"v{len(alleles) + len(new) + 1}", gene=gene,
                    contig=rec.contig, pos=rec.pos, ref=rec.ref, alt=alt,
                    zygosity=zygosity, filter=status,
                    phase_gt=_text(sample.get("PGT")), phase_id=_text(sample.get("PID")),
                    dp=_int(sample.get("DP")), gq=_int(sample.get("GQ")),
                ))
            if len(new) == 2 and all(n["zygosity"] == "het" for n in new):
                new[0]["trans_partner"] = new[1]["variant_id"]
                new[1]["trans_partner"] = new[0]["variant_id"]
            alleles.extend(Allele(**n) for n in new)
    return alleles


def effect_class(ann: annotate.Annotation) -> str:
    """``lof``, ``protein_altering`` or ``other``. Only protein-coding transcripts can be LoF."""
    if ann.biotype != "protein_coding":
        return "other"
    if LOF_EFFECTS & set(ann.effects):
        return "lof"
    if ann.impact in ("HIGH", "MODERATE"):
        return "protein_altering"
    return "other"


def best_on_gene(annotations, gene: str):
    """The most severe transcript-level annotation on ``gene``, or None.

    Ties break on transcript id, so the choice is deterministic.
    """
    pool = [a for a in annotations if a.gene == gene and a.feature_type == "transcript"]
    if not pool:
        return None
    return min(pool, key=lambda a: (-_CLASS_RANK[effect_class(a)],
                                    -annotate.IMPACT_RANK.get(a.impact, -1), a.transcript))


def _exon(ann: annotate.Annotation) -> tuple:
    """``(exon_rank, is_last_exon)``; ``("", None)`` where the rank counts introns."""
    if INTRONIC & set(ann.effects) or "/" not in ann.rank:
        return "", None
    number, total = ann.rank.split("/", 1)
    try:
        return ann.rank, int(number) == int(total)
    except ValueError:
        return ann.rank, None


def classify(allele: Allele, primary: dict, secondary: dict, *, policy: str,
             mane_release: str):
    """Turn one allele's two annotation passes into a :class:`VariantCall`, or None.

    None means no transcript of the allele's panel gene covers it (padding only).
    """
    p_anns, p_nmd = primary.get(allele.variant_id, ((), frozenset()))
    s_anns, s_nmd = secondary.get(allele.variant_id, ((), frozenset()))
    p, s = best_on_gene(p_anns, allele.gene), best_on_gene(s_anns, allele.gene)
    if p is None and s is None:
        return None

    tier = tx.lof_tier(
        primary_lof=p is not None and effect_class(p) == "lof",
        secondary_lof=s is not None and effect_class(s) == "lof",
    )
    if tier == "non_mane" or p is None:
        source, which, nmd_genes = s, "secondary", s_nmd
    else:
        source, which, nmd_genes = p, "primary", p_nmd

    rank, last = _exon(source)
    flags = []
    if tier != "not_lof" and TRUNCATING & set(source.effects):
        if last:
            flags.append("nmd_escape_last_exon")
        if allele.gene not in nmd_genes:
            flags.append("snpeff_predicts_no_nmd")
    if allele.filter != "PASS":
        flags.append("filtered")

    return tx.VariantCall(
        gene=allele.gene, transcript=source.transcript, annotation_pass=which,
        lof_tier=tier, mane_release=mane_release, transcript_policy=policy,
        consequence="&".join(source.effects), impact=source.impact,
        hgvs_c=source.hgvs_c, hgvs_p=source.hgvs_p, variant_id=allele.variant_id,
        zygosity=allele.zygosity, effect_class=effect_class(source), exon_rank=rank,
        last_exon=last, flags=tuple(flags),
    )


def phase_relation(a: Allele, b: Allele) -> str:
    """``cis``, ``trans`` or ``unphased`` for two heterozygous alleles.

    Two ALTs at one site (a 1/2 genotype) are trans by definition. Otherwise GATK
    physical phasing decides: this VCF's header defines ``PID`` as connecting "records
    within a phasing group" and ``PGT`` as "how the alternate alleles are phased in
    relation to one another", so alleles sharing a PID are cis when their PGT match and
    trans when they differ. Anything else is unphased -- the dataset is single-sample,
    with no parents to phase against.
    """
    if a.trans_partner == b.variant_id:
        return "trans"
    if a.phase_id and a.phase_id == b.phase_id and a.phase_gt and b.phase_gt:
        return "cis" if a.phase_gt == b.phase_gt else "trans"
    return "unphased"


def _configuration(kind: str, members: list, phase: str) -> dict:
    flags = {f for _, c in members for f in c.flags}
    if phase == "unphased":
        flags.add("unphased_pair")
    return {
        "class": kind,
        "variant_ids": [a.variant_id for a, _ in members],
        "phase": phase,
        "lof_tiers": [c.lof_tier for _, c in members],
        "effect_classes": [c.effect_class for _, c in members],
        "flags": sorted(flags),
    }


def attach_clinvar(genes: dict, crossref: dict) -> None:
    """Attach each configuration member's ClinVar summary, in place.

    A configuration is what gate G1 is decided on, so what the archive says about its
    alleles belongs beside it -- not only in a second file a reader has to join by hand.
    """
    if crossref.get("status") != "ok":
        return
    summaries = {m["variant_id"]: _clinvar.summarize(m) for m in crossref["alleles"]}
    for result in genes.values():
        for config in result["configurations"]:
            config["clinvar"] = [summaries[v] for v in config["variant_ids"] if v in summaries]


def attach_gnomad(genes: dict, frequencies: dict) -> None:
    """Attach each configuration member's gnomAD summary, in place -- see :func:`attach_clinvar`."""
    if frequencies.get("status") != "ok":
        return
    datasets = frequencies["release"]["datasets"]
    summaries = {e["variant_id"]: _gnomad.summarize(e, datasets)
                 for e in frequencies["alleles"]}
    for result in genes.values():
        for config in result["configurations"]:
            config["gnomad"] = [summaries[v] for v in config["variant_ids"] if v in summaries]


def resolve_gene(items: list) -> dict:
    """Biallelic configurations and a verdict for one gene's ``(Allele, VariantCall)`` pairs."""
    passing = [(a, c) for a, c in items if a.filter == "PASS"]
    lof = [(a, c) for a, c in passing if c.lof_tier != "not_lof"]
    altering = [(a, c) for a, c in passing
                if c.lof_tier == "not_lof" and c.effect_class == "protein_altering"]
    het_lof = [(a, c) for a, c in lof if a.zygosity == "het"]

    configs = []
    for a, c in lof:
        if a.zygosity == "hom_alt":
            kind = "biallelic_lof" if c.lof_tier == "mane_select" else "biallelic_lof_non_mane"
            configs.append(_configuration(kind, [(a, c)], "homozygous"))
    for (a, c), (b, d) in combinations(het_lof, 2):
        phase = phase_relation(a, b)
        if phase == "cis":
            continue
        kind = ("biallelic_lof" if c.lof_tier == d.lof_tier == "mane_select"
                else "biallelic_lof_non_mane")
        configs.append(_configuration(kind, [(a, c), (b, d)], phase))
    for a, c in het_lof:
        for b, d in altering:
            if b.zygosity == "hom_alt":
                phase = "trans"            # the second hit is on both copies
            elif b.zygosity == "het":
                phase = phase_relation(a, b)
            else:
                continue
            if phase == "cis":
                continue
            configs.append(_configuration("lof_plus_protein_altering", [(a, c), (b, d)], phase))

    order = {k: i for i, k in enumerate(GENE_VERDICTS)}
    configs.sort(key=lambda cfg: (order[cfg["class"]], cfg["variant_ids"]))
    if configs:
        verdict = configs[0]["class"]
    elif lof:
        verdict = "monoallelic_lof"
    else:
        verdict = "no_lof"
    return {
        "verdict": verdict,
        "configurations": configs,
        "n_alleles": len(items),
        "n_pass": len(passing),
        "n_lof_pass": len(lof),
        "n_protein_altering_pass": len(altering),
        "filtered_lof_variant_ids": [a.variant_id for a, c in items
                                     if a.filter != "PASS" and c.lof_tier != "not_lof"],
    }


def overall_verdict(genes: dict) -> tuple:
    """``(verdict, genes)`` across the panel; see :data:`VERDICT_NOTES`."""
    def having(v: str) -> list:
        return sorted(g for g, r in genes.items() if r["verdict"] == v)

    if having("biallelic_lof"):
        found = having("biallelic_lof")
        return ("single_candidate" if len(found) == 1 else "multiple_candidates"), found
    if having("biallelic_lof_non_mane"):
        return "non_mane_candidates_only", having("biallelic_lof_non_mane")
    if having("lof_plus_protein_altering"):
        return "lof_plus_protein_altering_only", having("lof_plus_protein_altering")
    return "no_biallelic_lof", []


def call(vcf_path, config: dict, out_dir, *, panel) -> dict:
    """Run the causal-gene half of L0 and write its two artifacts.

    Returns:
        ``{"verdict", "candidate_genes", "n_alleles", "clinvar", "gnomad"}`` -- categories
        and counts only.
    """
    policy = tx.transcript_policy(config)
    mane_release = config["annotator"]["mane_release"]
    pad = int(config.get("l0_region_pad", DEFAULT_PAD))

    regions = annotate.gene_regions(config, panel)
    alleles = read_alleles(vcf_path, regions, pad=pad)
    records = [(a.variant_id, a.contig, a.pos, a.ref, a.alt) for a in alleles]
    passes = {name: (annotate.annotate(config, records, flags) if records else {})
              for name, flags in tx.snpeff_passes(config)}

    calls, uncovered = [], []
    for allele in alleles:
        result = classify(allele, passes["primary"], passes["secondary"],
                          policy=policy, mane_release=mane_release)
        if result is None:
            uncovered.append(allele)
        else:
            calls.append((allele, result))

    genes = {g: resolve_gene([(a, c) for a, c in calls if a.gene == g]) for g in panel}
    verdict, candidates = overall_verdict(genes)

    out_dir = Path(out_dir)
    # Written before the cross-reference, which downloads a release and runs snpEff again:
    # a network failure there must not discard annotation work already finished, for the
    # same reason run.py writes the burden before this half starts.
    (out_dir / "variant_calls.json").write_text(json.dumps({
        "calls": [{"genotype": asdict(a), "annotation": asdict(c)} for a, c in calls],
        "outside_panel_transcripts": [asdict(a) for a in uncovered],
    }, indent=2), encoding="utf-8")

    if _clinvar.enabled(config):
        # A failure here is fatal on purpose. Writing the G1 artifact with the archive's
        # evidence missing invites the call to be made without it; re-running once the
        # release is reachable is the cheaper mistake. Turning the key off in config is
        # the deliberate way to run without it.
        crossref = _clinvar.cross_reference(
            config, alleles, {a.variant_id: c for a, c in calls}, regions, pad=pad)
        (out_dir / "clinvar_crossref.json").write_text(
            json.dumps(crossref, indent=2, sort_keys=True), encoding="utf-8")
        attach_clinvar(genes, crossref)
        clinvar_status = {"status": "ok", "release": crossref["release"]["file_date"],
                          "artifact": "clinvar_crossref.json"}
    else:
        # Explicit, because a missing cross-reference and an empty one read alike in an
        # artifact, and only one of them means "the archive knows nothing about these".
        clinvar_status = {"status": "disabled",
                          "reason": "config['clinvar']['enabled'] is false"}

    if _gnomad.enabled(config):
        # Fatal on failure for the same reason as the cross-reference above.
        frequencies = _gnomad.frequencies(config, alleles, regions, pad=pad)
        (out_dir / "gnomad_frequencies.json").write_text(
            json.dumps(frequencies, indent=2, sort_keys=True), encoding="utf-8")
        attach_gnomad(genes, frequencies)
        gnomad_status = {"status": "ok", "release": frequencies["release"]["version"],
                         "artifact": "gnomad_frequencies.json"}
    else:
        gnomad_status = {"status": "disabled",
                         "reason": "config['gnomad']['enabled'] is false"}

    (out_dir / "causal_gene_call.json").write_text(json.dumps({
        "verdict": verdict,
        "verdict_note": VERDICT_NOTES[verdict],
        "candidate_genes": candidates,
        "causal_gene_config": "not set by L0 -- gate G1 is a human decision made from this evidence",
        "clinvar": clinvar_status,
        "gnomad": gnomad_status,
        "per_gene": genes,
        "caveats": [*CAVEATS[:2], FREQUENCY_CAVEATS[gnomad_status["status"]], *CAVEATS[2:]],
        "provenance": {
            "annotator": config["annotator"]["tool"],
            "snpeff_version": annotate.snpeff_version(),
            "database": config["annotator"]["database"],
            "transcript_policy": policy,
            "passes": [[name, list(flags)] for name, flags in tx.snpeff_passes(config)],
            "mane_release": mane_release,
            "panel": list(panel),
            "regions": {g: list(r) for g, r in regions.items()},
            "region_pad_bp": pad,
            "lof_effects": sorted(LOF_EFFECTS),
            "filter_policy": "only FILTER=PASS alleles enter configurations",
            "phasing": "1/2 genotypes, then GATK PGT/PID; otherwise unphased",
            "decision": "mngmt/decisions.md D5",
            "citations": [
                "hanks2004 doi:10.1038/ng1449",
                "macarthur2012 doi:10.1126/science.1215040",
                "aboutayoun2018 doi:10.1002/humu.23626",
            ],
        },
        "seed": config["seed"],
    }, indent=2, sort_keys=True), encoding="utf-8")

    return {"verdict": verdict, "candidate_genes": candidates, "n_alleles": len(alleles),
            "clinvar": clinvar_status["status"], "gnomad": gnomad_status["status"]}
