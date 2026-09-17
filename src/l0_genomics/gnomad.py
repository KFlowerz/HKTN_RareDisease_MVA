"""L0 -- gnomAD population frequency: how often each panel allele occurs in the population.

Purpose
    Report, for each panel allele the subject carries, how often gnomAD observed it, so
    that gate G1 is not decided with the question "is this second hit simply common?"
    unanswered. gnomAD aggregates exome and genome sequencing across hundreds of thousands
    of individuals [chen2024] doi:10.1038/s41586-023-06045-0.

Inputs
    The alleles from :mod:`src.l0_genomics.causal`, the panel regions,
    ``config["gnomad"]`` and ``config["reference_dir"]``.

Outputs
    ``results_dir/l0_genomics/gnomad_frequencies.json`` -- per allele and per dataset
    (exomes, genomes): whether gnomAD observed it, its allele count, allele number,
    frequency, homozygote count and highest genetic-ancestry-group frequency, plus the
    allele number of nearby sites as a proxy for whether that position was sampled at all.
    A compact form is attached to each configuration in ``causal_gene_call.json``.

Method
    1. For each panel gene, the gene's whole padded span is read from gnomAD's public
       per-chromosome sites VCF with ``bcftools query``, and cached as a table outside the
       repository (:mod:`src.refcache`). Only the core count and frequency fields are read.
    2. Alleles match on ``(contig, pos, ref, alt)`` reduced to the minimal representation
       [tan2015] doi:10.1093/bioinformatics/btv112, the same rule as the ClinVar
       cross-reference, with the ``chr`` prefix normalised.
    3. Absence is only reported alongside the allele number of sites within
       :data:`NEARBY_BP` of the allele. A variant gnomAD did not observe at a position it
       barely sampled is not rare; it is unmeasured.

Guardrail
    **Whole gene spans, never the subject's coordinates** -- decision D8. The whole
    release is 24 GB for one chromosome, so unlike ClinVar it is not downloaded; instead
    htslib reads the byte ranges covering each *panel gene's* span. Those spans come from
    the panel in ``config`` and are identical for every proband run through it, so the
    request says which genes the project studies -- already public in this repository --
    and nothing about which variants this child carries. :func:`extract` takes a gene
    span and nothing else, and ``tests/test_l0_gnomad.py`` asserts the queried regions do
    not depend on the alleles.

    **Only CC0 fields are read.** gnomAD's primary data are dedicated to the public domain
    (CC0); some annotations shipped in the same files are not -- SpliceAI's are CC BY-NC
    4.0 [gnomad_policies]. :data:`FIELDS` names every field read; the rest are never
    parsed.

    **Frequency is reported, never interpreted.** Rarity is necessary but not sufficient
    for pathogenicity, and nothing here applies ACMG/AMP criterion PM2 or any other
    [richards2015] doi:10.1038/gim.2015.30.

    **The output identifies the child**: it pairs their alleles with frequencies. It is
    written under ``results_dir`` only and never logged.
"""

from __future__ import annotations

import json
import statistics
import subprocess
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from .. import refcache
from .clinvar import minimal

RELEASE = "4.1.1"
DATASETS = ("exomes", "genomes")
URL_TEMPLATE = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/{release}/"
                "vcf/{dataset}/gnomad.{dataset}.v{release}.sites.{contig}.vcf.bgz")
#: gnomAD policies, read 2026-09-17 [gnomad_policies].
LICENCE = "CC0 1.0 (core fields only; SpliceAI and other licensed annotations are not read)"

#: Every INFO field read. AF_grpmax is the highest frequency across genetic ancestry
#: groups; nhomalt the number of individuals homozygous for the allele.
FIELDS = ("AC", "AN", "AF", "nhomalt", "AF_grpmax")
#: Window for the nearby-sites sampling proxy.
NEARBY_BP = 50

CAVEATS = (
    "Rarity is necessary but not sufficient for pathogenicity: most very rare variants "
    "are benign. No ACMG/AMP criterion is applied here.",
    "Presence in gnomAD does not make an allele benign for a recessive disease: healthy "
    "heterozygous carriers of a pathogenic allele are expected in any large population "
    "sample. The homozygote count is the more informative number.",
    "Absence is informative only where the position was sampled. The allele number of "
    "nearby sites is reported as a proxy; gnomAD's per-base coverage files are not read, "
    "and a region with no nearby variant sites gives no proxy at all.",
    "Matching is on the minimal representation of (contig, position, ref, alt). Two "
    "sources that left-align the same indel to different positions would not match, and "
    "the variant would be reported as not observed.",
    "Frequencies differ between genetic ancestry groups. The subject's ancestry is not "
    "assessed and must not be inferred from these numbers; AF_grpmax reports the highest "
    "group frequency instead.",
    "Sites gnomAD itself did not pass are reported with their filter, not dropped, and "
    "are not counted as observed: site_in_release says the release holds the site, "
    "observed says at least one carrier was counted there.",
)


@dataclass(frozen=True)
class Site:
    """One allele at one gnomAD site."""

    contig: str
    pos: int
    ref: str
    alt: str
    filter: str
    ac: int | None
    an: int | None
    af: float | None
    nhomalt: int | None
    af_grpmax: float | None

    @property
    def key(self) -> tuple:
        return (bare_contig(self.contig), *minimal(self.pos, self.ref, self.alt))


def bare_contig(contig: str) -> str:
    """``chr15`` and ``15`` compare equal."""
    return contig[3:] if contig.startswith("chr") else contig


def remote_contig(contig: str) -> str:
    """gnomAD's GRCh38 files name contigs with the ``chr`` prefix."""
    return contig if contig.startswith("chr") else f"chr{contig}"


def _number(text: str, kind):
    """A field value, or None for bcftools' missing marker or a missing field."""
    if text in (".", ""):
        return None
    return kind(text)


def parse_sites(lines) -> list:
    """Parse ``bcftools query`` rows: CHROM POS REF ALT FILTER then :data:`FIELDS`.

    A multi-allelic row, whose per-allele fields are comma-separated, becomes one
    :class:`Site` per ALT; AN is shared. gnomAD's sites files are split already, so this
    is defensive -- which is why it is tested rather than assumed.
    """
    sites = []
    for line in lines:
        if not line.strip():
            continue
        contig, pos, ref, alts, filt, ac, an, af, nhomalt, grpmax = line.rstrip("\n").split("\t")
        alts = alts.split(",")

        def per_allele(text, i, kind):
            parts = text.split(",")
            return _number(parts[i] if len(parts) == len(alts) else ".", kind)

        for i, alt in enumerate(alts):
            sites.append(Site(
                contig=contig, pos=int(pos), ref=ref, alt=alt, filter=filt,
                ac=per_allele(ac, i, int), an=_number(an, int), af=per_allele(af, i, float),
                nhomalt=per_allele(nhomalt, i, int), af_grpmax=per_allele(grpmax, i, float),
            ))
    return sites


def gnomad_dir(config: dict) -> Path:
    """The gnomAD extracts' cache, inside ``reference_dir``."""
    path = refcache.reference_dir(config) / "gnomad"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _remote_headers(url: str) -> dict:
    """What the server says about the file, so an extract traces to a release object."""
    request = urllib.request.Request(url, method="HEAD", headers=refcache.USER_AGENT)
    with urllib.request.urlopen(request, timeout=60) as response:
        headers = response.headers
        return {"last_modified": headers.get("Last-Modified"),
                "etag": headers.get("ETag"),
                "bytes": int(headers.get("Content-Length") or 0) or None}


def extract(config: dict, dataset: str, gene: str, span: tuple) -> tuple:
    """Read one panel gene's whole span from gnomAD; return ``(path, provenance)``.

    ``span`` is ``(contig, start, end)`` of the *gene*, already padded. This function has
    no access to any allele, by design (decision D8).

    Raises:
        RuntimeError: If ``bcftools`` fails -- a network error must not read as "no
            variant in this gene".
    """
    release = (config.get("gnomad") or {}).get("release", RELEASE)
    contig, start, end = span
    url = URL_TEMPLATE.format(release=release, dataset=dataset, contig=remote_contig(contig))
    cache = gnomad_dir(config)
    stem = f"{gene}.{dataset}.v{release}.{remote_contig(contig)}_{start}_{end}"
    table, sidecar = cache / f"{stem}.tsv", cache / f"{stem}.json"

    if not (table.exists() and sidecar.exists()):
        remote = _remote_headers(url)
        region = f"{remote_contig(contig)}:{start}-{end}"
        fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t" + "\\t".join(
            f"%INFO/{f}" for f in FIELDS) + "\\n"
        part = table.with_name(table.name + ".part")
        # htslib downloads the remote index into the working directory; keep it in the
        # cache rather than wherever the pipeline was launched from.
        try:
            with open(part, "w", encoding="utf-8") as handle:
                result = subprocess.run(["bcftools", "query", "-r", region, "-f", fmt, url],
                                        stdout=handle, stderr=subprocess.PIPE, text=True,
                                        cwd=cache, timeout=1800, check=False)
        except BaseException:
            # A timeout or an interrupt must not leave a partial table to be taken for whole.
            part.unlink(missing_ok=True)
            raise
        finally:
            for index in cache.glob("*.tbi"):
                index.unlink()
        if result.returncode != 0:
            part.unlink(missing_ok=True)
            raise RuntimeError(f"bcftools could not read gnomAD {dataset} for {gene}: "
                               f"{result.stderr.strip()[-500:]}")
        part.rename(table)
        sidecar.write_text(json.dumps({
            "url": url, "region": region, "remote": remote,
            "retrieved": date.today().isoformat()}, indent=2), encoding="utf-8")

    provenance = json.loads(sidecar.read_text(encoding="utf-8"))
    provenance.update({"file": table.name, "sha256": refcache.sha256(table)})
    return table, provenance


def load_sites(path: Path, gene: str, dataset: str) -> list:
    """:func:`parse_sites` over one extract.

    Raises:
        ValueError: If the extract holds no site. Every panel gene has variant sites in
            gnomAD, so an empty table means a wrong contig name or a failed read.
    """
    with open(path, encoding="utf-8") as handle:
        sites = parse_sites(handle)
    if not sites:
        raise ValueError(f"gnomAD {dataset} extract for {gene} ({path.name}) holds no site "
                         "-- a contig-naming mismatch or an incomplete read?")
    return sites


def lookup(allele, sites: list) -> dict:
    """One allele against one dataset's sites for its gene."""
    key = (bare_contig(allele.contig), *minimal(allele.pos, allele.ref, allele.alt))
    exact = [s for s in sites if s.key == key]
    nearby = [s.an for s in sites
              if bare_contig(s.contig) == key[0] and abs(s.pos - allele.pos) <= NEARBY_BP
              and s.an is not None]
    an_max = max((s.an for s in sites if s.an is not None), default=None)
    median = statistics.median(nearby) if nearby else None
    # "Observed" means gnomAD counted at least one carrier. A site it kept but flagged, one
    # whose carriers all failed its filters (AC=0), and one whose release reports no AC at
    # all are each present without anyone having been seen to carry the allele; reporting
    # that as observed would read, at gate G1, as "the population carries this".
    carriers = exact[0].ac if exact else None
    entry = {
        "observed": bool(exact) and bool(carriers),
        "site_in_release": bool(exact),
        "n_sites_nearby": len(nearby),
        "an_nearby_median": median,
        "an_max_in_gene": an_max,
        "an_nearby_fraction_of_max": (round(median / an_max, 4)
                                      if median is not None and an_max else None),
    }
    if exact:
        site = exact[0]
        entry.update({"ac": site.ac, "an": site.an, "af": site.af, "nhomalt": site.nhomalt,
                      "af_grpmax": site.af_grpmax, "filter": site.filter})
    return entry


def frequencies(config: dict, alleles, regions: dict, *, pad: int) -> dict:
    """Look up every panel allele in each configured gnomAD dataset."""
    settings = config.get("gnomad") or {}
    datasets = tuple(settings.get("datasets", DATASETS))
    # Spans come from the panel regions alone -- never from the alleles (decision D8).
    spans = {gene: (contig, max(0, start - pad), end + pad)
             for gene, (contig, start, end) in sorted(regions.items())}

    sites: dict = {}
    extracts: dict = {}
    for dataset in datasets:
        for gene, span in spans.items():
            path, provenance = extract(config, dataset, gene, span)
            sites[(dataset, gene)] = load_sites(path, gene, dataset)
            extracts.setdefault(dataset, {})[gene] = {
                **provenance, "n_sites": len(sites[(dataset, gene)])}

    entries = []
    for allele in alleles:
        entry = {"variant_id": allele.variant_id, "gene": allele.gene}
        for dataset in datasets:
            entry[dataset] = lookup(allele, sites.get((dataset, allele.gene), []))
        entries.append(entry)

    return {
        "status": "ok",
        "release": {"source": "gnomAD", "version": settings.get("release", RELEASE),
                    "licence": LICENCE, "datasets": list(datasets), "extracts": extracts},
        "alleles": entries,
        "provenance": {
            "matching": "minimal representation of (contig, pos, ref, alt) [tan2015]",
            "fields_read": list(FIELDS),
            "query": "whole padded panel-gene spans via htslib, never subject coordinates "
                     "(decision D8)",
            "region_pad_bp": pad,
            "nearby_bp": NEARBY_BP,
            "decision": "mngmt/decisions.md D8",
            "citations": [
                "chen2024 doi:10.1038/s41586-023-06045-0",
                "tan2015 doi:10.1093/bioinformatics/btv112",
                "richards2015 doi:10.1038/gim.2015.30",
                "gnomad_policies",
            ],
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }


def summarize(entry: dict, datasets=DATASETS) -> dict:
    """The compact form attached to a configuration in ``causal_gene_call.json``."""
    out = {"variant_id": entry["variant_id"]}
    for dataset in datasets:
        d = entry.get(dataset)
        if d is None:
            continue
        # filter and af travel with the counts: a flagged site and a clean one are not the
        # same evidence, and the compact form is what a G1 reviewer reads first.
        out[dataset] = {k: d.get(k) for k in
                        ("observed", "site_in_release", "ac", "an", "af", "nhomalt",
                         "filter", "an_nearby_fraction_of_max")}
    return out


def enabled(config: dict) -> bool:
    """Whether the lookup runs. Default on; ``gnomad.enabled: false`` turns it off."""
    return bool((config.get("gnomad") or {}).get("enabled", True))
