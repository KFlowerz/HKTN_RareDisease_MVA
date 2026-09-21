"""L0 -- ClinVar cross-reference: what the public archive already says about these alleles.

Purpose
    Attach the existing clinical interpretation, if any, to each panel allele the subject
    carries, so that gate G1 is decided against what has already been submitted about
    these variants rather than against consequence prediction alone. ClinVar aggregates
    submitted interpretations with the evidence and review status behind them
    [landrum2018] doi:10.1093/nar/gkx1153.

Inputs
    The alleles and :class:`~src.l0_genomics.transcripts.VariantCall` records from
    :mod:`src.l0_genomics.causal`, the panel regions, ``config["clinvar"]`` and
    ``config["reference_dir"]``.

Outputs
    ``results_dir/l0_genomics/clinvar_crossref.json`` -- per allele: the exact ClinVar
    record if one exists, other records at the same position, and pathogenic records
    affecting the same protein residue; per gene: how many records the archive holds, as
    the denominator; plus the release's provenance and the caveats below. A compact form
    of the same evidence is attached to each configuration in ``causal_gene_call.json``.

Method
    1. The whole GRCh38 release is downloaded once and cached (:mod:`src.refcache`).
    2. Records inside the padded panel regions are parsed; the rest of the archive is
       discarded.
    3. Alleles match on ``(contig, pos, ref, alt)`` reduced to the minimal representation
       [tan2015] doi:10.1093/bioinformatics/btv112, so the two files' padding conventions
       cannot hide a match. Records at the same position with a different allele are
       reported separately -- they are context, not a match.
    4. Pathogenic and likely-pathogenic records in the panel regions are annotated by the
       same MANE Select pass decision D5 prescribes, and those changing the same protein
       residue as one of the subject's alleles are reported. Same-residue and same-change
       pathogenic records are the observations behind ACMG/AMP criteria PM5 and PS1
       [richards2015] doi:10.1038/gim.2015.30 -- reported here as what the archive holds,
       never applied as a criterion.

Guardrail
    **The whole release is downloaded; no variant is ever queried.** A per-variant lookup
    against ClinVar's API would put this child's coordinates on an NCBI server, which
    ``COMPLIANCE.md`` prohibits -- the same reason annotation is local
    (:mod:`src.l0_genomics.annotate`). Matching happens entirely on this machine.

    **This is not a variant classification.** ClinVar aggregates submissions of varying
    quality, and the review status travels with every record reported here for exactly
    that reason. Nothing in this module applies ACMG/AMP criteria, assigns pathogenicity,
    or decides gate G1; it reports what the archive holds so a person can weigh it.

    **Absence is not evidence of benignity.** A variant absent from ClinVar has not been
    submitted -- in a disease with roughly 50 patients worldwide, that is the expected
    state for a private allele, not a finding.

    **The output identifies the child**: it pairs their alleles with named conditions. It
    is written under ``results_dir`` only and never logged; the returned summary carries
    counts and categories.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from .. import refcache
from . import annotate
from . import transcripts as tx

#: A **dated** weekly release, never the rolling ``clinvar.vcf.gz``. That file is always
#: reachable and never reproducible: NCBI replaces its contents every week, so two runs a
#: fortnight apart cross-reference the subject's alleles against different archives and
#: nothing in the output would say so. The release actually used is pinned in
#: ``config/pipeline.yaml`` and recorded in ``docs/references.md``.
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260913.vcf.gz"
#: Where NCBI moves a dated release once the next one lands. Tried after the primary, so
#: the pin keeps working when the release rotates out of the top-level directory.
CLINVAR_ARCHIVE_URL = ("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
                       "archive_2.0/2026/clinvar_20260913.vcf.gz")
#: ClinVar is a US-Government work in the public domain; NCBI asks for citation, not
#: permission. Clean for a redistributed CC-BY-4.0 output (``../l4_validate/sources.md``).
LICENCE = "US-Gov public domain"

#: Review status to star rating, as ClinVar defines it [clinvar_reviewstatus]. An unknown
#: status maps to ``None`` -- a guessed 0 would silently understate a record's support.
REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_multiple_submitters": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_interpretation_for_the_single_variant": 0,
    "no_classification_for_the_single_variant": 0,
    "no_classifications_from_unflagged_records": 0,
}

#: Significance classes that count as "the archive calls this pathogenic".
PATHOGENIC_CLASSES = frozenset({"pathogenic", "likely_pathogenic"})

#: Protein position in an HGVS.p string, e.g. ``p.Arg550Ter`` -> 550.
_HGVS_P_RE = re.compile(r"^p\.[A-Za-z]{3}(\d+)")

CAVEATS = (
    "Absence from ClinVar is not evidence of benignity. It means no submitter has "
    "deposited the variant -- the expected state for a private allele in a disease with "
    "roughly 50 known patients.",
    "ClinVar aggregates submissions of varying quality. A one-star record is one "
    "submitter's opinion; weigh every match by its review status, which is reported with "
    "it.",
    "Matching is on the minimal representation of (contig, position, ref, alt). Two "
    "sources that left-align the same indel to different positions would not match, and "
    "the variant would be reported as absent.",
    "Same-residue pathogenic records are reported as an observation, not as an applied "
    "ACMG/AMP criterion. PM5 and PS1 carry conditions -- the same transcript, a "
    "comparable amino acid change, and splice effects excluded -- that are not evaluated "
    "here.",
    "The allele frequencies carried in ClinVar records come from ESP, ExAC and 1000 "
    "Genomes. They are legacy fields, absent from most records, and are no substitute "
    "for a gnomAD lookup.",
)


@dataclass(frozen=True)
class ClinVarRecord:
    """One ClinVar VCF record inside a panel region."""

    variation_id: str
    contig: str
    pos: int
    ref: str
    alt: str
    clnsig: str                  # CLNSIG, verbatim
    significance: str            # its class -- see :func:`significance_class`
    review: str                  # CLNREVSTAT, verbatim
    stars: int | None            # REVIEW_STARS[review], or None if unrecognised
    conditions: str = ""         # CLNDN
    variant_type: str = ""       # CLNVC
    molecular_consequence: str = ""   # MC, Sequence Ontology terms
    origin: str = ""             # ORIGIN
    gene_info: str = ""          # GENEINFO
    frequencies: dict = field(default_factory=dict)   # AF_ESP / AF_EXAC / AF_TGP

    @property
    def key(self) -> tuple:
        """The allele's minimal representation, for matching."""
        return (self.contig, *minimal(self.pos, self.ref, self.alt))


def minimal(pos: int, ref: str, alt: str) -> tuple:
    """Reduce ``(pos, ref, alt)`` to its minimal representation.

    Trailing then leading bases shared by REF and ALT are removed, which is what makes two
    VCFs' padding conventions comparable [tan2015] doi:10.1093/bioinformatics/btv112. It
    does *not* re-align a repeat: two sources that left-aligned the same indel differently
    still will not match, which is why that limit is a stated caveat.
    """
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def significance_class(clnsig: str) -> str:
    """Bucket a ``CLNSIG`` value.

    ClinVar packs several statements into one field: alternatives separated by ``/``,
    additional assertions (drug response, risk factor) after ``|``. The leading term is
    the pathogenicity statement; the rest is recorded verbatim on the record instead.
    """
    parts = {p.replace("_low_penetrance", "").strip("_")
             for p in clnsig.split("|")[0].split("/")}
    if parts & {"Conflicting_classifications_of_pathogenicity",
                "Conflicting_interpretations_of_pathogenicity"}:
        return "conflicting"
    if "Pathogenic" in parts:
        return "pathogenic"
    if "Likely_pathogenic" in parts:
        return "likely_pathogenic"
    if "Benign" in parts:
        return "benign"
    if "Likely_benign" in parts:
        return "likely_benign"
    if parts & {"Uncertain_significance", "Uncertain_risk_allele"}:
        return "uncertain"
    return "other"


def protein_position(hgvs_p: str):
    """The residue number in an HGVS.p string, or None if it names none."""
    match = _HGVS_P_RE.match(hgvs_p or "")
    return int(match.group(1)) if match else None


def _info(field_text: str) -> dict:
    return dict(kv.split("=", 1) for kv in field_text.split(";") if "=" in kv)


def parse_records(lines, regions: dict, *, pad: int) -> tuple:
    """Parse VCF text, keeping records inside the padded panel regions.

    Returns ``(file_date, records)``. Every ALT of a multi-allelic record becomes its own
    :class:`ClinVarRecord`, since matching is per allele. Takes lines rather than a path so
    it can be exercised without a VCF file existing anywhere.
    """
    # One interval per gene, not a merged span per contig: two panel genes on the same
    # contig (which config's l0_panel_extension allows) would otherwise pull in every
    # record between them -- tens of thousands, annotated and matched for nothing.
    spans: dict = {}
    for contig, start, end in regions.values():
        spans.setdefault(contig, []).append((max(0, start - pad), end + pad))

    file_date, records = "", []
    for line in lines:
        if line.startswith("#"):
            if line.startswith("##fileDate="):
                file_date = line.strip().split("=", 1)[1]
            continue
        fields = line.rstrip("\n").split("\t")
        intervals = spans.get(fields[0])
        if intervals is None:
            continue
        pos = int(fields[1])
        if not any(lo <= pos <= hi for lo, hi in intervals) or fields[4] in (".", ""):
            continue
        info = _info(fields[7])
        clnsig = info.get("CLNSIG", "")
        review = info.get("CLNREVSTAT", "")
        frequencies = {k.lower(): float(info[k]) for k in ("AF_ESP", "AF_EXAC", "AF_TGP")
                       if k in info}
        for alt in fields[4].split(","):
            records.append(ClinVarRecord(
                variation_id=fields[2], contig=fields[0], pos=pos, ref=fields[3], alt=alt,
                clnsig=clnsig, significance=significance_class(clnsig), review=review,
                stars=REVIEW_STARS.get(review), conditions=info.get("CLNDN", ""),
                variant_type=info.get("CLNVC", ""),
                molecular_consequence=",".join(
                    part.split("|")[1] for part in info.get("MC", "").split(",")
                    if "|" in part),
                origin=info.get("ORIGIN", ""), gene_info=info.get("GENEINFO", ""),
                frequencies=frequencies,
            ))
    return file_date, records


def load_records(path: Path, regions: dict, *, pad: int) -> tuple:
    """:func:`parse_records` over the cached release.

    Raises:
        ValueError: If no record falls in any panel region. Six well-studied genes always
            have ClinVar records, so an empty result means the contig naming does not
            match -- which would otherwise read as "none of these variants are known".
    """
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        file_date, records = parse_records(handle, regions, pad=pad)
    if not records:
        raise ValueError(
            f"no ClinVar record falls in any panel region of {path.name} -- a contig-naming "
            "mismatch? The panel genes all have ClinVar records."
        )
    return file_date, records


def gene_of(record: ClinVarRecord, regions: dict, *, pad: int) -> str:
    """Which panel gene's region a record sits in, or ``""``.

    Genes are tried in symbol order, so a record inside two overlapping panel regions is
    assigned deterministically rather than by dict order.
    """
    for gene, (contig, start, end) in sorted(regions.items()):
        if record.contig == contig and start - pad <= record.pos <= end + pad:
            return gene
    return ""


def annotate_pathogenic(config: dict, records) -> dict:
    """``{record.key: ((transcript, hgvs_p), ...)}`` on the call-making transcript(s).

    Keyed by allele, not by variation id: ClinVar gives every ALT of a multi-allelic
    record the same variation id, so keying by id would attach one allele's amino-acid
    change to another's. The panel's current records happen to be biallelic throughout,
    which is exactly why this would go unnoticed.

    Only pathogenic and likely-pathogenic records are annotated: a same-residue benign or
    uncertain record says nothing a person would act on, and the pass costs a snpEff run.
    The pass is the primary (MANE Select) one from decision D5, so the residue numbers
    compare against the subject's calls, which are made on the same transcript.
    """
    wanted = [r for r in records if r.significance in PATHOGENIC_CLASSES]
    if not wanted:
        return {}
    flags = dict(tx.snpeff_passes(config))["primary"]
    keyed = {f"cv{i}": r for i, r in enumerate(wanted)}
    annotations = annotate.annotate(
        config, [(k, r.contig, r.pos, r.ref, r.alt) for k, r in keyed.items()], flags)

    # Every transcript-level hit is kept, not the first: a record can be annotated on more
    # than one gene's MANE Select transcript where spans overlap, and keeping only the
    # first would silently drop the one the comparison needs.
    out = {}
    for name, record in keyed.items():
        anns, _ = annotations.get(name, ((), frozenset()))
        hits = tuple((a.transcript, a.hgvs_p) for a in anns
                     if a.feature_type == "transcript" and protein_position(a.hgvs_p) is not None)
        if hits:
            out[record.key] = hits
    return out


def _record_summary(record: ClinVarRecord) -> dict:
    return {"variation_id": record.variation_id, "clnsig": record.clnsig,
            "significance": record.significance, "review": record.review,
            "stars": record.stars, "conditions": record.conditions,
            "variant_type": record.variant_type,
            "molecular_consequence": record.molecular_consequence,
            "origin": record.origin, "frequencies": record.frequencies}


def match_allele(allele, call, records, protein: dict) -> dict:
    """Cross-reference one allele against the panel's ClinVar records.

    ``records`` are the records of that allele's gene; ``protein`` is
    :func:`annotate_pathogenic`'s map. Returns the allele's entry in the artifact, in
    which the three findings are disjoint: an exact match, other records at the same
    position, and other pathogenic records changing the same residue.
    """
    key = (allele.contig, *minimal(allele.pos, allele.ref, allele.alt))
    exact = [r for r in records if r.key == key]
    same_position = [r for r in records if r.pos == allele.pos and r.key != key]

    residue = protein_position(getattr(call, "hgvs_p", "")) if call is not None else None
    same_residue = []
    if residue is not None:
        for record in records:
            # The exact match is reported above; repeating it here would read as a second,
            # independent record supporting the same residue.
            if record.key == key:
                continue
            for transcript, hgvs_p in protein.get(record.key, ()):
                if transcript != call.transcript or protein_position(hgvs_p) != residue:
                    continue
                same_residue.append({**_record_summary(record), "transcript": transcript,
                                     "hgvs_p": hgvs_p,
                                     "same_change": hgvs_p == call.hgvs_p})
                break

    return {
        "variant_id": allele.variant_id,
        "gene": allele.gene,
        "match": "exact" if exact else ("same_position" if same_position else "none"),
        "records": [_record_summary(r) for r in exact],
        "same_position_records": [_record_summary(r) for r in same_position],
        "same_protein_position_pathogenic": sorted(
            same_residue, key=lambda r: (not r["same_change"], r["variation_id"])),
    }


def cross_reference(config: dict, alleles, calls: dict, regions: dict, *, pad: int) -> dict:
    """Cross-reference every panel allele. ``calls`` maps variant id to its VariantCall."""
    settings = config.get("clinvar") or {}
    url = settings.get("url", CLINVAR_URL)
    # Same release, other location. NCBI rotates a dated file into archive_2.0/ when the
    # next weekly lands, so the pin has to name both or it expires on its own.
    fallbacks = list(settings.get("fallback_urls") or
                     ([CLINVAR_ARCHIVE_URL] if url == CLINVAR_URL else []))
    release = refcache.reference_dir(config) / url.rsplit("/", 1)[-1]
    provenance = refcache.fetch_first([url, *fallbacks], release)
    file_date, records = load_records(release, regions, pad=pad)

    by_gene: dict = {gene: [] for gene in regions}
    for record in records:
        gene = gene_of(record, regions, pad=pad)
        if gene:
            by_gene[gene].append(record)
    protein = annotate_pathogenic(config, records)

    matches = [match_allele(a, calls.get(a.variant_id), by_gene.get(a.gene, ()), protein)
               for a in alleles]

    # Format drift is the one failure a test suite cannot catch: invented rows model
    # ClinVar's fields, so a renamed INFO key would pass every test and silently yield
    # records with no significance and no stars. Counting what did not map makes that
    # visible in the artifact instead of invisible.
    unmapped: dict = {}
    for record in records:
        if record.stars is None:
            unmapped[record.review] = unmapped.get(record.review, 0) + 1
    return {
        "status": "ok",
        "release": {"source": "ClinVar", "licence": LICENCE, "file_date": file_date,
                    **provenance},
        "gene_context": {
            gene: {"n_records": len(rs),
                   "n_pathogenic": sum(r.significance in PATHOGENIC_CLASSES for r in rs),
                   "n_annotated_on_call_transcript": sum(r.key in protein for r in rs)}
            for gene, rs in by_gene.items()},
        "alleles": matches,
        "format_check": {
            "n_records_parsed": len(records),
            "n_significance_other": sum(r.significance == "other" for r in records),
            "unmapped_review_statuses": dict(sorted(unmapped.items())),
            "interpretation": (
                "Counts of what this release's vocabulary did not map onto. A handful is "
                "normal -- some records carry no review status at all. A count near "
                "n_records_parsed means ClinVar's field names have changed and the "
                "mapping in src/l0_genomics/clinvar.py needs revisiting; the matches "
                "above would still be correct, but their significance would not be."
            ),
        },
        "provenance": {
            "matching": "minimal representation of (contig, pos, ref, alt) [tan2015]",
            "region_pad_bp": pad,
            "annotation_pass": "primary (MANE Select), decision D5",
            "snpeff_version": annotate.snpeff_version(),
            "database": config["annotator"]["database"],
            "citations": [
                "landrum2018 doi:10.1093/nar/gkx1153",
                "tan2015 doi:10.1093/bioinformatics/btv112",
                "richards2015 doi:10.1038/gim.2015.30",
                "clinvar_reviewstatus",
            ],
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }


def summarize(match: dict) -> dict:
    """The compact form attached to a configuration in ``causal_gene_call.json``."""
    top = match["records"][0] if match["records"] else None
    return {
        "variant_id": match["variant_id"],
        "match": match["match"],
        "significance": top["significance"] if top else None,
        "stars": top["stars"] if top else None,
        "n_same_position_records": len(match["same_position_records"]),
        "n_same_protein_position_pathogenic": len(match["same_protein_position_pathogenic"]),
    }


def enabled(config: dict) -> bool:
    """Whether the cross-reference runs. Default on; ``clinvar.enabled: false`` turns it off."""
    return bool((config.get("clinvar") or {}).get("enabled", True))
