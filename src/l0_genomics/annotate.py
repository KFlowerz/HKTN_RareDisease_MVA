"""L0 -- local variant annotation with snpEff.

Purpose
    Run snpEff over variant records and parse its ``ANN`` field, with every
    privacy-relevant flag enforced in one place so no caller can forget one.

Inputs
    ``config["annotator"]`` (``database``, ``java_heap``); variant records as
    ``(variant_id, contig, pos, ref, alt)`` tuples; gene symbols.

Outputs
    :class:`Annotation` records per variant id, snpEff's ``NMD`` tag genes, and gene
    regions from the local database.

Guardrail
    **Annotation is local.** Ensembl's VEP REST API would put this child's variants on a
    third-party server, which ``COMPLIANCE.md`` prohibits and no convenience justifies.
    Gene *symbols* are public identifiers; variants never leave the machine.

    Records reach snpEff on stdin, so no ``.vcf`` is written and no dataset filename
    appears in a process listing, and the ``ID`` column carries an internal id, never the
    dataset's own. Every invocation carries ``-noLog`` (snpEff otherwise reports usage
    statistics to its server) and ``-nodownload`` (otherwise a missing database is
    fetched mid-run); ``ann`` also carries ``-noStats`` (otherwise summary files describing
    the sample land in the working directory). snpEff runs from a throwaway directory, and
    its stderr is discarded because its warnings can echo the records it was given.
"""

from __future__ import annotations

import subprocess
import tempfile
import threading
from dataclasses import dataclass

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
#: Carried by every snpEff invocation; ``ann`` additionally gets ``-noStats``.
PRIVACY_FLAGS = ("-noLog", "-nodownload")


@dataclass(frozen=True)
class Annotation:
    """One snpEff ``ANN`` entry: one allele's consequence on one feature."""

    allele: str
    effects: tuple        # e.g. ("stop_gained", "splice_region_variant")
    impact: str           # HIGH | MODERATE | LOW | MODIFIER
    gene: str
    feature_type: str     # "transcript" for transcript-level consequences
    transcript: str
    biotype: str
    rank: str             # exon or intron rank, e.g. "3/10"
    hgvs_c: str
    hgvs_p: str


def snpeff_argv(config: dict, command: str, *args: str) -> list:
    """Build a snpEff command line carrying every mandatory flag."""
    extra = ("-noStats",) if command in ("ann", "eff") else ()
    return ["snpEff", f"-Xmx{config['annotator']['java_heap']}", command,
            *PRIVACY_FLAGS, *extra, *args]


def run_snpeff(config: dict, command: str, *args: str, stdin=None):
    """Yield snpEff's stdout lines. ``stdin`` is an iterable of text, fed from a thread.

    Raises:
        RuntimeError: If snpEff exits non-zero. Its stderr is not included, because it
            can echo variant records.
    """
    with tempfile.TemporaryDirectory(prefix="snpeff-") as cwd:
        proc = subprocess.Popen(
            snpeff_argv(config, command, *args),
            stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=cwd,
        )
        if stdin is not None:
            def feed() -> None:
                try:
                    for chunk in stdin:
                        proc.stdin.write(chunk)
                except BrokenPipeError:
                    pass
                finally:
                    try:
                        proc.stdin.close()
                    except BrokenPipeError:
                        pass

            threading.Thread(target=feed, daemon=True).start()
        yield from proc.stdout
        if proc.wait():
            raise RuntimeError(
                f"snpEff {command} exited with status {proc.returncode} "
                "(stderr suppressed: it can echo variant records)"
            )


def snpeff_version() -> str:
    """snpEff's version string, for provenance."""
    return subprocess.run(["snpEff", "-noLog", "-version"], capture_output=True,
                          text=True).stdout.strip()


def gene_regions(config: dict, genes) -> dict:
    """``{gene: (contig, start, end)}`` for gene symbols, from the local database.

    Raises:
        ValueError: If a symbol is absent -- a panel gene silently missing would read as
            "no variants in that gene".
    """
    db = config["annotator"]["database"]
    regions = {}
    for line in run_snpeff(config, "genes2bed", db, *genes):
        if line.startswith("#") or not line.strip():
            continue
        contig, start, end, name = line.rstrip("\n").split("\t")[:4]
        regions[name.split(";")[0]] = (contig, int(start), int(end))
    missing = [g for g in genes if g not in regions]
    if missing:
        raise ValueError(f"gene symbol(s) not found in {db}: {missing}")
    return regions


def _info_value(info: str, key: str) -> str:
    return next((kv[len(key) + 1:] for kv in info.split(";") if kv.startswith(key + "=")), "")


def parse_ann(info: str) -> tuple:
    """Parse the ``ANN`` field of a VCF INFO string into :class:`Annotation` records."""
    out = []
    for entry in _info_value(info, "ANN").split(","):
        f = entry.split("|")
        if len(f) < 11:
            continue
        out.append(Annotation(
            allele=f[0], effects=tuple(f[1].split("&")), impact=f[2], gene=f[3],
            feature_type=f[5], transcript=f[6], biotype=f[7], rank=f[8],
            hgvs_c=f[9], hgvs_p=f[10],
        ))
    return tuple(out)


def parse_tag_genes(info: str, tag: str) -> frozenset:
    """Genes named in a snpEff ``LOF``/``NMD`` tag, e.g. ``NMD=(GENE|ENSG...|1|1.00)``."""
    value = _info_value(info, tag)
    return frozenset(
        part.strip("()").split("|")[0] for part in value.split(",") if part.strip("()")
    )


def annotate(config: dict, records, flags) -> dict:
    """Annotate records in one pass: ``{variant_id: (annotations, nmd_genes)}``.

    Args:
        config: Pipeline configuration (``annotator``).
        records: ``(variant_id, contig, pos, ref, alt)`` tuples, one ALT each.
        flags: The pass's transcript flags, from
            :func:`~src.l0_genomics.transcripts.snpeff_passes`.
    """
    lines = [VCF_HEADER] + [
        f"{contig}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n"
        for vid, contig, pos, ref, alt in records
    ]
    out = {}
    for line in run_snpeff(config, "ann", "-lof", *flags, config["annotator"]["database"], "-",
                           stdin=lines):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        out[f[2]] = (parse_ann(f[7]), parse_tag_genes(f[7], "NMD"))
    return out
