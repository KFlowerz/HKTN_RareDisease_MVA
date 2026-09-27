"""Check the Track 1 submission CSV on disk, the way the organizers' parser will read it.

**Purpose.** ``build_track1_submission.py`` validates rows *before* it writes them. This
reads the written file back and validates it as a stranger would, which is the only way to
catch what happens afterwards: a hand edit, a rename, a spreadsheet round-trip that adds a
byte-order mark or reformats a number. A submission slot is spent whether or not the file
was correct when it was generated.

**Inputs.** A path to the CSV. Defaults to the packaged
``results/submissions/track1/track1_submission.csv``.

**Outputs.** A report on stdout and an exit status: ``0`` if the file is uploadable, ``1``
if an error would cost the slot. Warnings never fail the run; they are things only the
maintainer can judge.

**Guardrail.** This file is the subject's variant coordinates (decision D21), and a
terminal transcript is custody location ``C7``. So **nothing here prints a field value**
from a coordinate column. Every finding names a row number and a column, never what is in
it. ``proband_id`` is the one value printed, because confirming it is the point and it is
an assigned study identifier rather than a coordinate.

Usage::

    python scripts/check_track1_submission.py [path/to/submission.csv]
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT = REPO_ROOT / "results" / "submissions" / "track1" / "track1_submission.csv"

#: The template's header, in order. Must match build_track1_submission.py's HEADER; they
#: are stated twice on purpose, so a change to one is caught rather than followed.
HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
          "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]

MAX_ROWS = 10
#: The organizers' own template, when it is on disk. It is the authority for the header and
#: for the example values, and beats anything stated here.
TEMPLATE = REPO_ROOT / "docs" / "reference" / "track1_submission_template.csv"
#: The id the template uses in both of its example rows, and the builder's default. There is
#: exactly one proband in this challenge, so this is most likely the expected value rather
#: than a stand-in -- but it is the organizers' example, so a submission carrying it is
#: worth one deliberate look, not an automatic failure.
TEMPLATE_ID = "PROBAND01"
#: The values the template demonstrates in its finding_type column.
FINDING_TYPES = ("primary", "secondary")

CONTIGS = ["chr" + str(n) for n in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
BASES = re.compile(r"^[ACGTN]+$")

#: GRCh38 primary assembly lengths, for an out-of-range position. A position past the end
#: of its chromosome is almost certainly wrong -- but these are quoted from memory of the
#: assembly report rather than read from a .fai here, so an overrun is reported as a
#: WARNING and never fails the run.
GRCH38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569, "chrMT": 16569,
}

ALLELE_COLUMNS = ("chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2")


class Report(object):
    """Errors cost the submission slot; warnings need a person's judgement."""

    def __init__(self) -> None:
        self.errors = []
        self.warnings = []
        self.notes = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def note(self, message: str) -> None:
        self.notes.append(message)


def check_bytes(path: Path, report: Report) -> str:
    """Encoding and line endings, before anything tries to parse the file."""
    raw = path.read_bytes()
    if not raw.strip():
        report.error("the file is empty")
        return ""

    if raw.startswith(b"\xef\xbb\xbf"):
        report.error(
            "the file starts with a UTF-8 byte-order mark; the first column name will "
            "parse as '\\ufeffproband_id' and the row will not match the template. A "
            "spreadsheet almost certainly wrote this file -- rebuild it with the builder")
        raw = raw[3:]

    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        report.error("the file is not valid UTF-8 (%s); rebuild it" % exc.reason)
        return ""

    if "\r\n" in text:
        report.note("line endings are CRLF, which csv readers handle")
    if text.rstrip("\r\n") != text.rstrip():
        report.warn("the file ends with trailing whitespace beyond a final newline")
    return text


def check_template(report: Report) -> None:
    """If the organizers' template is on disk, it outranks this file's idea of the header."""
    if not TEMPLATE.is_file():
        report.note("template not found at %s; checking against this script's own copy of "
                    "the header" % TEMPLATE.relative_to(REPO_ROOT))
        return
    try:
        rows = list(csv.reader(io.StringIO(TEMPLATE.read_text(encoding="utf-8-sig"))))
    except (OSError, UnicodeDecodeError) as exc:
        report.warn("could not read the template (%s)" % exc)
        return
    if not rows:
        report.warn("the template is empty")
        return
    header = [cell.strip() for cell in rows[0]]
    if header != HEADER:
        report.error(
            "this script's header no longer matches the organizers' template. The "
            "template is the authority -- update HEADER in this script AND in "
            "scripts/build_track1_submission.py. Template says: %s" % ", ".join(header))
    else:
        report.note("header matches the organizers' template")


def check_header(rows: list, report: Report) -> bool:
    if not rows:
        report.error("no header row")
        return False
    header = [cell.strip() for cell in rows[0]]
    if header == HEADER:
        return True
    missing = [c for c in HEADER if c not in header]
    extra = [c for c in header if c not in HEADER]
    if missing:
        report.error("header is missing column(s): %s" % ", ".join(missing))
    if extra:
        report.error("header has unexpected column(s): %s" % ", ".join(extra))
    if not missing and not extra:
        report.error("header has the right columns in the wrong order; the template's "
                     "order is: %s" % ", ".join(HEADER))
    return False


def check_allele(row: dict, n: int, index: int, report: Report, required: bool) -> tuple:
    """One allele triple. Returns ``(contig, pos, ref, alt)`` as strings, or ``()``."""
    suffix = "_%d" % index
    contig = row["chrom" + suffix].strip()
    pos = row["pos" + suffix].strip()
    ref = row["ref" + suffix].strip().upper()
    alt = row["alt" + suffix].strip().upper()
    present = [c for c in (contig, pos, ref, alt) if c]

    if not present:
        if required:
            report.error("row %d: allele %d is empty, but every row needs a first allele"
                         % (n, index))
        return ()
    if len(present) != 4:
        report.error("row %d: allele %d is partly filled -- some of chrom/pos/ref/alt are "
                     "set and some are blank" % (n, index))
        return ()

    if not contig.startswith("chr"):
        report.error("row %d: chrom%s is not chr-prefixed. An unprefixed chromosome scores "
                     "zero while looking correct in the file" % (n, suffix))
    elif contig not in CONTIGS:
        report.error("row %d: chrom%s is not a primary GRCh38 contig" % (n, suffix))

    try:
        position = int(pos)
    except ValueError:
        report.error("row %d: pos%s is not an integer (a spreadsheet may have reformatted "
                     "it -- check for a decimal point or thousands separator)" % (n, suffix))
        position = None
    else:
        if position <= 0:
            report.error("row %d: pos%s is not positive" % (n, suffix))
        elif contig in GRCH38 and position > GRCH38[contig]:
            report.warn("row %d: pos%s is past the end of that chromosome in GRCh38"
                        % (n, suffix))

    for name, value in (("ref", ref), ("alt", alt)):
        if not BASES.match(value):
            report.error("row %d: %s%s is not a plain A/C/G/T/N sequence"
                         % (n, name, suffix))
    if ref and alt and ref == alt:
        report.error("row %d: ref%s and alt%s are identical" % (n, suffix, suffix))

    return (contig, pos, ref, alt)


def check_rows(rows: list, report: Report) -> None:
    records = [dict(zip(HEADER, row)) for row in rows[1:] if any(c.strip() for c in row)]

    if not records:
        report.error("no data rows")
        return
    if len(records) > MAX_ROWS:
        report.error("%d data rows; the template allows %d" % (len(records), MAX_ROWS))

    blank = [i for i, row in enumerate(rows[1:], start=2) if not any(c.strip() for c in row)]
    if blank:
        report.warn("blank row(s) at line %s" % ", ".join(str(i) for i in blank))

    widths = set(len(row) for row in rows[1:] if any(c.strip() for c in row))
    if widths and widths != {len(HEADER)}:
        report.error("not every row has %d fields (found widths: %s). An unescaped comma "
                     "in a notes field does this" % (len(HEADER), sorted(widths)))

    ids = set()
    epcrs = []
    signatures = []

    for n, row in enumerate(records, start=1):
        proband = row["proband_id"].strip()
        if not proband:
            report.error("row %d: proband_id is empty" % n)
        ids.add(proband)

        first = check_allele(row, n, 1, report, required=True)
        second = check_allele(row, n, 2, report, required=False)
        if first and second and first == second:
            report.error("row %d: both alleles are the same variant" % n)
        signatures.append((first, second))

        raw_epcr = row["epcr"].strip()
        try:
            epcr = float(raw_epcr)
        except ValueError:
            report.error("row %d: epcr is not a number" % n)
        else:
            epcrs.append(epcr)
            if not 0 < epcr <= 1:
                report.error("row %d: epcr is %s, outside (0, 1]" % (n, raw_epcr))

        finding = row["finding_type"].strip()
        if not finding:
            report.error("row %d: finding_type is empty" % n)
        elif finding.lower() not in FINDING_TYPES:
            report.warn("row %d: finding_type is not one of %s, the values the template "
                        "demonstrates" % (n, "/".join(FINDING_TYPES)))
        if not row["notes"].strip():
            report.warn("row %d: notes is empty. This project's limitations are supposed "
                        "to travel into this field" % n)

    if len(ids) > 1:
        report.error("rows carry %d different proband_id values; one proband per file"
                     % len(ids))
    elif ids:
        only = sorted(ids)[0]
        report.note("proband_id: %s" % only)
        if only == TEMPLATE_ID:
            report.note(
                "         ^ this is the value the organizers' own template uses in both "
                "example rows, and there is one proband in this challenge, so it is most "
                "likely correct. Confirm it once against the Space's Track 1 tab")
        else:
            report.warn(
                "proband_id is not %s, which is what the organizers' template shows in "
                "both example rows. If that is deliberate, ignore this; if not, a wrong id "
                "is rejected outright and costs a submission slot" % TEMPLATE_ID)

    if epcrs and epcrs != sorted(epcrs, reverse=True):
        report.error("rows are not in descending epcr order, so the best-supported "
                     "configuration is not at rank 1")

    seen = {}
    for n, signature in enumerate(signatures, start=1):
        if signature in seen and any(signature):
            report.error("row %d duplicates row %d" % (n, seen[signature]))
        else:
            seen[signature] = n

    report.note("rows: %d (max %d)" % (len(records), MAX_ROWS))
    if epcrs:
        report.note("epcr, in file order: %s" % ", ".join("%g" % e for e in epcrs))


def check_name(path: Path, report: Report) -> None:
    stem = path.stem
    if stem == "track1_submission":
        report.warn(
            "the file is still named '%s'. The Space asks for a name carrying your "
            "username and a short approach name -- rebuild with --name" % stem)
    if " " in stem:
        report.warn("the filename contains a space, which some upload forms reject")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path", nargs="?", default=str(DEFAULT),
                        help="the CSV to check (default: the packaged Track 1 submission)")
    args = parser.parse_args()

    path = Path(args.path)
    print("checking %s" % path)
    if not path.is_file():
        print("\nERROR  the file does not exist.")
        print("       Build it: python scripts/build_track1_submission.py --proband-id ...")
        return 1

    report = Report()
    check_template(report)
    check_name(path, report)
    text = check_bytes(path, report)

    if text:
        rows = list(csv.reader(io.StringIO(text)))
        if check_header(rows, report):
            check_rows(rows, report)

    print("  size: %d bytes" % path.stat().st_size)
    for note in report.notes:
        print("  %s" % note)

    if report.warnings:
        print("\n%d warning(s) -- your call:" % len(report.warnings))
        for message in report.warnings:
            print("  ~ %s" % message)

    if report.errors:
        print("\n%d error(s) -- these would cost the slot:" % len(report.errors))
        for message in report.errors:
            print("  x %s" % message)
        print("\nNOT READY to upload.")
        return 1

    print("\nREADY to upload." if not report.warnings
          else "\nReady to upload once the warnings above are settled.")
    print("results/ is gitignored (custody C3). Upload this file to the Space directly;")
    print("do not copy it into the repository. See decision D21.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
