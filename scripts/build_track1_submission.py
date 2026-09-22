"""Build the Track 1 submission CSV from L0's recorded causal-gene call.

The submission format is the subject's variant coordinates (decision D21), so this writes
**under ``results/``** — gitignored, custody location ``C3``, inside the purge scope — and
nowhere else. It prints counts and checks, never a coordinate: a session transcript is
custody location ``C7`` and does not need them in it.

What it emits, per the Space's template:

    proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes

Ranking strategy, checked against the organizers' own scorer rather than assumed:

* **Row 1 is the compound-heterozygous pair** L0 called in the causal gene. A full match at
  rank 1 scores the maximum 100 rank points and an F-max of 1.0.
* **Row 2 is the truncating allele alone, and it is scoring-neutral.** It was added as a
  hedge against the second allele being something a called VCF cannot see — a copy-number,
  structural or deep-intronic variant, which D9 records as an open caveat — and measurement
  showed it changes nothing: ``score_proband`` already awards partial credit from the pair
  row itself, because it tests ``row.variants & true_variants``. The row is kept for the
  human reviewer, whose judging is separate from the automated score, and its ``notes``
  field states the failure mode it stands for. It is not kept because it helps the metric.

**Pair first is the right order, and the alternative was measured.** Leading with the
truncating allele alone would raise F-max from 0.500 to 0.667 in the fallback scenario and
cost 50 rank points in the main one, because the full match would fall to rank 2. Rank
points run 0–100 and F-max 0–1, so that trade is clearly bad.

Usage::

    python scripts/build_track1_submission.py [--proband-id PROBAND01]

Run the pipeline's L0 first; this reads its artifacts and computes nothing.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "results"
L0 = RESULTS / "l0_genomics"
OUT_DIR = RESULTS / "track1"

#: The template's header, in order. Column names are the Space's, not ours.
HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
          "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]

#: Ranked first: the best-supported configuration. Not 0.95 or higher -- the pair is
#: unphased and the missense allele is untested (D9), and the number should say so.
EPCR_PAIR = 0.70
#: The hedge. Lower, because on its own it is an incomplete answer for a recessive disease.
EPCR_TRUNCATING = 0.45

PAIR_NOTE = (
    "Biallelic BUB1B, the reported MVA configuration: a MANE Select truncating allele "
    "(ClinVar pathogenic, two-star) paired with a rare kinase-domain missense. UNPHASED - "
    "no parental samples, so cis cannot be excluded; this is the largest uncertainty and "
    "no available data closes it. The missense allele is untested in vitro."
)
TRUNCATING_NOTE = (
    "Hedge: the truncating allele alone, in case the second hit is a copy-number, "
    "structural or deep-intronic variant that a called VCF cannot surface."
)


def load(path: Path):
    if not path.is_file():
        raise SystemExit(
            f"{path} not found. Run L0 first: "
            "PYTHONHASHSEED=42 python -m src.pipeline --only l0_genomics")
    return json.loads(path.read_text(encoding="utf-8"))


def genotype_index(calls: list) -> dict:
    """``{variant_id: genotype record}`` -- where contig, pos, ref and alt live."""
    return {call["genotype"]["variant_id"]: call["genotype"] for call in calls
            if call.get("genotype", {}).get("variant_id")}


def annotation_index(calls: list) -> dict:
    return {call["annotation"]["variant_id"]: call["annotation"] for call in calls
            if call.get("annotation", {}).get("variant_id")}


def with_chr(contig: str) -> str:
    """``7`` -> ``chr7``. The dataset's VCF uses bare contig names; the template does not.

    Asserted rather than assumed: an unprefixed chromosome scores zero while looking
    entirely correct in the file, which is the kind of mistake that is only discovered
    after the submission slot is spent.
    """
    contig = str(contig).strip()
    return contig if contig.startswith("chr") else f"chr{contig}"


def build_rows(proband_id: str) -> tuple:
    """``(rows, summary)`` -- the CSV rows and what may safely be printed about them."""
    call = load(L0 / "causal_gene_call.json")
    calls = load(L0 / "variant_calls.json")["calls"]
    genotypes = genotype_index(calls)
    annotations = annotation_index(calls)

    candidates = call.get("candidate_genes") or []
    if len(candidates) != 1:
        raise SystemExit(
            f"L0 proposes {len(candidates)} candidate gene(s); this builder encodes the "
            "single-gene compound-het case recorded at gate G1 (D9). Re-check the call "
            "before submitting.")
    gene = candidates[0]

    configurations = call["per_gene"][gene]["configurations"]
    if len(configurations) != 1:
        raise SystemExit(f"{gene} has {len(configurations)} configurations; expected one.")
    config = configurations[0]
    variant_ids = config["variant_ids"]
    if len(variant_ids) != 2:
        raise SystemExit(
            f"expected a two-variant compound-het configuration, got {len(variant_ids)}.")

    # Which of the pair is the truncating one? The hedge row needs it specifically.
    tiers = dict(zip(variant_ids, config["lof_tiers"]))
    effects = dict(zip(variant_ids, config["effect_classes"]))
    truncating = [v for v in variant_ids if effects.get(v) == "lof"]
    if len(truncating) != 1:
        raise SystemExit(
            f"expected exactly one loss-of-function allele in the pair, got "
            f"{len(truncating)}. Effect classes: {config['effect_classes']}.")
    lof_id = truncating[0]

    def fields(variant_id: str) -> tuple:
        genotype = genotypes.get(variant_id)
        if genotype is None:
            raise SystemExit(f"{variant_id} is not in variant_calls.json")
        return (with_chr(genotype["contig"]), int(genotype["pos"]),
                str(genotype["ref"]).upper(), str(genotype["alt"]).upper())

    first, second = (fields(variant_ids[0]), fields(variant_ids[1]))
    lof_chrom, lof_pos, lof_ref, lof_alt = fields(lof_id)

    rows = [
        {"proband_id": proband_id,
         "chrom_1": first[0], "pos_1": first[1], "ref_1": first[2], "alt_1": first[3],
         "chrom_2": second[0], "pos_2": second[1], "ref_2": second[2], "alt_2": second[3],
         "epcr": EPCR_PAIR, "finding_type": "primary", "notes": PAIR_NOTE},
        {"proband_id": proband_id,
         "chrom_1": lof_chrom, "pos_1": lof_pos, "ref_1": lof_ref, "alt_1": lof_alt,
         "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
         "epcr": EPCR_TRUNCATING, "finding_type": "primary", "notes": TRUNCATING_NOTE},
    ]

    summary = {
        "gene": gene,
        "configuration_class": config["class"],
        "phase": config["phase"],
        "effect_classes": config["effect_classes"],
        "lof_tiers": [tiers[v] for v in variant_ids],
        "flags": config["flags"],
        "n_rows": len(rows),
        "chr_prefixed": all(str(r["chrom_1"]).startswith("chr") for r in rows),
        "epcr_descending": [r["epcr"] for r in rows] == sorted(
            (r["epcr"] for r in rows), reverse=True),
        "distinct_variants": len({(first), (second)}),
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--proband-id", default="PROBAND01",
                        help="as provided in the dataset (default: PROBAND01)")
    parser.add_argument("--name", default="track1_submission",
                        help="basename; the Space asks for username and approach in it")
    args = parser.parse_args()

    rows, summary = build_rows(args.proband_id)

    for check, message in (
        (summary["chr_prefixed"], "chromosomes are not chr-prefixed"),
        (summary["epcr_descending"], "rows are not in descending EPCR order"),
        (summary["distinct_variants"] == 2, "the pair's two variants are not distinct"),
        (len(rows) <= 10, "more than 10 rows; the template allows 10"),
        (all(0 < r["epcr"] <= 1 for r in rows), "an EPCR is outside (0, 1]"),
    ):
        if not check:
            raise SystemExit(f"refusing to write the submission: {message}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{args.name}.csv"
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)

    # Deliberately no coordinates in this output -- see the module docstring.
    print(f"wrote {out.relative_to(REPO_ROOT)}  ({out.stat().st_size} bytes)")
    print("  gene               :", summary["gene"])
    print("  configuration      :", summary["configuration_class"])
    print("  phase              :", summary["phase"])
    print("  effect classes     :", summary["effect_classes"])
    print("  LoF tiers          :", summary["lof_tiers"])
    print("  flags              :", summary["flags"])
    print("  rows               :", summary["n_rows"], "(max 10)")
    print("  chr-prefixed       :", summary["chr_prefixed"])
    print("  EPCR descending    :", summary["epcr_descending"])
    print()
    print("results/ is gitignored (custody C3). Do not copy this file into the repository;")
    print("upload it to the Space directly. See decision D21.")


if __name__ == "__main__":
    main()
