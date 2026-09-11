"""Compare snpEff transcript policies for loss-of-function calling on the SAC panel.

Purpose
    Evidence for L0's transcript-policy decision (see
    ``docs/research/transcript-policy-lof.md``). For each SAC-panel gene, and for each of
    three snpEff policies -- every transcript, ``-canon``, ``-tag MANE_Select`` -- report:

    1. which transcript(s) the policy keeps, and whether that is the MANE Select transcript;
    2. how the policy classifies the gene's public ClinVar P/LP variants;
    3. how many base pairs the policy makes LoF-eligible (CDS or canonical splice site), and
       where the policies disagree -- the territory in which a *novel* variant's call would
       depend on the policy, which (2) cannot measure because ClinVar variants are mostly
       described on the clinical transcript already.

Inputs
    ``config/pipeline.yaml`` (``annotator.database``, ``annotator.java_heap``,
    ``results_dir``); the panel from ``src/l0_genomics/run.py``; NCBI's current MANE summary
    table and ClinVar's GRCh38 VCF, both public and streamed.

Outputs
    ``<results_dir>/transcript_policy/``: ``picks.tsv``, ``clinvar_variants.tsv``,
    ``clinvar_summary.tsv``, ``territory.tsv``, ``manifest.json`` (tool, database and source
    versions, so every number can be traced to what produced it). Deterministic: no RNG.
    Runs in WSL from the repo root, ``python scripts/transcript_policy_check.py``; about
    16 minutes on a 12-core machine, most of it the per-position territory scan.

Guardrail
    **No patient data.** Nothing here reads ``data_dir``. Every variant is either a public
    ClinVar record or an invented substitution, and reaches snpEff on stdin, so no ``.vcf``
    file is written anywhere. snpEff runs with ``-noLog -nodownload -noStats`` from a
    throwaway directory. Hypothesis-generation tooling only; nothing here is a clinical
    variant classification.
"""

from __future__ import annotations

import collections
import gzip
import json
import re
import subprocess
import sys
import tempfile
import threading
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.l0_genomics.run import SAC_PANEL  # noqa: E402
from src.pipeline import load_config  # noqa: E402

MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
UA = {"User-Agent": "mva-track2-transcript-policy/1.0 (https://github.com/KFlowerz/HKTN_RareDisease_MVA)"}

POLICIES = {"all": [], "canon": ["-canon"], "mane": ["-tag", "MANE_Select"]}
RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
# LoF classes after MacArthur et al. -- nonsense, frameshift, canonical splice
# [macarthur2012] doi:10.1126/science.1215040. Read off HGVS.c so the classification does
# not depend on an invented REF base matching the genome: a bare c.N is CDS, c.N+1/+2 and
# c.N-1/-2 are the canonical donor/acceptor dinucleotides; c.-N / c.*N (UTR) never match.
LOF_ELIGIBLE = re.compile(r"^c\.\d+([+-][12])?[A-Z]>")
VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _fetch(url: str):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120)


def snpeff(cfg: dict, command: str, *args: str, stdin=None, cwd: str):
    """Yield snpEff stdout lines. ``stdin`` is an iterable of text, fed from a thread."""
    argv = ["snpEff", f"-Xmx{cfg['annotator']['java_heap']}", command,
            "-noLog", "-nodownload", *args]
    proc = subprocess.Popen(argv, stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=cwd)
    if stdin is not None:
        def feed():
            for chunk in stdin:
                proc.stdin.write(chunk)
            proc.stdin.close()
        threading.Thread(target=feed, daemon=True).start()
    yield from proc.stdout
    if proc.wait():
        raise RuntimeError(f"snpEff {command} exited {proc.returncode}")


def annotations(info: str):
    """Split a VCF INFO field's ANN into per-transcript field lists."""
    ann = next((kv[4:] for kv in info.split(";") if kv.startswith("ANN=")), "")
    return [a.split("|") for a in ann.split(",") if a]


def gene_coords(cfg: dict, panel, cwd: str) -> dict:
    coords = {}
    for line in snpeff(cfg, "genes2bed", cfg["annotator"]["database"], *panel, cwd=cwd):
        if line.startswith("#"):
            continue
        chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
        coords[name.split(";")[0]] = (chrom, int(start), int(end))
    return coords


def mane_table(panel) -> tuple:
    listing = _fetch(MANE_DIR).read().decode()
    name = re.search(r"MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz", listing).group(0)
    rows = collections.defaultdict(list)
    with _fetch(MANE_DIR + name) as resp, gzip.open(resp, "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["symbol"] in panel:
                rows[r["symbol"]].append(r)
    return name, rows


def clinvar_plp(panel) -> tuple:
    """Stream ClinVar; keep P/LP records whose GENEINFO names a panel gene."""
    quick = re.compile(r"GENEINFO=(?:[^;\t]*\|)?(" + "|".join(panel) + r"):")
    date, keep = None, []
    with _fetch(CLINVAR_VCF) as resp, gzip.open(resp, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##fileDate="):
                    date = line.strip().split("=", 1)[1]
                continue
            m = quick.search(line)
            if not m:
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            sig = info.get("CLNSIG", "")
            if not sig.startswith(("Pathogenic", "Likely_pathogenic")) or f[4] in (".", ""):
                continue
            mc = ",".join(x.split("|")[1] for x in info.get("MC", "").split(",") if "|" in x)
            keep.append(dict(id=f[2], chrom=f[0], pos=f[1], ref=f[3], alt=f[4], gene=m.group(1),
                             clnsig=sig, clnvc=info.get("CLNVC", ""),
                             review=info.get("CLNREVSTAT", ""), mc=mc))
    return date, keep


def policy_picks(cfg, coords, cwd) -> dict:
    """Transcripts each policy reports at an invented SNV at each gene's midpoint."""
    lines = [VCF_HEADER] + [f"{c}\t{(s + e) // 2}\t{g}\tA\tC\t.\tPASS\t.\n" for g, (c, s, e) in coords.items()]
    picks = {}
    for pol, flags in POLICIES.items():
        picks[pol] = collections.defaultdict(set)
        for line in snpeff(cfg, "ann", "-noStats", *flags, cfg["annotator"]["database"], "-",
                           stdin=lines, cwd=cwd):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            for q in annotations(f[7]):
                if len(q) > 6 and q[3] == f[2] and q[5] == "transcript":
                    picks[pol][f[2]].add(q[6])
    return picks


def clinvar_by_policy(cfg, variants, cwd) -> dict:
    """Most severe consequence on the variant's own gene, per policy."""
    lines = [VCF_HEADER] + [f"{v['chrom']}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}\t.\tPASS\t.\n"
                            for v in variants]
    gene_of = {v["id"]: v["gene"] for v in variants}
    out = collections.defaultdict(dict)
    for pol, flags in POLICIES.items():
        for line in snpeff(cfg, "ann", "-noStats", *flags, cfg["annotator"]["database"], "-",
                           stdin=lines, cwd=cwd):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            best = (-1, "no_transcript", "")
            for q in annotations(f[7]):
                if len(q) > 6 and q[3] == gene_of[f[2]] and q[5] == "transcript" and RANK.get(q[2], -1) > best[0]:
                    best = (RANK[q[2]], q[1], q[6])
            out[f[2]][pol] = best
    return out


def territory(cfg, coords, cwd) -> dict:
    """LoF-eligible positions per policy, from an invented SNV at every position of each span."""
    def lines():
        yield VCF_HEADER
        for g, (c, s, e) in coords.items():
            for pos in range(s, e + 1):     # inclusive: covers either BED convention
                yield f"{c}\t{pos}\t{g}\tA\tC\t.\tPASS\t.\n"

    elig = {p: collections.defaultdict(set) for p in POLICIES}
    for pol, flags in POLICIES.items():
        for line in snpeff(cfg, "ann", "-noStats", "-noLof", "-no-upstream", "-no-downstream",
                           "-no-intergenic", *flags, cfg["annotator"]["database"], "-",
                           stdin=lines(), cwd=cwd):
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            for q in annotations(f[7]):
                if len(q) > 9 and q[3] == f[2] and q[7] == "protein_coding" and LOF_ELIGIBLE.match(q[9]):
                    elig[pol][f[2]].add(int(f[1]))
                    break
    return elig


def write_tsv(path: Path, header: list, rows) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def main() -> int:
    cfg = load_config()
    out = Path(cfg["results_dir"]) / "transcript_policy"
    out.mkdir(parents=True, exist_ok=True)
    panel = tuple(SAC_PANEL)
    work = tempfile.mkdtemp(prefix="snpeff-")
    inv = {v: k for k, v in RANK.items()}

    version = subprocess.run(["snpEff", "-version"], capture_output=True, text=True).stdout.strip()
    coords = gene_coords(cfg, panel, work)
    if set(panel) - coords.keys():
        raise SystemExit(f"genes missing from {cfg['annotator']['database']}: {set(panel) - coords.keys()}")

    # 1. What each policy keeps, against MANE ------------------------------------------
    mane_file, mane = mane_table(panel)
    picks = policy_picks(cfg, coords, work)
    rows = []
    print(f"\n[1] transcript kept per policy (MANE: {mane_file})")
    for g in panel:
        select = {r["Ensembl_nuc"] for r in mane[g] if r["MANE_status"] == "MANE Select"}
        plus = {r["Ensembl_nuc"] for r in mane[g] if r["MANE_status"] == "MANE Plus Clinical"}
        for pol in POLICIES:
            kept = picks[pol][g]
            rows.append([g, pol, len(kept), ",".join(sorted(kept)), ",".join(sorted(select)),
                         ",".join(sorted(plus)) or "-", bool(select & kept)])
        print(f"  {g:7s} MANE={','.join(select)}  canon={','.join(picks['canon'][g])}  "
              f"mane-tag={','.join(picks['mane'][g])}  plus_clinical={','.join(plus) or '-'}")
    write_tsv(out / "picks.tsv", ["gene", "policy", "n_transcripts", "transcripts", "mane_select",
                                  "mane_plus_clinical", "includes_mane_select"], rows)

    # 2. Public ClinVar P/LP variants under each policy -------------------------------
    cv_date, variants = clinvar_plp(panel)
    cv = clinvar_by_policy(cfg, variants, work)
    rows, summ = [], collections.defaultdict(collections.Counter)
    for v in sorted(variants, key=lambda v: (v["gene"], int(v["id"]))):
        r = cv[v["id"]]
        rows.append([v["id"], v["gene"], v["clnsig"], v["clnvc"], v["review"], v["mc"]] +
                    [x for p in POLICIES for x in (inv.get(r[p][0], "NONE"), r[p][1], r[p][2])])
        s = summ[v["gene"]]
        s["n"] += 1
        for p in POLICIES:
            s[f"high_{p}"] += r[p][0] == 3
            s[f"moderate_{p}"] += r[p][0] == 2
        s["discordant"] += len({r[p][0] for p in POLICIES}) > 1
    write_tsv(out / "clinvar_variants.tsv",
              ["clinvar_id", "gene", "clnsig", "clnvc", "review", "clinvar_mc"] +
              [f"{p}_{k}" for p in POLICIES for k in ("impact", "effect", "transcript")], rows)
    cols = ["n"] + [f"high_{p}" for p in POLICIES] + [f"moderate_{p}" for p in POLICIES] + ["discordant"]
    write_tsv(out / "clinvar_summary.tsv", ["gene"] + cols,
              [[g] + [summ[g][c] for c in cols] for g in panel])
    print(f"\n[2] ClinVar P/LP (fileDate {cv_date}): gene n | HIGH all/canon/mane | discordant")
    for g in panel:
        s = summ[g]
        print(f"  {g:7s} {s['n']:3d} | {s['high_all']:3d}/{s['high_canon']:3d}/{s['high_mane']:3d} | {s['discordant']}")

    # 3. LoF-eligible territory per policy ---------------------------------------------
    elig = territory(cfg, coords, work)
    rows = []
    print("\n[3] LoF-eligible bp: gene all/canon/mane | all-not-mane canon-not-mane mane-not-canon")
    for g in panel:
        a, c, m = elig["all"][g], elig["canon"][g], elig["mane"][g]
        rows.append([g, len(a), len(c), len(m), len(a - m), len(c - m), len(m - c)])
        print(f"  {g:7s} {len(a):5d}/{len(c):5d}/{len(m):5d} | {len(a - m):4d} {len(c - m):4d} {len(m - c):4d}")
    write_tsv(out / "territory.tsv", ["gene", "bp_all", "bp_canon", "bp_mane", "bp_all_not_mane",
                                      "bp_canon_not_mane", "bp_mane_not_canon"], rows)

    try:
        sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True,
                             text=True).stdout.strip() or None
    except OSError:
        sha = None
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    script="scripts/transcript_policy_check.py", git_sha=sha, seed=cfg["seed"],
                    snpeff=version, database=cfg["annotator"]["database"],
                    java_heap=cfg["annotator"]["java_heap"], policies=POLICIES, panel=panel,
                    mane_summary=mane_file, clinvar_url=CLINVAR_VCF, clinvar_fileDate=cv_date,
                    n_clinvar_plp=len(variants), snpeff_cwd_leftovers=len(list(Path(work).iterdir())))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nwritten: {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
