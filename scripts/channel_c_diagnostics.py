"""Reproduce the measurements behind decision D19 -- why Channel C nominates nothing.

Channel C is implemented, runs, and then declines to nominate. That is an unusual outcome
to assert in a report, so the evidence for it is a script rather than a paragraph. Every
number D19 quotes is printed here.

Four measurements, in the order they were made:

  1. **The confound.** How far the proxy signature lies along the generic L1000
     drug-response axis, and what removing it does.
  2. **Random-query nulls, and why both failed.** A structured biological query beats an
     arbitrary gene set for reasons unrelated to the gene, so neither null can fail.
  3. **The control-gene null.** The same pipeline driven by *other* genes' knockdowns --
     the only comparison that isolates the question Channel C is asking.
  4. **What the head actually contains**, for the real query and for controls.

Usage::

    conda activate mva-track2
    python scripts/channel_c_diagnostics.py [--controls 40] [--random-draws 40]

Runtime is a few minutes and it needs the LINCS cache (about 6 GB, fetched on first use
by ``src.l2_channels.lincs``). It writes nothing: it reads reference data and prints.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

from src.l2_channels import (  # noqa: E402
    channel_c_signature as channel, connectivity, enrichment, lincs,
)

LOGGER = logging.getLogger("channel_c_diagnostics")

#: Compounds watched through the diagnostics. Tubulin binders are the class that keeps
#: surfacing; the others are pharmacologically implausible reversers of a mitotic
#: knockdown and are watched precisely because they surface beside the tubulin agents.
TUBULIN = {"COLCHICINE", "PACLITAXEL", "DOCETAXEL", "CABAZITAXEL", "ALBENDAZOLE",
           "IXABEPILONE", "VINBLASTINE", "VINCRISTINE"}
IMPLAUSIBLE = {"TOBRAMYCIN", "INOSITOL", "MEPIVACAINE", "ANIRACETAM", "OXANDROLONE",
               "ZALEPLON", "GLICLAZIDE", "MOLSIDOMINE"}


def load(config: dict, settings: dict):
    """Everything the measurements run on, read once."""
    paths, _ = lincs.ensure(config)
    landmarks = lincs.landmark_ids(lincs.read_metadata(paths["compound"]["genes"]))
    names, molecules, _, _ = channel._resolve_compounds(config, settings)
    kd_meta = lincs.read_metadata(paths["knockdown"]["meta"])
    gene = config["causal_gene"]

    import h5py

    with h5py.File(paths["knockdown"]["matrix"], "r") as kd, \
            h5py.File(paths["compound"]["matrix"], "r") as cp:
        kd_signature_ids, kd_genes = lincs.axes(kd)
        cp_signature_ids, cp_genes = lincs.axes(cp)
        _, kd_index, cp_index = channel._gene_axis(kd_genes, cp_genes, landmarks)
        consensus, proxy = channel.build_proxy(kd, kd_signature_ids, kd_meta, gene,
                                               kd_index, settings)
        controls = channel.control_genes(kd_meta, gene, settings, int(config["seed"]))
        control_consensus = {
            name: channel.build_proxy(kd, kd_signature_ids, kd_meta, name, kd_index,
                                      settings)[0]
            for name in controls}
        selected, _ = channel._select_signatures(
            lincs.read_metadata(paths["compound"]["meta"]), cp_signature_ids, names,
            molecules, settings)
        z_scores = channel.read_signatures(cp, selected, cp_index, settings)
    return consensus, proxy, control_consensus, selected, z_scores, molecules


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--controls", type=int, default=40,
                        help="control genes in the null (default 40)")
    parser.add_argument("--random-draws", type=int, default=40,
                        help="random queries for the superseded nulls (default 40)")
    parser.add_argument("--config", default=str(REPO_ROOT / "config" / "pipeline.yaml"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    config.setdefault("results_dir", str(REPO_ROOT / "results"))
    settings = channel.settings_for(config)
    settings["n_control_genes"] = args.controls
    gene = config["causal_gene"]

    print(f"Channel C diagnostics -- proxy gene {gene}, seed {config['seed']}\n")
    started = time.monotonic()
    consensus, proxy, control_consensus, selected, z_scores, molecules = load(config,
                                                                              settings)
    print(f"loaded in {time.monotonic() - started:.0f}s: {z_scores.shape[0]} compound "
          f"signature(s) x {z_scores.shape[1]} landmark genes\n")

    groups = channel._groups(selected)
    keys = sorted(groups)
    cells = np.array([e["cell_line"] for e in selected])
    by_name = {key: molecules.get(key, {}).get("name", "").upper() for key in keys}
    order = connectivity.descending_order(z_scores)

    # ---------------------------------------------------------------- 1. the proxy
    print("=" * 78)
    print("1. THE PROXY, AND THE CONFOUND")
    print("=" * 78)
    concordance = proxy["cross_cell_line_concordance"]
    print(f"  cell lines pooled            {len(proxy['cell_lines'])} "
          f"({', '.join(proxy['cell_lines'])})")
    print(f"  below conventional cc_q75    {proxy['n_below_conventional_cc_q75']} of "
          f"{len(proxy['signatures'])} (threshold {channel.CONVENTIONAL_CC_Q75})")
    print(f"  cross-cell-line Spearman     median {concordance['median']}, "
          f"range {concordance['min']}..{concordance['max']}")

    deflated, deflation = channel.deflate_generic_axis(consensus, z_scores)
    print(f"\n  r(proxy, mean compound signature) before  {deflation['r_before']:+.4f}")
    print(f"  r(proxy, mean compound signature) after   {deflation['r_after']:+.4f}")
    print("  -> Half the proxy is the axis every perturbation shares. Uncorrected, the")
    print("     score largely measures whether a compound provokes that response at all.")

    # ---------------------------------------------------------------- 2. random nulls
    print("\n" + "=" * 78)
    print("2. THE TWO RANDOM-QUERY NULLS, AND WHY BOTH FAIL")
    print("=" * 78)
    rng = np.random.default_rng(int(config["seed"]))
    random_best, random_zero = [], []
    for _ in range(args.random_draws):
        up, down = connectivity.random_query_masks(z_scores.shape[1],
                                                   settings["query_size"], rng)
        scores, _, wtcs = channel.score_population(z_scores, order, cells, groups, keys,
                                                   up, down)
        random_best.append(scores.min())
        random_zero.append(float((wtcs == 0).mean()))

    for label, query in (("raw", consensus), ("deflated", deflated)):
        up, down = connectivity.query_masks(query, settings["query_size"])
        scores, _, wtcs = channel.score_population(z_scores, order, cells, groups, keys,
                                                   up, down)
        p = connectivity.empirical_p(scores.min(), random_best, tail="lower")
        print(f"  {label:<9} best NCS {scores.min():+.4f}  vs random median "
              f"{np.median(random_best):+.4f}  p={p:.4f}   "
              f"molecules<0 {int((scores < 0).sum())}/{len(keys)}  "
              f"zero-WTCS {float((wtcs == 0).mean()):.1%}")
    print(f"  random queries    zero-WTCS {np.median(random_zero):.1%}")
    print("  -> Both nulls pass for a structural reason: a random query leaves far more")
    print("     signatures at exactly zero, so the observed distribution is shifted")
    print("     against it whatever gene the query came from. A null that cannot fail is")
    print("     not a null. See D19.")

    # ---------------------------------------------------------------- 3. control genes
    print("\n" + "=" * 78)
    print(f"3. THE CONTROL-GENE NULL ({len(control_consensus)} unrelated genes)")
    print("=" * 78)
    axis = channel.generic_axis(z_scores)
    up, down = connectivity.query_masks(deflated, settings["query_size"])
    observed, _, _ = channel.score_population(z_scores, order, cells, groups, keys, up,
                                              down)
    control_scores, control_heads = {}, {}
    for name, vector in sorted(control_consensus.items()):
        query = channel.project_out(vector, axis)
        c_up, c_down = connectivity.query_masks(query, settings["query_size"])
        scores, _, _ = channel.score_population(z_scores, order, cells, groups, keys,
                                                c_up, c_down)
        control_scores[name] = scores
        control_heads[name] = [by_name[keys[i]]
                               for i in np.argsort(scores, kind="stable")[:10]]

    depth, null_report = channel.calibrate_depth(observed, control_scores, settings)
    print("  depth   observed   control median   control best        p   passes")
    for rung in null_report["ladder"]:
        print(f"  {rung['depth']:>5} {rung['observed']:>10.4f} {rung['null_median']:>16.4f}"
              f" {rung['null_best']:>14.4f} {rung['p']:>8.4f}   {rung['passes']}")
    print(f"\n  nomination depth: {depth}  (0 means the channel nominates nothing)")

    # ---------------------------------------------------------------- 4. the head
    print("\n" + "=" * 78)
    print("4. WHAT THE HEAD CONTAINS")
    print("=" * 78)
    head = [by_name[keys[i]] for i in np.argsort(observed, kind="stable")[:10]]
    print(f"  {gene} top 10: {', '.join(head)}")
    print(f"    tubulin binders {sum(1 for n in head if n in TUBULIN)}, "
          f"implausible reversers {sum(1 for n in head if n in IMPLAUSIBLE)}")
    print("\n  the same for five control genes:")
    for name in sorted(control_heads)[:5]:
        names = control_heads[name]
        print(f"    {name:<10} {', '.join(names[:4])}"
              f"   (tubulin {sum(1 for n in names if n in TUBULIN)})")
    shared = sum(1 for name, names in control_heads.items()
                 if any(n in TUBULIN for n in names))
    print(f"\n  {shared} of {len(control_heads)} control genes also surface a tubulin "
          "binder in their top 10.")
    print("  -> The tubulin agents are extreme signatures, not a BUB1B-specific finding.")
    print(f"\ntotal {time.monotonic() - started:.0f}s")


if __name__ == "__main__":
    main()
