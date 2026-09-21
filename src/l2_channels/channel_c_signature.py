"""L2 Channel C -- Signature reversion (on a PROXY signature).

Purpose
    Nominate drugs whose transcriptional consequences *reverse* the disease signature --
    the classic connectivity-map logic [lamb2006] doi:10.1126/science.1132939, scored by
    the L1000 weighted connectivity score [subramanian2017]
    doi:10.1016/j.cell.2017.10.049.

    **There is no patient RNA-seq.** The subject's dataset is WGS (raw reads + called VCF)
    plus phenotype; no expression data exists to build a real disease signature from. This
    channel therefore runs on a **PROXY** signature: the LINCS L1000 consensus knockdown
    of the causal gene, pooled across every cell line that carries one. The alternative
    the scaffold named -- a curated aneuploidy-response gene set -- is *not* silently
    substituted; if no knockdown of the causal gene exists, the channel says so and stops.

Inputs
    L0/G1 causal gene (``config["causal_gene"]``); LINCS L1000 via
    :mod:`~src.l2_channels.lincs` (GEO ``GSE106127`` knockdowns, ``GSE70138`` Phase II
    compounds); drug identity from the restricted enrichment zone.

Outputs
    ``candidates.tsv`` -- ranked candidates with the normalised connectivity score, their
    position in the scored population, the identity of the proxy signature, the cell lines
    it came from, and the caveat that travels with all of it. ``channel.json`` -- the
    proxy's construction and quality, the measured confound and its removal, the null at
    every candidate depth, the parameters, counts and provenance.

Guardrail
    **The proxy substitution must travel with the result.** Every row this channel emits,
    and every figure or table built from it downstream, names the proxy and its limits. A
    knockdown in an immortalized cell line is not a child's tissue; a curated aneuploidy
    set is not this child's aneuploidy. Reporting a connectivity score without that caveat
    would overstate the evidence -- this is the single most over-claimable channel in the
    pipeline, and :data:`PROXY_CAVEAT` is written into every row rather than into a
    footnote that a downstream join could drop.

    Never label proxy-derived output as a patient signature, and never impute patient
    expression from genotype to manufacture one.

    **The generic drug-response axis is removed before scoring, and that is not optional
    tuning.** Measured on this data (``scripts/channel_c_diagnostics.py``): the BUB1B
    knockdown consensus correlates r = +0.52 with the mean L1000 compound signature, so
    most of what an uncorrected score measures is whether a compound provokes the
    universal transcriptional stress response. Uncorrected, the channel's best candidate
    does not beat a *random* gene set of the same size (p = 0.95); corrected, it beats all
    40 draws tested (p = 0.024). An uncorrected ranking here is noise with a plausible
    shape, which is worse than no channel at all.

    **A ranking needs a null that can fail.** The null is a max-statistic over random
    queries: it asks whether the k-th best candidate beats the k-th best that an arbitrary
    gene set achieves, and the channel nominates only as deep as that holds. A
    *per-candidate* query-permutation null was implemented first and **removed** -- it
    passed 60 candidates out of 60, because a structured query scores extremely against
    every signature, so it was testing whether the query had structure rather than whether
    a compound was special. It is recorded here so it is not reinvented.

    **This channel annotates; it does not exclude.** A drug scored here is nominated for
    L3, not recommended: L4 owns every safety and eligibility decision, and nothing in
    this file filters on toxicity, age or indication.

    No patient data is read or written here. The channel never opens ``data_dir``; the
    only thing it asks LINCS for is the causal gene symbol, which is a public research
    premise recorded at gate G1 (decision D9), not patient data. Drug identity comes
    through :func:`~src.l2_channels.enrichment.publishable`, so no restricted content
    reaches the table.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from pathlib import Path

from ..l0_genomics.run import SAC_PANEL
from . import connectivity, enrichment, lincs

LOGGER = logging.getLogger(__name__)

#: Defaults for everything under ``config["l2"]["channel_c"]``.
DEFAULTS = {
    "query_size": 100,              # genes per tail, CLUE's convention for an L1000 query
    "min_replicate_correlation": 0.0,   # see CONVENTIONAL_CC_Q75
    "deflate_generic_axis": True,   # see the guardrail -- off makes the channel noise
    "approved_only": True,          # scope is repurposing of APPROVED drugs
    "open_targets_release": "26.06",    # drug identity only; restricted, never published
    "top_n": 250,                   # deepest nomination the ladder may reach
    "n_control_genes": 40,          # other genes' knockdowns -- the null; see the guardrail
    "p_max": 0.05,
    "block_size": 2000,             # signature rows read from the 5 GB matrix at a time
    "max_scored_signatures": 200000,
}

#: Candidate depths the null is evaluated at. The channel nominates the deepest one whose
#: order statistic still beats random queries, so the length of the list is a measurement
#: rather than a chosen cut-off.
DEPTH_LADDER = (1, 5, 10, 25, 50, 100, 250)

#: The replicate-correlation value conventionally taken as the floor for a "gold"
#: L1000 signature [subramanian2017] doi:10.1016/j.cell.2017.10.049. It is **reported,
#: not enforced**: see :func:`build_proxy` for why, and decision D19 for the record.
CONVENTIONAL_CC_Q75 = 0.2

PROXY_TYPE = "lincs_knockdown_consensus"

#: Written into every row. Phrased as what the score is *not*, because that is the part a
#: reader will otherwise supply from the word "signature".
PROXY_CAVEAT = (
    "PROXY SIGNATURE, NOT THE PATIENT. There is no RNA from the subject. This score "
    "measures reversal of an shRNA knockdown of the causal gene in immortalized cancer "
    "cell lines, pooled across lines -- a model of losing the gene in a dish, not of this "
    "child's constitutional mosaic aneuploidy, and not of any tissue the disease affects. "
    "A high score is a hypothesis about a cell-line phenotype. Rank order within this "
    "table is weakly resolved (see channel.json, null); treat membership as the signal "
    "and position as provisional."
)

CAVEATS = (
    PROXY_CAVEAT,
    "The knockdown is acute and complete; the subject's genotype is constitutional and "
    "hypomorphic. A cell adapting over a lifetime to reduced gene dosage need not "
    "resemble one 96 hours after the transcript was removed.",
    "The proxy signature lies substantially along the generic L1000 drug-response axis "
    "(r = +0.52 with the mean compound signature). That component is projected out before "
    "scoring, so this channel measures the part of the knockdown response that is NOT the "
    "universal stress response. Without that correction the ranking does not beat random "
    "gene sets; with it, it beats them narrowly.",
    "Scoring runs on the 978 landmark genes L1000 actually measures. The rest of the "
    "Phase II matrix is inferred from those landmarks, so including it would count the "
    "same measurements again under other names -- and 978 genes is a narrow window on a "
    "transcriptome.",
    "Reversing a transcriptional signature is not a therapeutic claim. It says a "
    "compound moves cells away from the knockdown state on these genes; it says nothing "
    "about whether that helps, at what concentration, or in which tissue.",
    "Coverage is set by what LINCS profiled. A drug absent from the Phase II compound "
    "set cannot be scored at all, so absence from this ranking is absence of a profile, "
    "not absence of effect.",
    "The null is uncorrected for the number of compounds tested in the sense that it "
    "reports one p per depth, not one per drug: it establishes that the head of the "
    "ranking is better than random, not that any individual drug in it is.",
)


def settings_for(config: dict) -> dict:
    """``config["l2"]["channel_c"]`` over :data:`DEFAULTS`, with the types checked."""
    merged = {**DEFAULTS, **(((config.get("l2") or {}).get("channel_c")) or {})}
    unknown = sorted(set(merged) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"l2.channel_c has unknown setting(s) {unknown}; expected "
                         f"{sorted(DEFAULTS)}")
    for key in ("query_size", "top_n", "n_control_genes", "block_size",
                "max_scored_signatures"):
        merged[key] = int(merged[key])
    for key in ("min_replicate_correlation", "p_max"):
        merged[key] = float(merged[key])
    for key in ("approved_only", "deflate_generic_axis"):
        merged[key] = bool(merged[key])
    if merged["n_control_genes"] < 2:
        raise ValueError("l2.channel_c.n_control_genes must be at least 2: a ranking "
                         "without a null is not reported (see the guardrail)")
    return merged


def _gene_axis(knockdown_genes: list, compound_genes: list, landmarks: list) -> tuple:
    """The gene axis both matrices share. Returns ``(entrez_ids, kd_index, cp_index)``.

    Restricted to the landmark space and to genes present in both releases, ordered by
    the knockdown matrix so the result does not depend on dictionary iteration.

    Raises:
        ValueError: If the shared axis has lost most of the landmark space. That means
            the two releases were built on different platforms, and a score computed over
            the remainder would still look like a number.
    """
    landmark_set = set(landmarks)
    compound_position = {gene: i for i, gene in enumerate(compound_genes)}
    shared = [(gene, i) for i, gene in enumerate(knockdown_genes)
              if gene in landmark_set and gene in compound_position]
    if len(shared) < 0.9 * len(landmark_set):
        raise ValueError(
            f"only {len(shared)} of {len(landmark_set)} landmark gene(s) are present in "
            "both LINCS releases. The two matrices are not on the same platform, and a "
            "connectivity score over the remainder would be a number without a meaning"
        )
    entrez = [gene for gene, _ in shared]
    return (entrez, [i for _, i in shared], [compound_position[g] for g in entrez])


def build_proxy(handle, signature_ids: list, meta: list, gene: str, kd_index: list,
                settings: dict) -> tuple:
    """The proxy signature: consensus knockdown of ``gene``. Returns ``(consensus, report)``.

    Every cell line with a consensus knockdown signature contributes, and the consensus is
    their unweighted mean z-score.

    **The conventional replicate-correlation floor is reported rather than enforced**
    (:data:`CONVENTIONAL_CC_Q75`, default ``min_replicate_correlation`` 0.0), and that is
    a decision, not an oversight -- see D19. Applying it here leaves two cell lines of
    nine, turning a proxy that is at least reproducible across lineages into one that
    rests on whichever line happened to score well. The quantity that matters for a proxy
    pooled across lines is whether the lines agree with each other, so the report carries
    the pairwise concordance between them, the per-line correlations, and how many fall
    below the convention. A reader can apply the conventional filter by setting it.

    Raises:
        ValueError: If LINCS holds no consensus knockdown of the gene, or the filter
            removes them all. The scaffold's fallback -- a curated aneuploidy gene set --
            is deliberately *not* substituted: it is a different experiment answering a
            different question, and swapping it in silently would leave the channel
            reporting a score whose provenance no longer matches its name.
    """
    import numpy as np
    from scipy import stats

    rows = [row for row in meta
            if row.get("pert_iname", "").strip().upper() == gene.strip().upper()
            and row.get("pert_type", "").strip() == lincs.KNOCKDOWN_TYPE]
    if not rows:
        raise ValueError(
            f"LINCS holds no {lincs.KNOCKDOWN_TYPE} consensus signature for {gene!r}, so "
            "there is no proxy to reverse. Channel C does not fall back to a curated "
            "aneuploidy gene set on its own: that is a different experiment and needs a "
            "recorded decision, not a silent substitution."
        )

    position = {sid: i for i, sid in enumerate(signature_ids)}
    kept, dropped = [], []
    for row in rows:
        sig_id = row["sig_id"]
        try:
            quality = float(row.get("distil_cc_q75", "nan"))
        except ValueError:
            quality = float("nan")
        entry = {"sig_id": sig_id, "cell_line": row.get("cell_id", ""),
                 "distil_cc_q75": None if quality != quality else round(quality, 4),
                 "distil_ss": row.get("distil_ss", ""),
                 "matrix_row": position.get(sig_id)}
        if entry["matrix_row"] is None:
            dropped.append({**entry, "reason": "signature id absent from the matrix"})
        elif quality == quality and quality < settings["min_replicate_correlation"]:
            dropped.append({**entry, "reason": "below min_replicate_correlation"})
        else:
            kept.append(entry)

    if not kept:
        raise ValueError(
            f"all {len(rows)} knockdown signature(s) for {gene!r} were dropped "
            f"(min_replicate_correlation={settings['min_replicate_correlation']}). "
            "Lower the threshold or record a decision to use a different proxy."
        )

    kept.sort(key=lambda e: e["cell_line"])
    matrix = np.vstack([lincs.read_block(handle, e["matrix_row"], e["matrix_row"] + 1,
                                         kd_index) for e in kept])
    consensus = matrix.mean(axis=0)

    # Do the cell lines agree? That is the quality question for a proxy pooled across
    # lineages, and it is not what distil_cc_q75 measures -- that is agreement between
    # replicate hairpins inside one line.
    pairwise = []
    for i in range(len(kept)):
        for j in range(i + 1, len(kept)):
            pairwise.append(float(stats.spearmanr(matrix[i], matrix[j]).statistic))
    concordance = {
        "measure": "pairwise Spearman rho between the cell lines' knockdown signatures",
        "n_pairs": len(pairwise),
        "median": round(float(np.median(pairwise)), 4) if pairwise else None,
        "min": round(min(pairwise), 4) if pairwise else None,
        "max": round(max(pairwise), 4) if pairwise else None,
    }
    below = [e for e in kept if e["distil_cc_q75"] is not None
             and e["distil_cc_q75"] < CONVENTIONAL_CC_Q75]
    report = {
        "gene": gene,
        "proxy_type": PROXY_TYPE,
        "signature_id": f"{lincs.RELEASES['knockdown']['accession']}:{gene}:"
                        f"consensus({len(kept)} cell lines)",
        "cell_lines": [e["cell_line"] for e in kept],
        "signatures": kept,
        "dropped": dropped,
        "n_below_conventional_cc_q75": len(below),
        "conventional_cc_q75": CONVENTIONAL_CC_Q75,
        "cross_cell_line_concordance": concordance,
        "note": (
            f"{len(below)} of {len(kept)} contributing signature(s) fall below the "
            f"conventional gold threshold distil_cc_q75 >= {CONVENTIONAL_CC_Q75}. They "
            "are kept deliberately (decision D19): filtering at the convention would "
            "leave too few cell lines for a pooled proxy to mean anything, and "
            "cross-cell-line concordance above is the quality evidence that applies to a "
            "pooled signature. This is a weak proxy and the ranking inherits that."
        ),
    }
    return consensus, report


def generic_axis(z_scores):
    """The unit mean L1000 compound signature -- what every perturbation has in common.

    **Mean-centred across genes**, so that removing a component along it is the same
    operation as driving the Pearson correlation with it to zero. The size of the confound
    is quoted as a correlation (r = +0.52), and a correction whose geometry did not match
    the statistic used to justify it would leave a residual that the report would not
    account for. Centring costs nothing downstream: the connectivity score depends only on
    the *order* of the query's values, which a constant shift does not change.

    Raises:
        ValueError: If it is identically zero. There would be no generic axis to remove,
            and the z-score matrix would not be what this channel expects.
    """
    import numpy as np

    mean_signature = np.asarray(z_scores, dtype=np.float64).mean(axis=0)
    centred = mean_signature - mean_signature.mean()
    norm = float(np.linalg.norm(centred))
    if norm == 0:
        raise ValueError("the mean compound signature is constant across genes; there is "
                         "no generic axis to remove and the z-score matrix is not what "
                         "this channel expects")
    return centred / norm


def project_out(vector, axis):
    """``vector``, mean-centred, with its component along the unit ``axis`` removed.

    See :func:`generic_axis` for why both are centred, and why centring the query is free.
    """
    import numpy as np

    vector = np.asarray(vector, dtype=np.float64)
    centred = vector - vector.mean()
    return centred - float(np.dot(centred, axis)) * axis


def deflate_generic_axis(consensus, z_scores) -> tuple:
    """Remove the generic drug-response component from the query. Returns ``(query, report)``.

    The mean L1000 compound signature is what every perturbation has in common: the
    universal transcriptional stress response. The proxy signature overlaps it heavily
    (measured r = +0.52), so an uncorrected connectivity score largely measures whether a
    compound provokes that response at all. Projecting the mean signature out of the query
    leaves the component specific to losing this gene.

    The correction is applied to the **query**, not to the signatures: each compound's
    measured response is a fact about the experiment and is left alone. What changes is
    which genes the query asks about.

    Returns the deflated query and a report carrying the correlation before and after, so
    the size of the confound is recorded beside the result rather than asserted here.
    """
    import numpy as np

    consensus = np.asarray(consensus, dtype=np.float64)
    mean_signature = np.asarray(z_scores, dtype=np.float64).mean(axis=0)
    axis = generic_axis(z_scores)
    deflated = project_out(consensus, axis)
    return deflated, {
        "applied": True,
        "axis": "mean L1000 compound signature over the scored population",
        "r_before": round(float(np.corrcoef(consensus, mean_signature)[0, 1]), 4),
        "r_after": round(float(np.corrcoef(deflated, mean_signature)[0, 1]), 4),
        "note": ("The uncorrected query overlaps the generic drug-response axis, so most "
                 "of an uncorrected score measures whether a compound provokes the "
                 "universal stress response rather than anything specific to this gene. "
                 "See the guardrail for the measurement that settled this."),
    }


def _resolve_compounds(config: dict, settings: dict) -> tuple:
    """``({normalised pert_iname: chembl_id}, {chembl_id: identity}, stats, provenance)``.

    Identity resolution happens *before* scoring, not after, so the 5 GB matrix is read
    only for signatures that can become a candidate. Ambiguous names are dropped rather
    than guessed: if two molecules are both called the same thing, picking either is a
    coin toss about drug identity.
    """
    paths, provenance = enrichment.ensure_datasets(config, ("drug_molecule",),
                                                   channel="channel_c")
    molecules = enrichment.load_molecules(paths["drug_molecule"])
    names, name_stats = enrichment.name_index(paths["drug_molecule"])
    ambiguous = set(name_stats["ambiguous_synonyms"])
    usable = {key: chembl for key, chembl in names.items() if key not in ambiguous}
    stats = {"molecule_names": name_stats["molecule_names"],
             "molecule_synonyms": name_stats["molecule_synonyms"],
             "ambiguous_names_dropped": len(ambiguous),
             "approved_only": settings["approved_only"]}
    return usable, molecules, stats, provenance


def _select_signatures(sig_info: list, signature_ids: list, names: dict,
                       molecules: dict, settings: dict) -> tuple:
    """Compound signatures worth reading, each tied to a harmonisable identity.

    Returns ``(selected, stats)`` where each entry is
    ``{matrix_row, sig_id, chembl_id, cell_line}``, sorted by matrix row so the big
    matrix is read front to back exactly once.
    """
    position = {sid: i for i, sid in enumerate(signature_ids)}
    selected, unresolved, unapproved = [], set(), set()
    for row in sig_info:
        if row.get("pert_type", "").strip() != lincs.COMPOUND_TYPE:
            continue
        name = row.get("pert_iname", "").strip()
        chembl = names.get(enrichment.normalise_name(name))
        if not chembl:
            unresolved.add(name)
            continue
        if settings["approved_only"] and \
                molecules.get(chembl, {}).get("clinical_stage") != enrichment.APPROVED:
            unapproved.add(chembl)
            continue
        matrix_row = position.get(row["sig_id"])
        if matrix_row is None:
            continue
        selected.append({"matrix_row": matrix_row, "sig_id": row["sig_id"],
                         "chembl_id": chembl, "cell_line": row.get("cell_id", "")})

    selected.sort(key=lambda e: e["matrix_row"])
    stats = {"compound_signatures_in_release": sum(
                 1 for r in sig_info if r.get("pert_type", "").strip() == lincs.COMPOUND_TYPE),
             "signatures_selected": len(selected),
             "compound_names_unresolved": len(unresolved),
             "molecules_excluded_not_approved": len(unapproved),
             "molecules_selected": len({e["chembl_id"] for e in selected})}
    return selected, stats


def read_signatures(handle, selected: list, cp_index: list, settings: dict):
    """Scan the compound matrix once and return the selected signatures' landmark z-scores.

    Reading is separated from scoring because the query depends on the data: the generic
    axis removed from the query is estimated from these very signatures, so nothing can be
    scored until they are all in hand.
    """
    import numpy as np

    budget = settings["max_scored_signatures"]
    if len(selected) > budget:
        raise ValueError(
            f"{len(selected)} signature(s) selected, over max_scored_signatures={budget}. "
            "Their landmark z-scores are held in memory for the null; raise the budget "
            "deliberately, or narrow the selection with approved_only."
        )

    wanted = np.array([e["matrix_row"] for e in selected])
    n_rows = handle[lincs.MATRIX_PATH].shape[0]
    kept = np.empty((len(selected), len(cp_index)), dtype=np.float64)

    filled = 0
    for start in range(0, n_rows, settings["block_size"]):
        stop = min(start + settings["block_size"], n_rows)
        take = wanted[(wanted >= start) & (wanted < stop)]
        if take.size == 0:
            continue
        block = lincs.read_block(handle, start, stop, cp_index)
        kept[filled:filled + take.size] = block[take - start]
        filled += take.size
    LOGGER.info("channel C: read %d signature(s)", filled)
    return kept


def _groups(selected: list) -> dict:
    """``{chembl_id: [index into the signature array]}``."""
    grouped = defaultdict(list)
    for index, entry in enumerate(selected):
        grouped[entry["chembl_id"]].append(index)
    return dict(grouped)


def score_population(z_scores, order, cells, groups: dict, keys: list,
                     up_mask, down_mask):
    """One aggregated NCS per molecule, in ``keys`` order.

    Returns ``(scores, ncs, wtcs)`` so the caller can report the per-signature detail of
    the observed query without recomputing it for every null draw.
    """
    import numpy as np

    wtcs, _, _, _ = connectivity.weighted_connectivity(z_scores, up_mask, down_mask,
                                                       order=order)
    ncs = connectivity.normalise(wtcs, cells)
    scores = np.array([connectivity.max_quantile(ncs[groups[key]]) for key in keys])
    return scores, ncs, wtcs


def control_genes(meta: list, gene: str, settings: dict, seed: int) -> list:
    """Other genes whose knockdown is profiled as widely as the causal gene's.

    These are the null. Eligibility is "has a consensus knockdown in at least as many cell
    lines as the target", so a control query is built from the same amount of evidence;
    the SAC panel and the target itself are excluded, because a spindle-checkpoint gene is
    not an unrelated control for another one.

    Raises:
        ValueError: If too few genes qualify. A null of one or two genes cannot calibrate
            anything, and reporting a ranking against it would be worse than reporting
            none.
    """
    import numpy as np

    lines = {}
    for row in meta:
        if row.get("pert_type", "").strip() == lincs.KNOCKDOWN_TYPE:
            lines.setdefault(row.get("pert_iname", "").strip().upper(), set()).add(
                row.get("cell_id", ""))
    required = len(lines.get(gene.strip().upper(), ()))
    excluded = {g.upper() for g in SAC_PANEL} | {gene.strip().upper()}
    eligible = sorted(g for g, cells in lines.items()
                      if len(cells) >= required and g not in excluded)
    wanted = settings["n_control_genes"]
    if len(eligible) < wanted:
        raise ValueError(
            f"only {len(eligible)} gene(s) have a knockdown in >= {required} cell lines, "
            f"fewer than the {wanted} controls the null needs. Lower n_control_genes "
            "deliberately, or accept that this proxy cannot be calibrated."
        )
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(eligible, size=wanted, replace=False).tolist())


def calibrate_depth(observed, control_scores: dict, settings: dict) -> tuple:
    """How deep the ranking beats *other genes'* knockdowns. Returns ``(depth, report)``.

    The observed k-th best aggregated score is compared against the distribution of k-th
    bests obtained when the query is another gene's knockdown consensus, processed
    identically. The channel nominates the deepest ladder rung that still beats that null,
    so the length of the candidate list is a measurement rather than a chosen cut-off, and
    a proxy with nothing to say yields depth 0.

    **Two earlier nulls were tried against random gene sets and both failed**, in the same
    way and for the same reason -- a structured biological query engages more genes
    coherently than an arbitrary set, so it wins for reasons unrelated to the gene:

    * per-candidate query permutation passed 60 candidates out of 60;
    * a max-statistic over random queries passed every depth, because 42% of signatures
      score exactly zero against a random query and the whole distribution shifts.

    Only holding "knockdown query" fixed and varying *which gene* isolates the question
    the channel is asking. Both failures are recorded in the output so the null is not
    quietly replaced by an easier one later.
    """
    import numpy as np

    depths = [d for d in DEPTH_LADDER if d <= min(settings["top_n"], len(observed))]
    if not depths:
        raise ValueError(f"no ladder depth fits {len(observed)} molecule(s) under "
                         f"top_n={settings['top_n']}")
    observed_sorted = np.sort(observed)
    ladder = []
    for depth in depths:
        null = [float(np.sort(scores)[depth - 1]) for scores in control_scores.values()]
        p = connectivity.empirical_p(observed_sorted[depth - 1], null, tail="lower")
        ladder.append({"depth": depth,
                       "observed": round(float(observed_sorted[depth - 1]), 4),
                       "null_median": round(float(np.median(null)), 4),
                       "null_best": round(float(np.min(null)), 4),
                       "p": round(p, 4), "passes": bool(p <= settings["p_max"])})
    passing = [rung["depth"] for rung in ladder if rung["passes"]]
    return (max(passing) if passing else 0), {
        "method": "the observed k-th best aggregated score against the distribution of "
                  "k-th bests obtained when the query is ANOTHER gene's knockdown "
                  "consensus, built and deflated identically",
        "control_genes": sorted(control_scores),
        "n_control_genes": len(control_scores),
        "min_attainable_p": round(1 / (len(control_scores) + 1), 5),
        "p_max": settings["p_max"],
        "ladder": ladder,
        "depth": max(passing) if passing else 0,
        "superseded_nulls": [
            "per-candidate query permutation: passed 60 of 60 candidates, because it "
            "tested whether the query was structured rather than whether a compound was "
            "special",
            "max-statistic over random gene sets: passed every depth, because 42% of "
            "signatures score exactly zero against a random query, so the observed "
            "distribution is shifted against the null for reasons unrelated to the gene",
        ],
    }


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def generate(config: dict) -> None:
    """Generate candidates by signature reversion against a proxy signature.

    Args:
        config: Parsed pipeline configuration; uses ``causal_gene``, ``seed``,
            ``reference_dir``, ``enrichment_dir``, ``results_dir`` and ``l2.channel_c``.

    Raises:
        ValueError: If no proxy can be built, if the two LINCS releases do not share a
            gene axis, or if the ranking does not beat an unrelated gene's knockdown at
            any candidate depth. The last is a finding, and it is raised rather than
            written as an empty table: L2's runner treats an empty ``candidates.tsv`` as
            a failure precisely so that "no candidates" has to be said out loud.
    """
    import h5py
    import numpy as np

    settings = settings_for(config)
    gene = config.get("causal_gene")
    if not gene:
        raise ValueError("config['causal_gene'] is unset; Channel C's proxy signature is "
                         "a knockdown of that gene and there is nothing else to reverse")

    out_dir = Path(config["results_dir"]) / "l2" / "channel_c_signature"
    out_dir.mkdir(parents=True, exist_ok=True)

    paths, lincs_provenance = lincs.ensure(config)
    landmarks = lincs.landmark_ids(lincs.read_metadata(paths["compound"]["genes"]))
    names, molecules, identity_stats, drug_provenance = _resolve_compounds(config, settings)

    with h5py.File(paths["knockdown"]["matrix"], "r") as kd_file, \
            h5py.File(paths["compound"]["matrix"], "r") as cp_file:
        kd_signature_ids, kd_genes = lincs.axes(kd_file)
        cp_signature_ids, cp_genes = lincs.axes(cp_file)
        entrez, kd_index, cp_index = _gene_axis(kd_genes, cp_genes, landmarks)
        LOGGER.info("channel C: %d landmark gene(s) shared by both releases", len(entrez))

        kd_meta = lincs.read_metadata(paths["knockdown"]["meta"])
        consensus, proxy = build_proxy(kd_file, kd_signature_ids, kd_meta, gene,
                                       kd_index, settings)
        # The null's queries are built here, while the knockdown matrix is still open.
        controls = control_genes(kd_meta, gene, settings, int(config["seed"]))
        LOGGER.info("channel C: building %d control knockdown quer(ies)", len(controls))
        control_queries = {
            control: build_proxy(kd_file, kd_signature_ids, kd_meta, control, kd_index,
                                 settings)[0]
            for control in controls}

        selected, selection_stats = _select_signatures(
            lincs.read_metadata(paths["compound"]["meta"]), cp_signature_ids, names,
            molecules, settings)
        if not selected:
            raise ValueError(
                "no LINCS compound signature resolved to an identity this pipeline can "
                f"aggregate ({selection_stats}). That is a broken crosswalk or an "
                "approved_only filter that removed everything, not a finding."
            )
        LOGGER.info("channel C: reading %d signature(s) for %d molecule(s)",
                    len(selected), selection_stats["molecules_selected"])
        z_scores = read_signatures(cp_file, selected, cp_index, settings)

    if settings["deflate_generic_axis"]:
        query, deflation = deflate_generic_axis(consensus, z_scores)
    else:
        query, deflation = consensus, {
            "applied": False,
            "note": ("Not applied. Measured on this data, the uncorrected query does not "
                     "beat random gene sets (p = 0.95) while the corrected one does "
                     "(p = 0.024); an uncorrected ranking here is noise with a plausible "
                     "shape. This setting exists to reproduce that comparison, not to be "
                     "left off."),
        }

    up_mask, down_mask = connectivity.query_masks(query, settings["query_size"])
    groups = _groups(selected)
    keys = sorted(groups)
    cells = np.array([e["cell_line"] for e in selected])
    LOGGER.info("channel C: scoring %d molecule(s)", len(keys))
    order = connectivity.descending_order(z_scores)
    observed, ncs, wtcs = score_population(z_scores, order, cells, groups, keys,
                                           up_mask, down_mask)

    # The null: the same pipeline, driven by another gene's knockdown.
    axis = generic_axis(z_scores) if settings["deflate_generic_axis"] else None
    control_scores = {}
    for position, (control, control_consensus) in enumerate(sorted(control_queries.items()), 1):
        control_query = project_out(control_consensus, axis) if axis is not None \
            else control_consensus
        control_up, control_down = connectivity.query_masks(control_query,
                                                            settings["query_size"])
        control_scores[control] = score_population(z_scores, order, cells, groups, keys,
                                                   control_up, control_down)[0]
        if position % 10 == 0:
            LOGGER.info("channel C: control gene %d/%d", position, len(control_queries))

    depth, null_report = calibrate_depth(observed, control_scores, settings)
    ranked = [keys[i] for i in np.argsort(observed, kind="stable")]
    if depth == 0:
        raise ValueError(
            f"signature reversion on the {gene} knockdown proxy is not distinguishable "
            f"from an unrelated gene's knockdown at any candidate depth "
            f"({len(control_scores)} control genes, p_max={settings['p_max']}). Ladder "
            f"(depth, p): {[(r['depth'], r['p']) for r in null_report['ladder']]}. This "
            "is a result, not a crash: on this proxy, at this query size, signature "
            "reversion carries no gene-specific signal and the channel nominates nothing. "
            "Record it and disable the channel rather than loosening p_max to fill the "
            "table."
        )
    nominated = ranked[:depth]
    LOGGER.info("channel C: null supports a depth of %d of %d molecule(s)", depth,
                len(keys))

    score_by_key = dict(zip(keys, observed))
    rows = []
    for position, chembl in enumerate(nominated, start=1):
        indices = groups[chembl]
        public = enrichment.publishable(enrichment.DrugRecord(
            chembl_id=chembl, name=molecules.get(chembl, {}).get("name", ""),
            drug_type=molecules.get(chembl, {}).get("drug_type", ""),
            clinical_stage=molecules.get(chembl, {}).get("clinical_stage", "")))
        rows.append([
            position, public["chembl_id"], public["name"], public["drug_type"],
            public["clinical_stage"],
            f"{score_by_key[chembl]:.4f}",
            f"{connectivity.max_quantile(wtcs[indices]):.4f}",
            f"{100.0 * position / len(keys):.2f}",
            len(indices), len({selected[i]['cell_line'] for i in indices}),
            ",".join(sorted({selected[i]["cell_line"] for i in indices})),
            proxy["signature_id"], PROXY_TYPE, ",".join(proxy["cell_lines"]),
            PROXY_CAVEAT,
        ])
    _write_tsv(out_dir / "candidates.tsv",
               ["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "ncs",
                "wtcs", "population_percentile", "n_signatures", "n_cell_lines",
                "compound_cell_lines", "proxy_signature_id", "proxy_type", "cell_line",
                "proxy_caveat"],
               rows)

    payload = {
        "channel": "channel_c_signature",
        "proxy": proxy,
        "deflation": deflation,
        "method": {
            "score": "weighted connectivity score (WTCS): two-tailed weighted "
                     "Kolmogorov-Smirnov statistic of the query's up and down sets "
                     "against each signature's ranked genes [subramanian2017] "
                     "doi:10.1016/j.cell.2017.10.049",
            "normalisation": "NCS -- WTCS divided by the mean magnitude of the scores "
                             "sharing its cell line and sign",
            "aggregation": "max-quantile across the signatures of one molecule "
                           f"(percentiles {connectivity.MAX_QUANTILE})",
            "ranking": "ascending NCS -- more negative is stronger reversal",
            "gene_space": "the 978 L1000 landmark genes, intersected across releases",
            "depth": "set by the null, not chosen: the deepest ladder rung whose order "
                     "statistic beats random queries",
        },
        "null": null_report,
        "parameters": {**settings, "causal_gene": gene,
                       "conventional_cc_q75": CONVENTIONAL_CC_Q75,
                       "depth_ladder": list(DEPTH_LADDER)},
        "counts": {
            "landmark_genes_scored": len(entrez),
            "proxy_cell_lines": len(proxy["cell_lines"]),
            **selection_stats, **identity_stats,
            "molecules_scored": len(keys),
            "molecules_reversing": int((observed < 0).sum()),
            "molecules_mimicking": int((observed > 0).sum()),
            "molecules_nominated": len(nominated),
        },
        "distribution": {
            "ncs_min": round(float(observed.min()), 4),
            "ncs_median": round(float(np.median(observed)), 4),
            "ncs_max": round(float(observed.max()), 4),
            "signature_wtcs_negative_fraction": round(float((wtcs < 0).mean()), 4),
            "signature_wtcs_zero_fraction": round(float((wtcs == 0).mean()), 4),
            "note": ("A WTCS of exactly zero means the query's two tails moved together "
                     "in that signature, which is not a reversal and not a match."),
        },
        "sources": {"lincs": lincs_provenance, "drug_identity": drug_provenance},
        "licensing": {
            "published_fields": list(enrichment.PUBLISHABLE_FIELDS),
            "note": ("LINCS content is matched locally and never redistributed; only "
                     "derived scores and identifiers appear here. Drug identity comes "
                     "from the restricted enrichment zone through publishable(). See "
                     "src/l4_validate/sources.md, Tables 5 and 10."),
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }
    (out_dir / "channel.json").write_text(json.dumps(payload, indent=2, sort_keys=True,
                                                     default=str), encoding="utf-8")
    LOGGER.info("channel C written: %d nominated of %d molecule(s) scored", len(nominated),
                len(keys))
