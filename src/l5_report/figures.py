"""L5 -- figures.

Purpose
    Draw the figures that carry the scientific argument: what the safety triage refused,
    how the channels contributed, and how the pipeline scored on its own blinded
    benchmark.

Inputs
    Values already read from ``results/`` by :mod:`src.l5_report.run`.

Outputs
    PNG files under ``results/l5/figures/``, and the caption each one must be shown with.

Guardrail
    **No figure of this child's genome.** The per-chromosome aneuploidy burden is the
    obvious plot to draw and it is the one thing this module must not draw: in a
    ~50-patient population, which chromosomes are involved and at what mosaic fraction is
    close to an identifier. See :mod:`src.l5_report.publish_guard`.

    **A ranking is never drawn without its uncertainty.** The benchmark figure shows the
    bootstrap interval, not just the point estimate -- a bar chart of AUROC computed on
    eight positives would imply a precision that is not there.

    Rendering runs headless through the Agg backend, so the pipeline never tries to open a
    window on a machine with no display.
"""

from __future__ import annotations

import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

#: Print-safe, colour-blind-safe, and legible when a slide is projected badly.
INK = "#1a1a1a"
MUTED = "#5b5b5b"
RULE = "#d8d8d8"
STOP = "#8a1c1c"
WARN = "#c1761f"
COOL = "#2f6690"
OK = "#4a7a46"

FIGURE_DPI = 160


def _pyplot():
    """Matplotlib with the Agg backend, configured once."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": FIGURE_DPI, "savefig.dpi": FIGURE_DPI,
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": RULE, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
        "axes.spines.right": False, "figure.facecolor": "white",
        "savefig.bbox": "tight", "savefig.facecolor": "white",
    })
    return plt


def _save(fig, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    _pyplot().close(fig)
    LOGGER.info("figure written: %s", path.name)
    return path.name


def exclusions_by_reason(counts: dict, path: Path) -> tuple:
    """Why 1,880 of 1,963 candidates were refused. Returns ``(filename, caption)``.

    The headline figure of the whole report (D18): what the pipeline declined to propose
    is its strongest and most reproducible result. ``insufficient_evidence`` is drawn in a
    different colour because it is a *coverage* fact, not a safety finding, and reading it
    as "these drugs are dangerous" would misrepresent the triage.
    """
    plt = _pyplot()
    ordered = sorted(counts.items(), key=lambda kv: kv[1])
    labels = [k.replace("_", " ") for k, _ in ordered]
    values = [v for _, v in ordered]
    colours = [MUTED if key == "insufficient_evidence" else STOP for key, _ in ordered]

    fig, ax = plt.subplots(figsize=(7.2, 0.32 * len(ordered) + 1.1))
    bars = ax.barh(labels, values, color=colours, height=0.72)
    ax.bar_label(bars, padding=3, fontsize=8, color=MUTED, fmt="%d")
    ax.set_xlabel("candidates excluded")
    ax.set_xlim(0, max(values) * 1.16)
    ax.set_title("What the safety triage refused to propose", loc="left", weight="bold")
    ax.tick_params(length=0)
    ax.grid(axis="x", color=RULE, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    return _save(fig, path), (
        "Every candidate the pipeline declined, by the reason that fired. Grey is "
        "insufficient_evidence — a coverage fact, not a safety finding: openFDA "
        "describes drugs marketed in the United States and the channels nominate from a "
        "wider pool, so a candidate with no label is excluded because its status cannot "
        "be established. Red is a positive finding in the label text."
    )


def channel_contribution(per_channel: dict, by_convergence: dict, path: Path) -> tuple:
    """How many candidates each channel ranked, and how much they agreed.

    Drawn together on purpose. The left panel is what each channel contributed; the right
    is the uncomfortable half -- almost nothing is supported by two channels that each
    rank a narrow slice, which is the result the architecture hoped to avoid.
    """
    plt = _pyplot()
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.6, 3.0),
                                      gridspec_kw={"width_ratios": [1.15, 1]})

    names = sorted(per_channel, key=lambda k: per_channel[k], reverse=True)
    values = [per_channel[k] for k in names]
    bars = left.barh(names, values, color=COOL, height=0.6)
    left.bar_label(bars, padding=3, fontsize=8, color=MUTED, fmt="%d")
    left.set_xscale("log")
    left.set_xlabel("candidates ranked (log scale)")
    left.set_title("What each channel contributed", loc="left", weight="bold")
    left.tick_params(length=0)
    left.set_xlim(right=max(values) * 3)

    order = ["discriminating", "mixed", "broad_only", "none"]
    present = [(k, by_convergence.get(k, 0)) for k in order if k in by_convergence]
    colours = {"discriminating": OK, "mixed": WARN, "broad_only": MUTED, "none": RULE}
    bars = right.bar([k.replace("_", " ") for k, _ in present],
                     [v for _, v in present],
                     color=[colours[k] for k, _ in present], width=0.62)
    right.bar_label(bars, padding=3, fontsize=8, color=MUTED, fmt="%d")
    right.set_yscale("symlog")
    right.set_ylabel("candidates (symlog)")
    right.set_title("Cross-channel convergence", loc="left", weight="bold")
    right.tick_params(axis="x", length=0, labelrotation=20)
    for label in right.get_xticklabels():
        label.set_horizontalalignment("right")

    fig.tight_layout()
    return _save(fig, path), (
        "Left: how many candidates each channel ranked. Proximity ranks essentially every "
        "approved drug with a target in the interactome, which is why it is treated as a "
        "broad channel and its agreement alone is not reported as convergence. Right: how "
        "many candidates two or more channels agreed on. ZERO are discriminating: "
        "no candidate is supported by two channels that each rank a narrow slice. That is "
        "the weakest point in this submission and it is shown rather than described."
    )


def benchmark_recovery(rows, path: Path) -> tuple:
    """Blinded benchmark AUROC per scope, with the bootstrap interval.

    A channel whose positives are its own seed file is marked circular and drawn hatched:
    it cannot be scored on this benchmark, and leaving it off the figure entirely would
    hide that the pipeline knows it.
    """
    plt = _pyplot()
    scored = [r for r in rows if r.get("auroc") not in (None, "")]
    circular = [r for r in rows if r.get("circular")]

    fig, ax = plt.subplots(figsize=(6.6, 0.55 * len(scored) + 1.7))
    labels, centres, lows, highs = [], [], [], []
    for row in scored:
        labels.append(row["scope"])
        centres.append(float(row["auroc"]))
        lows.append(float(row["ci_low"]))
        highs.append(float(row["ci_high"]))

    positions = range(len(labels))
    for y, centre, low, high in zip(positions, centres, lows, highs):
        ax.plot([low, high], [y, y], color=COOL, linewidth=2.4, solid_capstyle="round")
        ax.plot([centre], [y], marker="D", color=COOL, markersize=7, zorder=3)
        ax.annotate(f"{centre:.2f}  [{low:.2f}, {high:.2f}]", (high, y),
                    textcoords="offset points", xytext=(8, 0), va="center",
                    fontsize=8, color=MUTED)

    ax.axvline(0.5, color=STOP, linewidth=1.1, linestyle="--")
    ax.annotate("chance", (0.5, -0.72), fontsize=8, color=STOP, ha="center")
    ax.set_yticks(list(positions), labels)
    ax.set_xlim(0.30, 1.30)
    ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("AUROC (bar = 95% bootstrap interval)")
    ax.set_ylim(-1.1, len(labels) - 0.4)
    ax.set_title("Blinded benchmark: does the pipeline recover known compounds?",
                 loc="left", weight="bold")
    ax.tick_params(length=0)
    ax.grid(axis="x", color=RULE, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    excluded = ", ".join(r["scope"] for r in circular)
    note = (f"Excluded as circular: {excluded}." if excluded else "")
    fig.tight_layout()
    return _save(fig, path), (
        "Recovery of a frozen positive set of published aneuploidy-stress compounds. "
        f"{note} The interval matters more than the point estimate: eight positives reach "
        "the ranking, so every estimate here is wide. Crossing 0.5 would mean the ranking "
        "is no better than chance. This measures whether the pipeline finds what the "
        "literature already contains — not whether it finds what would help this child."
    )


def evidence_tiers(tier_counts: dict, path: Path) -> tuple:
    """How the 83 survivors split by the *kind* of evidence behind them (D18)."""
    plt = _pyplot()
    labels = {"literature": "Tier 1\nliterature-supported",
              "network_only": "Tier 2\nnetwork proximity only"}
    keys = [k for k in ("literature", "network_only") if k in tier_counts]
    values = [tier_counts[k] for k in keys]

    fig, ax = plt.subplots(figsize=(5.6, 2.5))
    bars = ax.bar([labels.get(k, k) for k in keys], values,
                  color=[OK if k == "literature" else MUTED for k in keys], width=0.5)
    ax.bar_label(bars, padding=3, fontsize=10, color=INK, fmt="%d", weight="bold")
    ax.set_ylabel("surviving candidates")
    ax.set_ylim(0, max(values) * 1.22)
    ax.set_title("Survivors by kind of evidence, not by strength", loc="left",
                 weight="bold")
    ax.tick_params(length=0)
    ax.grid(axis="y", color=RULE, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    return _save(fig, path), (
        "Tier 1 rests on a verified, graded citation a reader can check. Tier 2 rests on a "
        "proximity score this pipeline computed, with nothing published linking the drug "
        "to this disease. Presenting them as one ranked list of 83 would invite a reader "
        "to treat those as the same kind of claim, so the report never does (decision "
        "D18)."
    )
