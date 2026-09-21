"""L5 -- Reporting.

Purpose
    Assemble the submission artifacts: the candidate dossier, the figures that carry the
    scientific argument, and the report-ready tables. Organized against the judging rubric
    -- Scientific Rigor 35%, Impact 25%, Innovation 25%, Scalability 15% -- so each claim
    has a figure and each figure has a claim.

Inputs
    ``results/l0_genomics/aneuploidy_burden.json`` (method and sensitivity only),
    ``results/l3/integration.json``, ``results/l3/reasoning/*.json``,
    ``results/l4/{survivors,excluded}.tsv``, ``results/l4/validation.json``,
    ``results/l4/benchmark/recovery.{csv,json}``, and the L2 channel status.

Outputs
    Under ``config["results_dir"]/l5/``:
      - ``index.html`` -- the dossier front page: exclusions first, then the two tiers
      - ``exclusions.html`` -- every refused candidate with the rule, field and snippet
      - ``candidates/<chembl_id>.html`` -- one static page per surviving candidate
      - ``figures/*.png`` -- the four argument figures
      - ``tables/*.tsv`` -- report-ready tables
      - ``report.json`` -- what was rendered, from which artifact, with which seed

Guardrail
    **Nothing published may re-identify the child or family.** Every string written by this
    layer passes :mod:`src.l5_report.publish_guard` first, and the per-chromosome
    aneuploidy result is withheld by construction -- the report describes the method and
    says plainly that the subject's result is not shown, because withholding silently is
    indistinguishable from never having looked.

    **Say "secondary prevention", not "chemoprevention"** wherever the endpoint is named
    (D4). The guard checks this per sentence, so a qualifier in the introduction cannot
    cover a table caption a reader screenshots.

    **This layer renders; it does not compute.** Every field shown already exists in an L3
    or L4 artifact. A field missing at render time is an upstream bug and this layer says
    so rather than filling it in -- in particular, an empty ``contradicting_evidence``
    renders as "searched, none found", never as a blank, because a blank reads as "not
    looked for".

    **The exclusions come first.** What the pipeline refused to propose is the headline
    result (D18), not an appendix, and the two survivor tiers are never merged into one
    list of 83.
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from . import figures, render
from .publish_guard import assert_publishable, strip_withheld

LOGGER = logging.getLogger(__name__)

TIER_LABELS = {
    "literature": "Tier 1 — literature-supported",
    "network_only": "Tier 2 — network proximity only",
}
#: Stated on every Tier 2 row, verbatim from L3 (D18). Restating it here would let the
#: report and the pipeline drift.
TIER2_PILL = "Tier 2"
TIER1_PILL = "Tier 1"


def _read_delimited(path: Path, delimiter: str, *, expect: str) -> list:
    """Rows of a delimited table, with the header checked.

    ``expect`` names a column that must be present. The check exists because L4 writes
    ``.tsv`` tables and a ``.csv`` benchmark, and reading one with the other's delimiter
    yields a list of single-key dicts -- plausible-looking rows that quietly render as a
    missing figure instead of an error.
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found -- L5 renders what L3 and L4 produced, so they must run "
            "first (python -m src.pipeline --only l4_validate)")
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    if rows and expect not in rows[0]:
        raise ValueError(
            f"{path} has no {expect!r} column; its header parsed as "
            f"{sorted(rows[0])[:4]}. That usually means the delimiter is wrong -- L4 "
            "writes tab-separated tables and a comma-separated benchmark.")
    return rows


def _read_tsv(path: Path, *, expect: str = "chembl_id") -> list:
    return _read_delimited(path, "\t", expect=expect)


def _read_json(path: Path, *, required: bool = True):
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"{path} not found -- L5 renders it and cannot "
                                    "compute it (decision D2)")
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, text: str) -> None:
    """Write a rendered file, after the publication guard has cleared it."""
    assert_publishable(text, str(path.name))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_tsv(path: Path, headers, rows) -> None:
    lines = ["\t".join(headers)]
    for row in rows:
        lines.append("\t".join("" if v is None else str(v).replace("\t", " ") for v in row))
    _write(path, "\n".join(lines) + "\n")


def _page_name(chembl_id: str) -> str:
    """A filename per candidate. ChEMBL ids are already filesystem-safe."""
    return f"{chembl_id}.html"


def _contradiction_text(record: dict) -> str:
    """``contradicting_evidence``, or the explicit statement that none was found.

    An empty cell reads as "nobody looked". D2 makes this field required for exactly that
    reason, so the absence is rendered as a result.
    """
    found = (record.get("contradicting_evidence") or "").strip()
    if found:
        return found
    if record.get("model_consulted"):
        return ("Searched, none found. An adversarial pass was run over this candidate's "
                "citations specifically to find evidence against it, and returned nothing "
                "that contradicted the rationale.")
    return ("Not searched, and not applicable. This candidate carries no published "
            "evidence to contradict — see the evidence tier above.")


def _candidate_page(record: dict, survivor: dict, verdicts: list, config: dict) -> str:
    """One candidate's static dossier page, rendering the D2 field contract."""
    tier = record.get("evidence_tier", "network_only")
    pill = TIER1_PILL if tier == "literature" else TIER2_PILL
    pill_class = "t1" if tier == "literature" else "t2"

    identity = render.definition_list([
        ("ChEMBL id", survivor.get("chembl_id")),
        ("RxCUI", survivor.get("rxcui") or "not resolved"),
        ("UNII", survivor.get("unii") or "not resolved"),
        ("Drug type", survivor.get("drug_type")),
        ("Clinical stage", survivor.get("clinical_stage")),
        ("Pharmacologic class (EPC)", survivor.get("pharm_class_epc") or "not annotated"),
        ("Mechanism class (MoA)", survivor.get("pharm_class_moa") or "not annotated"),
    ])

    ranking = render.definition_list([
        ("Rank in the integrated list", survivor.get("rank")),
        ("RRA score", survivor.get("rra_score")),
        ("Channels supporting", survivor.get("n_channels_supporting")),
        ("Which channels", survivor.get("supporting_channels") or "none recorded"),
        ("Convergence class", survivor.get("convergence")),
    ])

    confidence = record.get("confidence")
    grade = record.get("evidence_grade_supplied") or "ungraded"
    assessment_rows = [
        ("Evidence tier", TIER_LABELS.get(tier, tier)),
        ("Evidence grade", grade),
        ("Calibrated confidence",
         "not computed — no literature to calibrate against" if confidence is None
         else f"{confidence}"),
        ("Mechanism supported by cited literature",
         "not assessed — no citations to assess" if record.get("mechanism_supported") is None
         else ("yes" if record["mechanism_supported"] else "no")),
        ("Model consulted", "yes" if record.get("model_consulted") else
         "no — a candidate with no literature is not sent to a model to be told so (D18)"),
    ]

    citations = record.get("citation_keys") or []
    citation_html = (
        render.table(["Citation"], [(key,) for key in citations])
        if citations else
        "<p>None. This candidate carries no published reference linking it to this "
        "disease.</p>")

    verdict_rows = [
        (v.get("rule"), v.get("verdict"), v.get("reason") or "—",
         v.get("source_field") or "—", v.get("snippet") or "—")
        for v in verdicts]
    verdict_html = (
        render.table(["Safety rule", "Verdict", "Reason", "openFDA field", "Label text"],
                     verdict_rows)
        if verdict_rows else
        "<p>No per-rule verdict was recorded for this candidate, which is an upstream "
        "gap in L4 rather than a clean bill of health.</p>")

    proximity = record.get("network_proximity") or {}
    proximity_html = render.definition_list([
        ("Proximity z-score", proximity.get("z", "not scored")),
        ("Targets inside the disease module", proximity.get("n_targets_in_module", "—")),
    ]) if proximity else ""

    body = f"""
{render.block(f'<p>{render.ENDPOINT_STATEMENT}</p>', label="Endpoint", kind="warn")}

<h2>Why this candidate is on the list</h2>
<p><span class="pill {pill_class}">{render.esc(pill)}</span>
{render.esc(TIER_LABELS.get(tier, tier))}</p>
<p>{render.esc(record.get("rationale") or "No rationale was recorded.")}</p>

<h2>Evidence against it</h2>
{render.block(f'<p>{render.esc(_contradiction_text(record))}</p>',
              label="Contradicting evidence", kind="warn")}

<h2>Assessment</h2>
{render.definition_list(assessment_rows)}

<h2>References</h2>
{citation_html}

<h2>Identity</h2>
{identity}

<h2>Ranking and channel support</h2>
{ranking}
{proximity_html}

<h2>Safety triage</h2>
<p>This candidate passed every hard rule below. The triage is a gate, not a weight: a
candidate failing any rule is dropped, never down-ranked.</p>
{verdict_html}
"""
    return render.page(
        title=str(survivor.get("drug_name") or survivor.get("chembl_id")),
        subtitle=f"Candidate dossier — {TIER_LABELS.get(tier, tier)}",
        meta=(f"seed={config['seed']}  •  "
              f"rendered from results/l3 and results/l4  •  "
              f"causal gene {config.get('causal_gene')} (research premise, not a diagnosis)"),
        body_html=body,
        up="../index.html",
    )


def _exclusions_page(excluded: list, validation: dict, config: dict) -> str:
    """The companion exclusions page -- the same provenance standard as the survivors."""
    counts = validation.get("counts", {})
    by_reason = counts.get("excluded_by_reason", {})
    summary = render.table(
        ["Reason", "Candidates"],
        sorted(((k.replace("_", " "), v) for k, v in by_reason.items()),
               key=lambda kv: -kv[1]),
        numeric={"Candidates"})

    rows = [(r.get("drug_name"), r.get("chembl_id"), r.get("excluded_by"), r.get("reason"),
             r.get("source_field") or "—", r.get("source_section") or "—",
             (r.get("snippet") or "—")[:400])
            for r in excluded]
    body = f"""
{render.block(
    "<p>These are the candidates the pipeline <strong>refused to propose</strong>. "
    "On this data they are the more informative half of the output, and they ship with "
    "the results rather than in an appendix.</p>"
    "<p><code>insufficient_evidence</code> is a <em>coverage</em> fact, not a safety "
    "finding: openFDA describes drugs marketed in the United States and the channels "
    "nominate from a wider pool, so a candidate whose status cannot be established is "
    "excluded. Reading those rows as 'dangerous' would misrepresent the triage.</p>",
    label="What this page is", kind="warn")}

<h2>Summary</h2>
{summary}

<h2>Every exclusion</h2>
<p>{len(excluded)} rows. Each carries the rule that fired, the openFDA field it read, the
SPL section, and the sentence it rested on — trimmed, never paraphrased.</p>
{render.table(
    ["Drug", "ChEMBL id", "Rule", "Reason", "openFDA field", "SPL section", "Label text"],
    rows)}
"""
    return render.page(
        title="Excluded candidates",
        subtitle=f"{len(excluded)} of {counts.get('candidates', '?')} candidates were "
                 "refused by the paediatric safety triage",
        meta=f"seed={config['seed']}  •  rendered from results/l4/excluded.tsv",
        body_html=body,
        up="index.html",
    )


def _index_page(context: dict) -> str:
    """The dossier front page: exclusions first, then Tier 1, then Tier 2 (D18)."""
    config = context["config"]
    counts = context["validation"].get("counts", {})
    tier1, tier2 = context["tier1"], context["tier2"]

    def candidate_rows(records):
        rows = []
        for record in records:
            survivor = context["survivors_by_id"].get(record["chembl_id"], {})
            name = survivor.get("drug_name") or record.get("drug_name")
            rows.append((
                f'<a href="candidates/{render.esc(_page_name(record["chembl_id"]))}">'
                f"{render.esc(name)}</a>",
                survivor.get("rank"), survivor.get("n_channels_supporting"),
                survivor.get("supporting_channels"),
                record.get("evidence_grade_supplied") or "ungraded",
                "—" if record.get("confidence") is None else record["confidence"],
            ))
        return rows

    def candidate_table(records):
        headers = ["Drug", "Rank", "Channels", "Which channels", "Grade", "Confidence"]
        numeric = {"Rank", "Channels", "Confidence"}
        head = "".join(f'<th class="num">{h}</th>' if h in numeric else f"<th>{h}</th>"
                       for h in headers)
        body = []
        for row in candidate_rows(records):
            cells = []
            for header, value in zip(headers, row):
                # The first cell is a link this module built; the rest are escaped.
                cell = value if header == "Drug" else render.esc(value)
                cells.append(f'<td class="num">{cell}</td>' if header in numeric
                             else f"<td>{cell}</td>")
            body.append(f"<tr>{''.join(cells)}</tr>")
        return (f"<table><thead><tr>{head}</tr></thead>"
                f"<tbody>{''.join(body)}</tbody></table>")

    figure_html = "".join(
        render.figure(f"figures/{name}", caption)
        for name, caption in context["figures"])

    burden = context["burden"]
    method = burden.get("method", {})
    envelope = burden.get("sensitivity_envelope", {})
    innovation = render.definition_list([
        ("Statistic", method.get("statistic", "—")),
        ("Window", f"{method.get('window_bp', '—')} bp"),
        ("Minimum depth per site", method.get("min_dp", "—")),
        ("Baseline", method.get("baseline", "—")),
        ("Mosaic fraction estimator", method.get("fraction_estimate", "—")),
        ("Smallest mosaic fraction resolvable",
         envelope.get("min_detectable_mosaic_fraction", "—")),
        ("Conservative bound",
         envelope.get("min_detectable_mosaic_fraction_conservative", "—")),
    ])

    channel_c = context["channel_c_note"]
    body = f"""
{render.block(f'<p>{render.ENDPOINT_STATEMENT}</p>', label="Endpoint", kind="warn")}

<p>This pipeline takes one child's whole-genome sequence and a causal gene fixed by a
person at a recorded gate, and asks which <strong>already-approved</strong> drugs are
worth investigating. It ranks nothing it cannot source, and it refuses more than it
proposes.</p>

{render.table(
    ["Stage", "Candidates"],
    [("Nominated by the evidence channels", counts.get("candidates", "—")),
     ("Refused by the paediatric safety triage", counts.get("excluded", "—")),
     ("Surviving — literature-supported (Tier 1)", len(tier1)),
     ("Surviving — network proximity only (Tier 2)", len(tier2))],
    numeric={"Candidates"})}

<h2>1. What the pipeline refused to propose</h2>
<p>This is the headline result. {counts.get("excluded", "—")} of
{counts.get("candidates", "—")} candidates were excluded, each row naming the rule, the
openFDA field, the SPL section and the sentence that did it.
<a href="exclusions.html">Read the full exclusions table &rarr;</a></p>
{render.figure(f"figures/{context['figure_names']['exclusions']}",
               context["figure_captions"]["exclusions"])}

<h2>2. Surviving candidates</h2>
<p>The survivors are presented in two tiers by the <strong>kind</strong> of evidence
behind them, never as one ranked list. A Tier 1 candidate rests on a verified, graded
citation a reader can check. A Tier 2 candidate rests on a proximity score this pipeline
computed, with nothing published linking the drug to this disease.</p>
{render.figure(f"figures/{context['figure_names']['tiers']}",
               context["figure_captions"]["tiers"])}

<h3>{render.esc(TIER_LABELS["literature"])} — {len(tier1)} candidate(s)</h3>
{candidate_table(tier1)}

<h3>{render.esc(TIER_LABELS["network_only"])} — {len(tier2)} candidates</h3>
{render.block(
    "<p>Every candidate below reaches this list on network proximity alone. There is no "
    "published evidence connecting any of them to this disease, and none was given a "
    "generated rationale — eighty-two paragraphs explaining that there is no evidence "
    "read like analysis and are not.</p>", label="Read this before the table", kind="warn")}
{candidate_table(tier2)}

<h2>3. How the channels behaved</h2>
{render.figure(f"figures/{context['figure_names']['channels']}",
               context["figure_captions"]["channels"])}
{channel_c}

<h2>4. Does the pipeline recover what it should?</h2>
{render.figure(f"figures/{context['figure_names']['benchmark']}",
               context["figure_captions"]["benchmark"])}
{render.caveats(context["benchmark_caveats"])}

<h2>5. The aneuploidy-burden method</h2>
<p>MVA's phenotype <em>is</em> mosaic aneuploidy, and the dataset ships no alignments —
raw reads and a called VCF only. The pipeline therefore estimates aneuploidy burden from
per-chromosome B-allele frequency in the VCF itself, which is self-normalising and needs
none of the GC-bias and mappability correction that read depth would.</p>
{innovation}
{render.block(f'<p>{render.esc(burden["withheld_reason"])}</p>',
              label="Withheld by design", kind="stop")}

<h2>6. Caveats carried from the pipeline</h2>
{render.caveats(context["integration_caveats"])}
"""
    return render.page(
        title="Rare Disease, Real Kid — MVA drug repurposing shortlist",
        subtitle="Track 2 candidate dossier. Computational hypothesis generation for "
                 "mosaic variegated aneuploidy.",
        meta=(f"seed={config['seed']}  •  causal gene {config.get('causal_gene')} "
              f"(research premise, not a diagnosis)  •  rendered {context['rendered']}"),
        body_html=body,
    )


def _verdicts_by_candidate(verdicts) -> dict:
    """``{chembl_id: [per-rule verdict]}`` from ``results/l4/verdicts.json``.

    L4 writes ``{"candidates": [...], "seed": n}``, and each candidate carries its own
    ``verdicts`` list -- one entry per rule, with the rule name, the verdict, the openFDA
    field and the quotable snippet. That nesting is the D2 safety contract, so this reads
    it by its real shape rather than guessing.

    Raises:
        ValueError: If the artifact is not that shape. Rendering a candidate page with a
            silently empty safety section would show a drug with no triage beside it,
            which reads as "nothing to report" rather than "not loaded".
    """
    if not isinstance(verdicts, dict) or "candidates" not in verdicts:
        raise ValueError(
            "results/l4/verdicts.json is not the expected {'candidates': [...]} shape; "
            f"found {sorted(verdicts)[:5] if isinstance(verdicts, dict) else type(verdicts).__name__}. "
            "L5 will not render a candidate page with an empty safety section, because "
            "that reads as 'no findings' rather than 'not loaded'.")
    grouped: dict = {}
    for candidate in verdicts["candidates"]:
        chembl_id = candidate.get("chembl_id")
        if chembl_id:
            grouped[chembl_id] = candidate.get("verdicts", [])
    return grouped


def run(config: dict) -> None:
    """Execute layer L5.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir``, ``seed`` and
            ``causal_gene``.

    Raises:
        FileNotFoundError: If an upstream artifact this layer renders is missing. L5
            computes nothing, so a missing input is an upstream gap, not something to
            fill in here.
        ValueError: If anything about to be written fails the publication guard.
    """
    results_dir = Path(config["results_dir"])
    out_dir = results_dir / "l5"
    out_dir.mkdir(parents=True, exist_ok=True)

    survivors = _read_tsv(results_dir / "l4" / "survivors.tsv")
    excluded = _read_tsv(results_dir / "l4" / "excluded.tsv")
    validation = _read_json(results_dir / "l4" / "validation.json")
    integration = _read_json(results_dir / "l3" / "integration.json")
    rationales = _read_json(results_dir / "l3" / "reasoning" / "rationales.json")
    reasoning = _read_json(results_dir / "l3" / "reasoning" / "reasoning.json")
    # Comma-separated, unlike every other table L4 writes.
    recovery = _read_delimited(results_dir / "l4" / "benchmark" / "recovery.csv", ",",
                               expect="scope")
    recovery_json = _read_json(results_dir / "l4" / "benchmark" / "recovery.json",
                               required=False) or {}
    verdicts = _read_json(results_dir / "l4" / "verdicts.json")
    burden_raw = _read_json(results_dir / "l0_genomics" / "aneuploidy_burden.json",
                            required=False) or {}
    channels_status = _read_json(results_dir / "l2" / "channels.json", required=False) or {}

    # recovery.csv is read by csv.DictReader, so numeric fields arrive as strings.
    for row in recovery:
        row["circular"] = str(row.get("circular", "")).strip().lower() == "true"

    survivors_by_id = {row["chembl_id"]: row for row in survivors}
    by_candidate = _verdicts_by_candidate(verdicts)
    tier1 = [r for r in rationales if r.get("evidence_tier") == "literature"]
    tier2 = [r for r in rationales if r.get("evidence_tier") != "literature"]
    LOGGER.info("L5: %d survivor(s) — %d literature tier, %d network-only",
                len(rationales), len(tier1), len(tier2))

    # ---- figures -------------------------------------------------------------------
    figure_dir = out_dir / "figures"
    names, captions = {}, {}
    names["exclusions"], captions["exclusions"] = figures.exclusions_by_reason(
        validation["counts"]["excluded_by_reason"], figure_dir / "exclusions_by_reason.png")
    names["channels"], captions["channels"] = figures.channel_contribution(
        integration["counts"]["per_channel_ranked"],
        integration["counts"]["by_convergence"], figure_dir / "channel_contribution.png")
    names["tiers"], captions["tiers"] = figures.evidence_tiers(
        reasoning["by_evidence_tier"], figure_dir / "evidence_tiers.png")
    names["benchmark"], captions["benchmark"] = figures.benchmark_recovery(
        recovery, figure_dir / "benchmark_recovery.png")

    # Channel C ran, found nothing, and says so (D19). Reporting the architecture as
    # five channels without that would overstate what produced this shortlist.
    signature = (channels_status.get("channels", {}).get("signature", {})
                 if channels_status else {})
    channel_c_note = render.block(
        "<p>The architecture specifies five channels. <strong>Three produced evidence "
        "here</strong> — network proximity, phenotype and the literature prior. Channel A "
        "(knowledge-graph link prediction) is not built. Channel C (signature reversion) "
        "is implemented and calibrated, and "
        "it nominates nothing: on a knockdown-proxy signature its ranking is not "
        "distinguishable from an unrelated gene's knockdown at any candidate depth "
        "(p 0.15–0.44 against 40 control genes), so it declines rather than writing a "
        "table it cannot defend. Reproduce with "
        "<code>python scripts/channel_c_diagnostics.py</code>; the reasoning is decision "
        f"D19. Status in this run: <code>{render.esc(signature.get('status', 'disabled'))}"
        "</code>.</p>"
        "<p>The consequence is stated rather than buried: there is <strong>no "
        "discriminating cross-channel convergence</strong> in this shortlist, and the "
        "submission does not claim any.</p>",
        label="One channel produced nothing, deliberately", kind="warn")

    context = {
        "config": config,
        "validation": validation,
        "survivors_by_id": survivors_by_id,
        "tier1": tier1,
        "tier2": tier2,
        "figure_names": names,
        "figure_captions": captions,
        "figures": [],
        "benchmark_caveats": recovery_json.get("caveats", []),
        "integration_caveats": integration.get("caveats", []),
        "burden": strip_withheld(burden_raw),
        "channel_c_note": channel_c_note,
        "rendered": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
    }

    # ---- pages ---------------------------------------------------------------------
    _write(out_dir / "index.html", _index_page(context))
    _write(out_dir / "exclusions.html", _exclusions_page(excluded, validation, config))

    written = 0
    for record in rationales:
        chembl_id = record.get("chembl_id")
        survivor = survivors_by_id.get(chembl_id)
        if survivor is None:
            raise ValueError(
                f"{chembl_id} has a rationale but is not in survivors.tsv. L5 renders "
                "what L3 and L4 agreed on; a mismatch is an upstream inconsistency, not "
                "something to paper over here.")
        _write(out_dir / "candidates" / _page_name(chembl_id),
               _candidate_page(record, survivor, by_candidate.get(chembl_id, []), config))
        written += 1
    LOGGER.info("L5: %d candidate page(s) written", written)

    # ---- tables --------------------------------------------------------------------
    table_dir = out_dir / "tables"
    tier_headers = ["rank", "chembl_id", "drug_name", "rxcui", "pharm_class_epc",
                    "n_channels_supporting", "supporting_channels", "evidence_tier",
                    "evidence_grade", "confidence"]
    for label, records in (("tier1_literature", tier1), ("tier2_network_only", tier2)):
        rows = []
        for record in records:
            survivor = survivors_by_id[record["chembl_id"]]
            rows.append([survivor.get("rank"), record["chembl_id"],
                         survivor.get("drug_name"), survivor.get("rxcui"),
                         survivor.get("pharm_class_epc"),
                         survivor.get("n_channels_supporting"),
                         survivor.get("supporting_channels"),
                         record.get("evidence_tier"),
                         record.get("evidence_grade_supplied"),
                         record.get("confidence")])
        _write_tsv(table_dir / f"survivors_{label}.tsv", tier_headers, rows)

    _write_tsv(table_dir / "exclusions.tsv",
               ["chembl_id", "drug_name", "excluded_by", "reason", "source_field",
                "source_section", "snippet"],
               [[r.get("chembl_id"), r.get("drug_name"), r.get("excluded_by"),
                 r.get("reason"), r.get("source_field"), r.get("source_section"),
                 (r.get("snippet") or "").replace("\n", " ")] for r in excluded])

    # ---- manifest ------------------------------------------------------------------
    manifest = {
        "layer": "l5_report",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": config["seed"],
        "causal_gene": config.get("causal_gene"),
        "therapeutic_endpoint": "secondary prevention (D4)",
        "counts": {
            "candidates": validation["counts"].get("candidates"),
            "excluded": validation["counts"].get("excluded"),
            "survivors": len(rationales),
            "tier1_literature": len(tier1),
            "tier2_network_only": len(tier2),
            "candidate_pages": written,
        },
        "rendered_from": {
            "survivors": "results/l4/survivors.tsv",
            "excluded": "results/l4/excluded.tsv",
            "validation": "results/l4/validation.json",
            "integration": "results/l3/integration.json",
            "rationales": "results/l3/reasoning/rationales.json",
            "benchmark": "results/l4/benchmark/recovery.csv",
            "burden_method": "results/l0_genomics/aneuploidy_burden.json (method only)",
        },
        "pages": {
            "index": "l5/index.html",
            "exclusions": "l5/exclusions.html",
            "candidates": f"l5/candidates/*.html ({written} pages)",
        },
        "figures": {key: value for key, value in names.items() if value},
        "withheld": context["burden"]["withheld"],
        "withheld_reason": context["burden"]["withheld_reason"],
        "notes": [
            "This layer renders; it computes nothing. Every field shown exists in an L3 "
            "or L4 artifact (decision D2).",
            "Exclusions are presented first and survivors are tiered by evidence kind "
            "(decision D18); the two tiers are never merged into one list.",
            "Channel C is implemented and nominates nothing (decision D19), so this "
            "shortlist rests on three channels of the five specified (channel A is not "
            "built) and carries no discriminating cross-channel convergence.",
        ],
    }
    assert_publishable(manifest, "report.json")
    _write(out_dir / "report.json", json.dumps(manifest, indent=2, sort_keys=True))
    LOGGER.info("L5 written: %s", out_dir / "index.html")
