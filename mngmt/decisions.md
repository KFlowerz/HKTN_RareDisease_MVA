# Decisions

Choices that shape what gets built, with the reasoning that produced them. Recorded so a later
session — or a judge — can see why the pipeline looks the way it does, and so a decision is not
silently reversed by someone who never saw the argument against it.

Format: what was decided, when, why, and what it obliges. Superseded decisions stay, marked.

---

## D1 — The deliverable is a batch pipeline, not an application (2026-09-02)

`python -m src.pipeline` reads a config, runs L0→L5, and writes artifacts. No server, no persistent
state, no interactive UI.

**Why.** The rubric weights Scientific Rigor 35%, Impact 25%, Innovation 25%, Scalability 15%.
Nothing there rewards an interface. An interactive explorer would cost roughly a week of a five-week
solo build and earn none of it.

There is also a safety argument, and it is the stronger one: an interactive drug-filtering tool
*presents* as a clinical decision aid regardless of the disclaimer attached to it. This project is
explicitly hypothesis generation for a cancer-predisposed child, and
[COMPLIANCE.md](../COMPLIANCE.md) puts clinical recommendation out of scope. A static, dated,
citation-bearing document is read as an argument. A filterable app is read as advice.

**Obliges.** No web framework enters `environment.yml`. Reproducibility is carried by the pinned
environment, the seed, and the config — not by a hosted instance.

---

## D2 — Human-readable output is a static candidate dossier (2026-09-02)

L5 renders one page per surviving candidate, plus a companion exclusions table. Machine-readable
tables are still emitted; the dossier is the human-facing form.

**Why.** A CSV is a poor answer to Impact 25% and to the explainability requirement, and it gives the
3-minute video nothing to show but a terminal scroll. The dossier is also the honest shape for this
work: it forces every candidate to arrive with its contradicting evidence and its safety verdict
attached, rather than as a rank in a list.

**This is a W4 obligation, not a W5 one.** The dossier can only render fields that L3 and L4 emit.
Discovering a missing field in W5 means re-running the Claude reasoning step, which costs API spend
and breaks the checkpoint discipline in [src/l3_integrate/run.py](../src/l3_integrate/run.py).

### Per-candidate field contract

L3 and L4 must emit all of the following for every surviving candidate:

| Field | Source | Note |
|---|---|---|
| `rxcui`, `preferred_name` | L3 harmonize | RxNorm backbone; harmonize before aggregating |
| `pharm_class`, `moa_class` | L4 (RxClass / MED-RT) | CC-BY-clean core only |
| `rank_aggregate`, `score` | L3 aggregate | RRA/Borda |
| `n_channels_supporting` | L3 aggregate | Cross-channel convergence is the headline signal |
| `per_channel_rank` | L3 aggregate | One column per channel; absent ≠ zero |
| `rationale` | L3 Claude | May introduce no mechanism no channel cited |
| `contradicting_evidence` | L3 Claude | **Required, not optional.** Empty must mean "searched, found none" — recorded as such, never a blank cell |
| `confidence`, `confidence_basis` | L3 Claude | Calibrated, with what drove it |
| `evidence_grade` | L2 channel E / L3 | in-vitro / in-vivo / clinical — travels with the claim |
| `references[]` | L2 channel E | PMID or DOI per claim, each resolved against a real record |
| `safety_verdict[]` | L4 triage | One row per rule: rule name, verdict, source field, quotable label snippet |
| `provenance` | all | Artifact paths under `results/`, `seed`, config hash — so the row can be regenerated |
| `caveats[]` | L5 | Hypothesis-only disclaimer; Channel C proxy caveat where that channel contributed |

The exclusions table carries `rxcui`, `preferred_name`, the rule that fired, the source field, and
the snippet — the same provenance standard as the survivors. What the pipeline refused to propose is
a headline result, not an appendix.

**Obliges.** L5 renders; it does not compute. Any field the dossier shows must already exist in an
L3 or L4 artifact, so the dossier can be regenerated without re-running anything upstream.

---

## D3 — Audience is judges first, clinicians as a stated aspiration (2026-09-02)

The report argues the method to judges. The dossier demonstrates what a clinician-facing output
*would* look like, without claiming to be one.

**Why.** [CLAUDE.md](../CLAUDE.md) names treating clinicians as operators, but this is n=1 in a
population of roughly 50 patients, and no output here is a treatment recommendation. Writing for
clinicians as the primary audience would imply a readiness the evidence does not support. Writing
for judges while keeping the output legible is honest about both.

**Obliges.** Every candidate page states plainly that it is a computational hypothesis. Rationales
are written to be readable, not to be actionable. No dosing, no clinical guidance, no efficacy
language anywhere in the dossier, the report, or the video.

---

## Open

- **`therapeutic_endpoint`** (gate G2) — `chemoprevention` | `symptomatic` | `mitotic_fidelity`.
  Still `null`. This is a product decision as much as a scientific one: it determines what the ranked
  list *means*. Leaving it open costs roughly double the L4 annotation work, since both target sets
  are carried forward.
- **Hackathon close date** — unknown, so the 30-day deletion deadline cannot be computed. See
  [docs/data_custody.md](../docs/data_custody.md).
