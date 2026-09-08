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

## D4 — Therapeutic endpoint is `chemoprevention`, meaning *secondary* prevention (2026-09-08)

Gate G2. `therapeutic_endpoint: chemoprevention` in
[config/pipeline.yaml](../config/pipeline.yaml).

**The framing is corrected, not just chosen.** The scaffold defined this endpoint as
"aneuploidy-buffering / remove pre-malignant clones" — primary prevention, stopping a first cancer.
The proband has already had one: the phenotype document records a malignancy as the oncological
event that triggered investigation. So the endpoint is **secondary** prevention — recurrence and
second-primary risk reduction in a child with a demonstrated cancer-predisposition phenotype. Every
report, figure, and candidate page must say so; "chemoprevention" unqualified overstates what this
is.

**Why not the alternatives.**

- `mitotic_fidelity` — no approved-drug route to restoring checkpoint function, and if L0 calls
  BUB1B or TRIP13 the checkpoint is already impaired, making SAC-directed inhibitors mechanistically
  counterproductive. As a ranked output it would be empty or unfounded. Carried as discussion.
- `symptomatic` — cheapest to run and Channel D is ready, but the documented manifestations are
  already actively managed clinically, so a computational shortlist adds least there, and it does
  not address what threatens this child. These candidates still surface through Channel D tagged
  `endpoint="symptomatic"`, so choosing chemoprevention does not suppress them.

**This decision does not presuppose L0.** The cancer risk here is established empirically by the
phenotype — the tumour occurred — not inferred from a genotype. So G2 could be settled before G1
without guessing the causal gene.

**Obliges.**

- L1 treats the **downstream** target set (aneuploidy-stress buffering) as primary; the upstream set
  is still emitted and discussed.
- L4's genotoxic gate binds hard: much of the mechanistically interesting literature is cytotoxic,
  and a cancer-predisposed child is the worst possible recipient of a genotoxin. The resulting
  exclusions table is a headline result, not a loss.
- L4's `bbb_penetration` criterion is **not** required. The documented phenotype does not establish a CNS requirement so the scaffold's assumption of CNS
  involvement does not hold. Record the criterion as not-applicable rather than scoring it.
- L5 states secondary prevention explicitly, and never uses efficacy language.
- The L4 benchmark's positive set stays coherent: "does the pipeline recover known
  aneuploidy/SAC-relevant compounds?" is a chemoprevention-shaped question. A different endpoint
  would have required rebuilding it — see the contamination warning in
  [docs/research/aneuploidy-selective-compounds.md](../docs/research/aneuploidy-selective-compounds.md).

---

## Open

- **`causal_gene`** (gate G1) — still `null`. A finding from L0, never a setting.
- **Hackathon close date** — unknown, so the 30-day deletion deadline cannot be computed. See
  [docs/data_custody.md](../docs/data_custody.md).
- **APA title casing** in [docs/references.md](../docs/references.md) — Crossref preserves publisher
  casing; a sentence-case pass is owed before submission.
