# MVA Track 2 — Drug Repurposing

Scaffold for **Track 2 (Drug Repurposing)** of *Rare Disease, Real Kid: MVA Hackathon 2026*.

**Disease.** Mosaic variegated aneuploidy (MVA, OMIM 257300) — autosomal-recessive, biallelic
loss-of-function in mitotic spindle-assembly-checkpoint (SAC) genes (commonly **BUB1B**, **CEP57**,
**TRIP13**; also *BUB1*, *BUB3*, *CEP192*). Chromosome missegregation produces mosaic aneuploidy, with
growth restriction, microcephaly, developmental delay/seizures, and cancer predisposition. **No
disease-modifying therapy exists**; management is symptomatic plus cancer surveillance.

**What this repo produces.** Given the causal gene/pathway and phenotype, a **ranked, evidence-backed,
safety-filtered** shortlist of *approved* drugs worth investigating — each with a mechanistic rationale
and full provenance.

---

## ⚠️ Hypothesis generation only — not a clinical recommendation

Every output of this pipeline is a **computational hypothesis**. This is n=1 in a population of roughly 50
patients worldwide. Nothing here is a treatment recommendation, a claim of efficacy, or a substitute for
clinical judgement.

**Explicitly out of scope:**

- Wet-lab validation
- De-novo molecule design
- Any clinical recommendation
- Any recontact with the subject, family, or MVA Society contacts
- Any efficacy claim

See [COMPLIANCE.md](COMPLIANCE.md) for the full checklist and [DATA.md](DATA.md) for data handling.

---

## Pipeline overview (L0 → L5)

The design is **multi-channel**: five weak, independent evidence generators run in parallel and are
combined by consensus. That is deliberate — any single repurposing method is fragile on a hyper-rare
disease, so convergence across methods carries the signal rather than any one score.

| Layer | Package | Responsibility |
|---|---|---|
| **L0** | [src/l0_genomics/](src/l0_genomics/) | Ingest the single-sample WGS VCF, annotate variants, focus on biallelic SAC-gene variants, reconcile with the Track-1 validated causal variant. *Innovation hook:* per-chromosome **B-allele frequency** from `FORMAT/AD` quantifies **aneuploidy burden** and estimates mosaic fraction. The dataset ships no BAM — see [DATA.md](DATA.md). |
| **L1** | [src/l1_target/](src/l1_target/) | Causal gene → protein/complex → interactome **disease module**. Defines an **upstream** target set (restore mitotic fidelity) and a **downstream** set (buffer aneuploidy stress / chemoprevention). |
| **L2** | [src/l2_channels/](src/l2_channels/) | Five parallel candidate generators — see below. |
| **L3** | [src/l3_integrate/](src/l3_integrate/) | Harmonize drug identities on the **RxNorm RxCUI** backbone; rank-aggregate (RRA/Borda) preferring cross-channel convergence; a **Claude-in-the-loop** step that synthesizes each rationale, actively searches for contradicting evidence, and emits a calibrated confidence. |
| **L4** | [src/l4_validate/](src/l4_validate/) | In-silico validation + **pediatric safety triage**, including **hard exclusion of genotoxic / cancer-risk-increasing agents** (MVA is cancer-predisposing). Plus a **blinded internal benchmark** — there is no external ground truth for Track 2. |
| **L5** | [src/l5_report/](src/l5_report/) | Renders the **candidate dossier** — one static page per surviving candidate, a companion exclusions page, the rubric-aligned figures and report-ready tables. Renders only; it computes nothing (D2), and every string passes a publication guard before it reaches a file (D20). |

### L2 channels

| Channel | Module | Method | Why it is here |
|---|---|---|---|
| **A** | `channel_a_kg.py` | Knowledge-graph link prediction (Open Targets / PrimeKG / Hetionet + GNN, explainable meta-paths) | Anchored at the **gene** level, not the disease level — MVA is absent or near-empty in most KGs (cold-start defense). |
| **B** | `channel_b_proximity.py` | Network proximity (Guney/Barabási) between drug-target sets and the disease module | Does **not** require the disease to exist in a KG at all. |
| **C** | `channel_c_signature.py` | Signature reversion against a **PROXY** signature | There is **no patient RNA-seq**, so the query is a LINCS L1000 knockdown of the causal gene. **Built, calibrated, and it nominates nothing** — its ranking is not distinguishable from an unrelated gene's knockdown (p 0.15–0.44 against 40 control genes), so it declines rather than writing a table it cannot defend. Disabled in the shipped config; see [D19](mngmt/decisions.md) and `scripts/channel_c_diagnostics.py`. |
| **D** | `channel_d_phenotype.py` | Phenotype/HPO-driven (Monarch / Orphanet / Open Targets) | Reaches symptomatic candidates that mechanism-first channels miss. |
| **E** | `channel_e_prior.py` | Aneuploidy-stress literature prior + agentic literature mining (RareAgent-style) | Captures published aneuploidy-tolerance biology no structured resource encodes. |

[src/pipeline.py](src/pipeline.py) is the orchestrator: it loads `config/pipeline.yaml`, seeds the RNG, and
runs L0 → L5 in order, writing to `results/` (gitignored).

---

## Status

**All six layers are built.** L0 → L5 run end to end and produce the dossier.

Four of the five L2 channels produce evidence. Channel A is unbuilt. **Channel C is built and
produces nothing**, deliberately: it is calibrated against other genes' knockdowns, fails that test,
and refuses to nominate ([D19](mngmt/decisions.md)). The consequence is stated rather than buried —
this shortlist carries **no discriminating cross-channel convergence**, and the submission does not
claim any.

What the pipeline produced on this dataset:

| | |
|---|---|
| Candidates nominated by the channels | 1,963 |
| **Refused by the paediatric safety triage** | **1,880** |
| Surviving, literature-supported (Tier 1) | 1 |
| Surviving, network proximity only (Tier 2) | 82 |

The exclusions are the headline result, not an appendix ([D18](mngmt/decisions.md)): each refused
row names the rule, the openFDA field, the SPL section and the sentence that did it.

### Reading the report

```bash
python -m src.pipeline --only l5_report
# then open results/l5/index.html
```

`results/` is gitignored, so the dossier is built rather than committed. It is static HTML with no
JavaScript and no template engine — a deliberate choice ([D2](mngmt/decisions.md)): an interactive
drug filter presents as a clinical decision aid whatever disclaimer is attached to it.

**The per-chromosome aneuploidy result is not published**, in any form, including aggregate counts
([D20](mngmt/decisions.md)). The method and the smallest mosaic fraction it can resolve are; the
child's result from it is not. In a ~50-patient population that profile is close to an identifier.

## Quickstart

```bash
conda env create -f environment.yml
conda activate mva-track2

# L0's local variant annotator database (public reference data, not patient data).
# -noLog stops snpEff reporting usage statistics to its server -- pass it on every call.
snpEff download -noLog GRCh38.115

# data/ and results/ are gitignored, so they do not survive a clone — create them:
mkdir -p data results

# download the gated dataset into ./data — see DATA.md
python -m src.pipeline            # currently raises NotImplementedError at L0

pytest -q                         # smoke test: every layer imports, run() is callable
```

## Configuration

[config/pipeline.yaml](config/pipeline.yaml) drives everything. Two fields are gate decisions, set by
a person from recorded evidence and **never hardcoded or guessed**:

- `causal_gene` — `BUB1B`, decided at gate G1 on 2026-09-17 from L0's candidate call
  ([decision D9](mngmt/decisions.md)). A research premise, not a diagnosis; its uncertainties are
  recorded with it.
- `therapeutic_endpoint` — `chemoprevention`, meaning *secondary* prevention, decided at gate G2 on
  2026-09-08 ([decision D4](mngmt/decisions.md)). The admissible values:
  - `chemoprevention` — aneuploidy-buffering; remove pre-malignant clones
  - `symptomatic` — phenotype-directed
  - `mitotic_fidelity` — restore the SAC; largely undruggable, discuss-only

Individual channels are toggled under `channels:`.

## Design requirements this structure serves

- **Fairness / generalizability** — n=1 means hypothesis generation only; the pipeline must demonstrably
  generalize to other rare monogenic diseases (Scalability, 15% of the rubric).
- **Safety / reliability** — pediatric filter and genotoxic exclusion are hard gates in L4, not weights.
- **Transparency / explainability** — provenance and a written rationale accompany every candidate.
- **Security / privacy** — a real minor's genome: local only, never committed, no re-identification.
- **Reproducibility** — pinned environment, seeded RNG, config-driven, containerizable, and runnable by a
  judge against their own re-downloaded copy of the gated data.

## Repository layout

```
├── config/pipeline.yaml       # seed, paths, causal_gene (G1), endpoint (G2), channel toggles
├── docs/architecture.md       # full architecture write-up
├── src/
│   ├── pipeline.py            # orchestrator
│   ├── l0_genomics/  l1_target/  l2_channels/
│   ├── l3_integrate/          # RxCUI harmonization, rank aggregation, Claude reasoning
│   ├── l4_validate/           # safety triage, benchmark, sources.md (annotation catalog + licenses)
│   └── l5_report/
├── data/                      # gitignored — patient data, never committed
├── results/                   # gitignored — all pipeline output
├── notebooks/
└── tests/test_smoke.py
```

## Licensing of outputs

Outputs are **CC-BY-4.0**. Drug annotations are therefore built from license-clean sources only (openFDA,
RxNorm/RxClass/MeSH/MED-RT, UNII/NDC). NC and ShareAlike sources — DrugBank, SIDER, DDInter, ChEMBL,
PharmGKB, WHO ATC bulk — are **segregated into a non-redistributed enrichment zone**, joined by RxCUI.
See [src/l4_validate/sources.md](src/l4_validate/sources.md).

## Attribution

Sage Bionetworks · MVA Society · Hugging Face · BEACON · AWS · Anthropic · and above all the family and
the child whose data makes this work possible.
