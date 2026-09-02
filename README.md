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
| **L5** | [src/l5_report/](src/l5_report/) | Rubric-aligned report figures, candidate tables, and video assets. |

### L2 channels

| Channel | Module | Method | Why it is here |
|---|---|---|---|
| **A** | `channel_a_kg.py` | Knowledge-graph link prediction (Open Targets / PrimeKG / Hetionet + GNN, explainable meta-paths) | Anchored at the **gene** level, not the disease level — MVA is absent or near-empty in most KGs (cold-start defense). |
| **B** | `channel_b_proximity.py` | Network proximity (Guney/Barabási) between drug-target sets and the disease module | Does **not** require the disease to exist in a KG at all. |
| **C** | `channel_c_signature.py` | Signature reversion against a **PROXY** signature | There is **no patient RNA-seq**. Uses LINCS L1000 knockdown of the causal gene or a curated aneuploidy gene set; the substitution and its limits are documented in-module. |
| **D** | `channel_d_phenotype.py` | Phenotype/HPO-driven (Monarch / Orphanet / Open Targets) | Reaches symptomatic candidates that mechanism-first channels miss. |
| **E** | `channel_e_prior.py` | Aneuploidy-stress literature prior + agentic literature mining (RareAgent-style) | Captures published aneuploidy-tolerance biology no structured resource encodes. |

[src/pipeline.py](src/pipeline.py) is the orchestrator: it loads `config/pipeline.yaml`, seeds the RNG, and
runs L0 → L5 in order, writing to `results/` (gitignored).

---

## Status

**Scaffold only.** Every layer's `run()` is a documented stub that raises `NotImplementedError`. No
algorithms are implemented, no APIs are called, and no data is present.

## Quickstart

```bash
conda env create -f environment.yml
conda activate mva-track2

# data/ and results/ are gitignored, so they do not survive a clone — create them:
mkdir -p data results

# download the gated dataset into ./data — see DATA.md
python -m src.pipeline            # currently raises NotImplementedError at L0

pytest -q                         # smoke test: every layer imports, run() is callable
```

## Configuration

[config/pipeline.yaml](config/pipeline.yaml) drives everything. Two fields are deliberately `null` and
**must not be hardcoded**:

- `causal_gene` — set only after L0 / Track-1 reconciliation.
- `therapeutic_endpoint` — the primary open scientific decision:
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
├── config/pipeline.yaml       # seed, paths, causal_gene (null), endpoint (null), channel toggles
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
