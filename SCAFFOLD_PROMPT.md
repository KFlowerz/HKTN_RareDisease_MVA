# VS Code Scaffolding Prompt — MVA Hackathon 2026, Track 2 (Drug Repurposing)

**How to use this file**
- Paste everything under **PROMPT** into your VS Code AI agent (Claude Code, GitHub Copilot Chat, Cline, Cursor) and run it against the open `HKTN_RareDisease_MVA` folder.
- Or persist it as project context: save the knowledge-base sections as `CLAUDE.md` (Claude Code) or `.github/copilot-instructions.md` (Copilot), so every future request inherits the architecture and constraints.
- The prompt is idempotent: it must not overwrite existing governance files and must never create or commit patient data.

---

## PROMPT (copy from here to the end)

You are a senior ML / bioinformatics engineer scaffolding a **reproducible, compliance-first** repository for **Track 2 (Drug Repurposing)** of the *Rare Disease, Real Kid: MVA Hackathon 2026*. Your job in this task is **structure and stubs only** — create the directory tree, well-documented module stubs (clear docstrings, type hints, and `TODO`s tied to each layer's responsibility, raising `NotImplementedError`), configuration, and governance docs. **Do not implement algorithms, and never download, synthesize, or commit patient data.** Work inside the current repository root; if a file listed below already exists, leave it untouched.

### Knowledge base — Needs Analysis

- **Disease.** MVA (mosaic variegated aneuploidy, OMIM 257300): autosomal-recessive, biallelic loss-of-function in mitotic spindle-assembly-checkpoint genes (commonly **BUB1B**, **CEP57**, **TRIP13**; also **BUB1/BUB3/CEP192**) → chromosome missegregation → mosaic aneuploidy. Phenotype: growth restriction, microcephaly, developmental delay/seizures, cancer predisposition. No disease-modifying therapy exists; management is symptomatic + cancer surveillance.
- **Beneficiaries / stakeholders.** The one real child (subject); the ~50-patient MVA population and families; MVA Society (advocacy); treating clinicians (operators).
- **Operational need.** A credible, explainable, in-silico-prioritized shortlist of **approved** drugs worth investigating — hypothesis generation with genuine value given no existing therapy.
- **Primary function.** Given the causal gene/pathway + phenotype, output a **ranked, evidence-backed, safety-filtered** list of repurposing candidates, each with a mechanistic rationale and provenance.
- **Scope.** Computational nomination → prioritization → safety triage of approved drugs, with explanations.
- **Out of scope (enforce in docs and code comments).** Wet-lab validation; de-novo molecule design; any clinical recommendation; any recontact with the family; any efficacy claim.
- **Primary open decision (leave configurable, do not hardcode).** The therapeutic endpoint: `chemoprevention` (aneuploidy-buffering / remove pre-malignant clones) | `symptomatic` (phenotype-directed) | `mitotic_fidelity` (restore SAC — largely undruggable, discuss-only). Default `null` in config.
- **Judging rubric this repo serves.** Scientific Rigor 35%, Impact 25%, Innovation 25%, Scalability 15%. Single submission; deliverables are a written report + this repo + a 3-minute video.

### Knowledge base — Architecture (L0–L5)

Multi-channel evidence design (parallel weak signals → consensus), chosen because any single method is fragile on a hyper-rare disease. Each `src/l*` package gets a `run()` entry point and a docstring stating purpose / inputs / outputs / guardrail.

- **L0 `l0_genomics`** — Ingest WGS **VCF** (+ **BAM**); annotate variants (VEP/OpenCRAVAT-style), focus on biallelic SAC-gene variants; reconcile with the Track-1 validated causal variant. *Innovation hook:* use per-chromosome **BAM read-depth** to quantify aneuploidy burden as a feature (MVA's phenotype literally is mosaic aneuploidy). Output: causal gene + variant effect + affected pathway module.
- **L1 `l1_target`** — Map causal gene → protein/complex → interactome **disease module** (Reactome/STRING/Open Targets). Define an **upstream** target set (restore mitotic fidelity) and a **downstream** set (buffer aneuploidy stress / chemoprevention).
- **L2 `l2_channels`** — Parallel candidate generators, each a submodule:
  - `channel_a_kg` — knowledge-graph link prediction (Open Targets / PrimeKG / Hetionet + GNN; explainable meta-paths). Anchor at **gene** level, not disease level (cold-start defense).
  - `channel_b_proximity` — network proximity between drug-target sets and the disease module (Guney/Barabási). Does not require the disease to exist in a KG.
  - `channel_c_signature` — signature reversion on a **PROXY** signature (LINCS L1000 knockdown of the causal gene, or a curated aneuploidy gene set), because there is **no patient RNA-seq**. Document the substitution and its limits.
  - `channel_d_phenotype` — phenotype/HPO-driven (Monarch/Orphanet/Open Targets) for symptomatic candidates.
  - `channel_e_prior` — aneuploidy-stress literature prior + agentic literature mining (RareAgent-style).
- **L3 `l3_integrate`** — Harmonize drug identities on the **RxNorm RxCUI** backbone; rank-aggregate across channels (RRA/Borda), prefer cross-channel convergence; a **Claude-in-the-loop** reasoning step that synthesizes each candidate's rationale, actively searches contradicting evidence, and emits a calibrated confidence (uses the `anthropic` package).
- **L4 `l4_validate`** — In-silico validation + **pediatric safety triage**: mechanistic plausibility, druggability, pediatric use, **hard exclusion of genotoxic / cancer-risk-increasing agents** (MVA is cancer-predisposing), CNS/BBB penetration where wanted, clinical status. Plus a **blinded internal benchmark** (does the pipeline recover known aneuploidy/SAC-relevant compounds?) — there is no external ground truth for Track 2.
- **L5 `l5_report`** — Assemble rubric-aligned report figures, candidate tables, and video assets.
- **`src/pipeline.py`** — Orchestrator that loads `config/pipeline.yaml`, seeds RNG, and runs L0→L5 in order, writing to `results/` (gitignored).

### Knowledge base — Drug annotation sourcing (wire into L4)

Resolve every candidate to an **RxCUI** (RxNorm), then cross-map to DrugBank/ChEMBL/UNII/ATC/PubChem/NDC (UniChem). Build the redistributed layer only from **license-clean** sources and segregate everything else.

- **CC-BY-4.0-clean core (use freely):** openFDA drug label + FAERS + drugsfda (CC0); RxNorm / RxClass / MeSH / MED-RT (US-Gov public domain); UNII / NDC. This covers therapeutic class, MoA/PK class, indications, contraindications, warnings/boxed, precautions, adverse reactions, dosage/administration, dosage forms, and approval status — mostly as **prose** needing extraction (the Claude layer's job).
- **Must segregate (NC / ShareAlike — do NOT bundle in CC-BY-4.0 outputs):** SIDER (CC BY-NC-SA), DrugBank (CC BY-NC), DDInter (CC BY-NC-SA), ChEMBL (CC BY-SA), PharmGKB (CC BY-SA), WHO ATC bulk. Keep these in a non-redistributed `enrichment/` zone joined by RxCUI.
- **Currency gotchas to encode as comments:** FDA pregnancy A/B/C/D/X letters are retired (PLLR → narrative Section 8); the NLM RxNav drug-interaction API was discontinued 2024-01-02 (no structured, CC-BY-clean DDI-with-severity source exists — plan a licensed source or SPL-NLP).

### Non-functional requirements (reflect in structure and docstrings)

Usefulness/efficacy; **fairness/generalizability** (n=1 → hypothesis generation only, and show the pipeline generalizes = Scalability 15%); **safety/reliability** (pediatric filter, genotoxic exclusion); **transparency/explainability** (provenance + rationale per candidate); **security/privacy** (real minor's genome — local only, no re-identification); **reproducibility** (pinned env, seeds, config-driven, containerizable, runs on the judge's own re-downloaded gated data).

### Hard constraints (MUST)

1. **Never** create, download, or commit patient data or identifiable intermediates. All data lives under `data/` and `results/`, which are gitignored.
2. Do **not** use Git LFS for patient data.
3. Do **not** add NC/ShareAlike data sources to any redistributed (CC-BY-4.0) output path.
4. Governance files (`.gitignore`, `DATA.md`, `COMPLIANCE.md`) must exist with the exact content below if missing; never weaken them. For `.gitignore` specifically: if it already exists, verify it contains every required pattern and add any that are missing before proceeding — do not assume an existing `.gitignore` is correct.
5. Leave `causal_gene` and `therapeutic_endpoint` as `null` in config — do not invent them.

### Files to create with exact content (only if missing)

**`.gitignore`** — *highest-priority file.* Create it if missing. If it already exists, **do not skip it**: verify every pattern below is present and add any that are missing (this is the one file where a wrong version can leak a real minor's genome onto GitHub). Never remove or weaken existing ignore rules.
```
data/
results/
*.vcf
*.vcf.gz
*.bam
*.bai
*.cram
*.crai
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.g.vcf*
*.parquet
*.h5
*.npz
.env
*.token
hf_token*
.venv/
__pycache__/
*.pyc
.ipynb_checkpoints/
.DS_Store
```

**`DATA.md`**
```
# Data
No patient data lives in this repo, ever. The dataset is gated (~85 GB, 11 files,
WGS VCF/BAM + phenotype): https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data
Download into ./data (kept out of git by .gitignore).
TODO: enumerate the 11 files; confirm whether any expression data exists.
If none, Channel C runs on a PROXY signature — document that choice.
```

**`COMPLIANCE.md`**
```
# Compliance checklist
- No patient data / identifiable intermediates in git (see .gitignore). No Git LFS for data.
- No recontact of the subject, family, or MVA Society contacts.
- Delete all data within 30 days of Hackathon close; email confirmation to
  RarediseaserealkidMVAhackathon2026@synapse.org.
- Manuscript embargo until organizers post their summary/preprint; code/outputs shareable anytime.
- Required attribution block (Sage Bionetworks, MVA Society, Hugging Face, BEACON, AWS, Anthropic, family).
- Re-identification-avoidance: publish nothing that could re-identify the child/family.
- Outputs CC-BY-4.0; avoid NC / ShareAlike sources in redistributed outputs.
```

**`config/pipeline.yaml`**
```
seed: 42
data_dir: ./data
results_dir: ./results
causal_gene: null            # set after L0 / Track-1 reconciliation
therapeutic_endpoint: null   # chemoprevention | symptomatic | mitotic_fidelity
channels: {kg: true, proximity: true, signature: true, phenotype: true, prior: true}
```

**`environment.yml`**
```
name: mva-track2
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - pandas
  - numpy
  - scipy
  - networkx
  - scikit-learn
  - pyyaml
  - pysam
  - bcftools
  - jupyterlab
  - pip
  - pip:
      - anthropic
      - huggingface_hub
```

### Directory tree to produce

```
.
├── README.md                 # purpose, pipeline overview, "hypothesis generation only" disclaimer
├── DATA.md
├── COMPLIANCE.md
├── environment.yml
├── config/pipeline.yaml
├── docs/architecture.md      # placeholder: "paste the architecture analysis here"
├── src/
│   ├── __init__.py
│   ├── pipeline.py           # orchestrator: load config, seed, run L0→L5
│   ├── l0_genomics/          # __init__.py + run.py (VCF/BAM → causal gene; BAM depth → aneuploidy burden)
│   ├── l1_target/            # gene → interactome disease module
│   ├── l2_channels/          # __init__.py + channel_a_kg.py, channel_b_proximity.py,
│   │                         #   channel_c_signature.py, channel_d_phenotype.py, channel_e_prior.py
│   ├── l3_integrate/         # id harmonization (RxCUI) + rank aggregation + claude_reasoning.py
│   ├── l4_validate/          # safety_triage.py (exclude genotoxic) + benchmark.py + sources.md (annotation catalog)
│   └── l5_report/            # figures/report assets
├── data/.gitkeep             # gitignored dir, keep placeholder
├── results/.gitkeep          # gitignored dir, keep placeholder
├── notebooks/
└── tests/test_smoke.py       # imports each layer; asserts run() is callable / raises NotImplementedError
```

Each `src/l*` module: module docstring (purpose / inputs / outputs / guardrail), typed `run(config: dict) -> None` raising `NotImplementedError` with a `# TODO` describing the real implementation. `src/l4_validate/sources.md` should contain the drug-annotation source table (field → vocabulary → source → license → CC-BY-4.0 compatibility) from the knowledge base above.

### Acceptance criteria (self-check before finishing)

- The tree matches above; every `src` package imports without error.
- `python -c "import src.pipeline"` succeeds; `pytest -q` runs (smoke test passes).
- `.gitignore` exists and contains every required pattern (`data/`, `results/`, and all `*.vcf/*.bam/*.cram/*.fastq*` variants, `.env`, tokens); `DATA.md` and `COMPLIANCE.md` exist with the required content.
- No file exists under `data/` or `results/` except `.gitkeep`.
- `git status` shows only source/docs/config — **run `git ls-files | grep -E '\.(vcf|bam|cram|fastq)'` and confirm it returns nothing.**
- `config/pipeline.yaml` has `causal_gene: null` and `therapeutic_endpoint: null`.

### Do NOT

- Do not fetch, generate, or commit any genomic/clinical data.
- Do not implement channel algorithms, model training, or API calls in this pass — stubs only.
- Do not hardcode the causal gene or therapeutic endpoint.
- Do not add DrugBank/SIDER/DDInter/ChEMBL/PharmGKB-derived data to any redistributed output path.
- Do not overwrite existing governance files.

(End of prompt.)
