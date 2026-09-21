# A safety-triaged repurposing shortlist for mosaic variegated aneuploidy

**Rare Disease, Real Kid: MVA Hackathon 2026 — Track 2 (Drug Repurposing)**

Seed 42 · causal gene `BUB1B` · endpoint *secondary* prevention · pipeline L0–L5
Every number below names the artifact under `results/` that produced it.

---

## Summary

Mosaic variegated aneuploidy has no disease-modifying therapy, roughly fifty patients
worldwide, and a cancer-predisposition phenotype that makes the usual repurposing
shortcuts dangerous. We built a six-layer pipeline that starts from one child's
whole-genome VCF and ends at a static dossier: 1,963 approved drugs nominated by parallel
evidence channels, **1,880 of them refused** by a paediatric safety triage, and 83
survivors presented in two tiers by the *kind* of evidence behind them.

The headline result is the refusal. Each of the 1,880 exclusions names the rule that
fired, the openFDA field it read, the SPL section, and the sentence it rested on. That is
the part of this work we would defend hardest, and it is the part a clinician could audit
line by line.

The second result is uncomfortable and we report it anyway. **Of the 83 survivors, 82
have no published evidence linking them to this disease at all.** They reach the shortlist
on network proximity — a number this pipeline computed, not a finding anyone has
published. Exactly one candidate carries a verified citation. We present those as two
tiers rather than one ranked list of 83, because merging them would invite a reader to
treat a computed score and a published result as the same kind of claim.

The third result is a null. We built the signature-reversion channel, calibrated it
against forty unrelated genes' knockdowns, found it carried **no gene-specific signal**,
and shipped it disabled. It nominates nothing.

**This is hypothesis generation.** Nothing here has been tested in a laboratory or a
clinic, and no candidate should reach a patient except through a clinician and a trial.

---

## 1. The problem, and what makes it hard

MVA is caused by biallelic loss of function in mitotic spindle-assembly-checkpoint genes —
commonly *BUB1B* (Hanks et al., 2004), also *CEP57*, *TRIP13*, *BUB1*, *BUB3* and
*CEP192*. Chromosome missegregation produces mosaic aneuploidy across tissues, with growth
restriction, microcephaly, developmental delay and a substantial cancer predisposition
(Scott et al., 2006). Management is symptomatic plus surveillance. There is no
disease-modifying therapy to repurpose *against*, and no trial population to validate in.

Four properties shaped every design decision:

**The disease is nearly absent from structured resources.** A knowledge graph keyed on the
disease node returns little or nothing. Any method that anchors on "MVA" as an entity will
cold-start.

**There is no patient RNA.** The dataset is raw reads plus a called VCF and a phenotype
document. Signature-based methods have no disease signature to reverse.

**There is no ground truth.** No drug has ever been shown to work here, so there is no
external benchmark and no way to measure whether a shortlist is *right* — only whether it
is defensible.

**The patient is a child with a cancer-predisposition syndrome.** A genotoxic agent is not
merely a poor candidate; it is the opposite of the therapeutic goal. Safety cannot be a
ranking weight.

### The endpoint

Gate G2 fixed the endpoint as **secondary prevention**: reducing recurrence and
second-primary risk in someone already at elevated risk. Not prevention of a first cancer,
and not treatment of an existing one. The distinction is load-bearing and the report states
it wherever the endpoint is named.

---

## 2. Architecture

Six layers, orchestrated by `src/pipeline.py`, config-driven, seeded, writing to
`results/`.

| Layer | Does |
|---|---|
| **L0** genomics | Ingest the VCF, annotate locally, call the causal gene panel, quantify aneuploidy burden from allele ratios |
| **L1** target | Causal gene → protein → interactome disease module; split into upstream and downstream target sets |
| **L2** channels | Five parallel candidate generators, deliberately independent |
| **L3** integrate | Harmonise drug identity, rank-aggregate, reason over the head with a local model |
| **L4** validate | Paediatric safety triage as hard gates, plus a blinded internal benchmark |
| **L5** report | Render the dossier; compute nothing |

The design bet is **multi-channel convergence**: any single repurposing method is fragile
on a hyper-rare disease, so several weak independent signals agreeing should carry more
than any one score. Section 6 reports what happened to that bet.

### What the channels are, and what they produced

| Channel | Method | Status | Ranked |
|---|---|---|---|
| **A** knowledge graph | Link prediction anchored at the *gene*, not the disease | **not built** | — |
| **B** proximity | Network proximity of drug targets to the disease module against a degree-matched null (Guney et al., 2016) | complete | 1,956 |
| **C** signature | Reversion of a proxy signature | **built, nominates nothing** | 0 |
| **D** phenotype | HPO-driven, via Monarch and Open Targets indications | complete | 25 |
| **E** prior | Curated aneuploidy-stress literature prior, verified against Europe PMC | complete | 8 |

*Source: `results/l3/integration.json`, `counts.per_channel_ranked` (seed 42).*

**Three of five channels produced evidence.** Channel A is unbuilt — the honest reason is
scope, not a finding. Channel C is built and produces nothing on purpose (Section 5). A
report claiming five channels would overstate what made this shortlist.

---

## 3. What the pipeline refused — the headline result

L4 applies safety as **hard gates, never weights**. A candidate failing any rule is
dropped, not down-ranked, because a ranked list is read as a recommendation however it is
captioned.

**1,880 of 1,963 candidates were excluded.**

| Rule | Excluded |
|---|---|
| Genotoxic or cancer-risk-increasing | 1,623 |
| Paediatric use not established | 256 |
| Not currently marketed | 1 |

Broken down by the finding that fired:

| Reason | n | What it means |
|---|---|---|
| `insufficient_evidence` | 1,231 | **A coverage fact, not a safety finding** |
| `safety_not_established` | 225 | Label does not establish paediatric safety |
| `clastogenic_finding` | 125 | Positive chromosomal-aberration finding in the label |
| `increased_tumours` | 91 | Carcinogenicity study reported increased tumours |
| `mutagenic_finding` | 54 | Positive mutagenicity finding |
| `malignancy_risk` | 51 | Label warns of malignancy risk |
| `positive_genotoxicity_assay` | 25 | Named assay positive |
| `cytotoxic_class` | 23 | Pharmacologic class is cytotoxic |
| `secondary_malignancy` | 19 | Secondary-malignancy warning |
| `aneugenic_finding` | 14 | **Aneugenic** — directly contraindicated here |
| others | 21 | carcinogenic, genotoxic statements, not indicated |

*Source: `results/l4/validation.json`, `counts.excluded_by_reason` (seed 42). Rendered as
`results/l5/figures/exclusions_by_reason.png` and the full table at
`results/l5/exclusions.html`.*

### Two things we want a judge to check

**`insufficient_evidence` is the largest bucket, and it is not a safety claim.** openFDA
describes drugs marketed in the United States; the channels nominate from a wider pool. A
candidate with no label is excluded because its status *cannot be established* — the
pipeline **fails closed**. Reading those 1,231 rows as "dangerous" would misrepresent the
triage, so the dossier says so on the exclusions page itself.

**Every exclusion is quotable.** Each row carries the openFDA field, the SPL section and
the sentence, trimmed but never paraphrased. This is only possible because openFDA is a
US-Government work in the public domain; a licence-restricted annotation source would have
made the exclusions unauditable. Of 262,883 label records, 54,930 matched a candidate
identifier.

That auditability is the strongest, most reproducible thing this pipeline does. It is also
the part most directly useful to a clinician, which is why it leads.

---

## 4. The survivors, in two tiers

83 candidates survived. They are presented **never as one ranked list**, because the
evidence behind them is not of one kind.

### Tier 1 — literature-supported: 1 candidate

**Trametinib.** Evidence grade *in vitro*, one verified citation
(doi:10.1038/s41467-024-52176-x), supported by two channels (prior and proximity),
computed confidence **0.132**.

The adversarial pass — which searches specifically for evidence *against* each candidate —
returned this, and we print it in full on the candidate's page:

> The evidence for trametinib's repurposing in mosaic variegated aneuploidy is primarily in
> vitro, and the in vitro evidence is limited to one record… the evidence grade is weak…
> Therefore, trametinib should not be shortlisted for this indication.

We left that in. A pipeline that searches for contradicting evidence and then suppresses it
when it lands is not doing the thing it claims.

### Tier 2 — network proximity only: 82 candidates

Every one reaches the shortlist because its targets sit near the disease module in an
interactome. **No published evidence connects any of them to this disease.**

These get **no generated rationale**. Eighty-two paragraphs explaining that there is no
evidence read like analysis and are not; the honest form is the absence stated once, as a
fact, on every row.

*Source: `results/l3/reasoning/rationales.json`, `results/l4/survivors.tsv` (seed 42).*

### Why this matters more than the ranking

The model answered `mechanism_supported: false` for all 83. That looked at first like a
degenerate field. It is the correct answer: there is no literature to support a mechanism
for 82 of them, and no contradicting record can be found in an empty evidence set either.

A shortlist of 83 drugs presented as one ranked table would be a more impressive-looking
deliverable and a less honest one.

---

## 5. The null result: channel C

Channel C scores drugs for *reversing* a disease signature. With no patient RNA, the query
must be a proxy: the LINCS L1000 consensus shRNA knockdown of *BUB1B* (GEO `GSE106127`),
pooled across all nine cell lines that carry one, scored against 46,601 Phase II compound
signatures (GEO `GSE70138`) covering 897 approved molecules, by the weighted connectivity
score (Subramanian et al., 2017).

We built it, calibrated it, and it nominates nothing. Four measurements, all reproducible
with `python scripts/channel_c_diagnostics.py`:

**The proxy is weak but real.** Nine cell lines, cross-cell-line Spearman ρ median 0.20
(range 0.02–0.47). Seven of nine fall below the conventional gold threshold
`distil_cc_q75 ≥ 0.2`. We keep and count them rather than filtering, because the
conventional filter leaves two cell lines — trading a weak pooled proxy for a weaker
single-lineage one.

**The proxy is half a confound.** It correlates **r = +0.52** with the mean L1000 compound
signature — the generic stress response every perturbation provokes. Uncorrected, the best
candidate does not beat a *random gene set* (p = 0.95), and the top hits are inert
compounds. Projecting that axis out is therefore required, not tuning.

**Two nulls failed before we found one that could.** A per-candidate query-permutation null
passed **60 of 60** candidates. A max-statistic over random queries passed **every** depth,
because a random query leaves 50.7% of signatures at exactly zero and the whole
distribution shifts. Both failed the same way: a structured biological query beats an
arbitrary gene set for reasons that have nothing to do with the gene. *A null that cannot
fail is not a null.*

**The null that isolates the question kills the channel.** Hold "knockdown query" fixed and
vary only *which gene*: 40 unrelated genes, each built and deflated identically.

| Depth | Observed | Control median | p |
|---|---|---|---|
| 1 | −1.2467 | −1.2063 | 0.366 |
| 10 | −1.1305 | −1.0704 | 0.244 |
| 50 | −1.0672 | −0.9137 | 0.195 |
| 250 | −0.9626 | 0.0000 | 0.146 |

*BUB1B is unremarkable at every depth.* Fifteen of the forty unrelated genes produce a
*stronger* best reversal. And the one pattern that looked like biology — colchicine,
albendazole, paclitaxel and ixabepilone at the top — is not gene-specific either: **15 of
40 control genes also surface a tubulin binder in their top 10**, including *POLE2* (a DNA
polymerase subunit) and *CRK* (an adaptor protein). Tubulin agents have extreme signatures;
that is all it was.

`generate()` raises rather than writing a table it cannot defend, and `channels.signature`
is `false` in the shipped config. The full reasoning is decision D19.

**What this costs us.** Channel C was the strongest remaining argument for convergence.
That argument is now closed negatively — see Section 6.

---

## 6. Does the central bet hold? No, and we show it

The architecture rests on cross-channel convergence. Here is what the channels actually
agreed on:

| Convergence class | Candidates |
|---|---|
| **discriminating** (≥2 channels that each rank a narrow slice) | **0** |
| mixed | 26 |
| none | 1,937 |

*Source: `results/l3/integration.json`, `counts.by_convergence` (seed 42).*

**Zero candidates are supported by two discriminating channels.** Proximity ranks 1,956 of
1,963 candidates — essentially every approved drug with a target in the interactome — so it
is classified as a *broad* channel and its agreement alone is not reported as convergence.
That threshold (`broad_channel_fraction: 0.5`) is in the config, not in the code, so a
reader can move it and see what changes.

This is the weakest point in the submission, and it is shown as a figure rather than
described in a footnote (`results/l5/figures/channel_contribution.png`). The honest
statement is: **this shortlist carries no discriminating cross-channel convergence, and we
do not claim any.**

---

## 7. Scientific rigor

### The blinded internal benchmark

There is no external ground truth, so we test whether the pipeline recovers a frozen
positive set of published aneuploidy-stress compounds.

| Scope | AUROC | 95% bootstrap CI | Positives found |
|---|---|---|---|
| Combined | 0.8147 | [0.5809, 0.9533] | 8 of 1,963 ranked |
| Proximity | 0.8003 | [0.6335, 0.9400] | 7 |
| Prior | — | — | **excluded as circular** |

*Source: `results/l4/benchmark/recovery.{csv,json}`, 2,000 bootstrap samples, seed 42.*

**The interval matters more than the point estimate.** Eight positives reach the ranking, so
0.81 is a number with a very wide envelope that touches near-chance at its lower bound. We
report the interval on the figure rather than a bar chart of the mean, because a bar chart
would imply a precision that does not exist.

**Channel E is excluded from its own benchmark.** Its positives *are* its seed file;
scoring it would measure whether a channel can find the list it was given. The headline is
channel B's recovery, because channel B never saw the set. Nothing here has been used to
re-tune a channel — a benchmark optimised against is a description, not a test.

**What it does and does not show.** It shows the ranking finds what the literature already
contains. It does not show the ranking finds what would help this child, and the positive
set is itself biased toward well-studied compounds.

### Evidence discipline

Every claim in this project is one of two admissible kinds: a **pipeline data result**
(artifact path + seed + config, regenerable) or a **cited source** in APA 7th with a
resolvable DOI. Every DOI was resolved against Crossref and checked for retraction and
correction notices; the verification log is in `docs/references.md` with the date each
check ran. An LLM-produced citation is treated as a hypothesis about the literature until
resolved — and dropped if it does not resolve.

### Reproducibility

Re-running a layer on unchanged inputs reproduces it **byte for byte**. Verified on L4 —
survivors, exclusions and the benchmark identical across runs, including the 2,000-sample
bootstrap interval — and on L5, where all four figures match by checksum.

| Pinned | How |
|---|---|
| Libraries | `environment.lock.yml` (exact linux-64 solve) and `environment.yml` (portable, fully pinned) |
| ClinVar | a **dated** weekly release, not the rolling `clinvar.vcf.gz`, with an archive fallback |
| STRING, Reactome, gnomAD, MANE, HPO, Monarch, Open Targets, LINCS | version or release date in `config/pipeline.yaml` |
| Every cached download | SHA-256 recorded beside the result |
| Model step | `temperature 0`, fixed seed, loopback-only endpoint |

`results/_manifest.json` records the seed, the SHA-256 of the config file, a digest of the
resolved settings with machine-specific paths excluded, the versions of the libraries that
decide the numbers, and whether hash randomization was actually fixed.

**Two things sit outside that envelope by nature, and we say so rather than claiming them
away.** Channel E queries Europe PMC live, so its literature reflects the index on the day
it ran — cached with query date and response hash. And the local model step is
deterministic in its *settings*; llama.cpp does not guarantee identical output across
builds and hardware, so that is a setting, not a promise.

---

## 8. Innovation

### Aneuploidy burden from allele ratios, with no alignments

MVA's phenotype *is* mosaic aneuploidy, so quantifying it is the natural innovation hook.
The dataset ships **no BAM** — raw reads plus a called VCF, verified at revision
`59e322d2…` on 2026-09-02. Depth-based methods have no input.

We estimate burden from **per-chromosome B-allele frequency in the VCF itself**: the mean
|BAF − 0.5| at heterozygous biallelic SNVs, in 10 Mb windows, at minimum depth 10, against
a baseline of the sample's own autosomal windows (median + scaled MAD). An allele ratio at
a single locus is self-normalising, so it needs none of the GC-bias and mappability
correction that depth would (Conlin et al., 2010; Loh et al., 2018).

The raw statistic sits at a noise floor when the true shift is zero and responds
quadratically to small shifts, so the mosaic-fraction estimator inverts the folded-normal
mean for the true BAF shift before converting — without that step, small fractions are
systematically overstated.

**Resolving power on this sample:** smallest detectable mosaic fraction **≈ 0.098**, with a
conservative bound of **≈ 0.245**, at z = 5.

### What we deliberately do not show

**The per-chromosome result is withheld — in every form, including aggregate counts.**

MVA case series publish karyotypes and the worldwide population is around fifty patients.
Which chromosomes are involved, how many, and at what mosaic fraction is close to a
fingerprint. A count of affected chromosomes is not a safe aggregate either: in a
population that small it narrows the field by itself.

This costs us the single most compelling figure the project could produce, in the criterion
where it would have counted most. We publish the **method and its resolving power** — which
is what the innovation claim is about, since the claim is *"allele ratios alone resolve
mosaic aneuploidy without alignments, down to ~10%"* — and not this child's result from it.

The withholding is **named, not silent**: `results/l5/report.json` lists the withheld fields
and the dossier prints the reason. Omitting it quietly would be indistinguishable from
never having built it. Decision D20.

### Enforced, not promised

`src/l5_report/publish_guard.py` fails closed on every string L5 writes — HPO term ids,
genomic coordinates, HGVS expressions, disease identifiers, karyotype strings, dbSNP ids.
Identifier patterns are *imported* from the L3 boundary rather than restated, so the two
cannot drift apart. Per-chromosome fields are checked structurally against mapping keys,
not as text — a first text-matching version rejected the report's own list of what it was
withholding.

The same guard enforces the endpoint wording **per sentence**: any sentence naming the
endpoint must qualify it as *secondary* prevention, because a qualifier in an introduction
does not travel with a table caption that a reader screenshots.

### Adversarial reasoning, locally

The L3 reasoning step runs a local open-weights model over a loopback socket — no patient-
derived content leaves the machine, enforced by a guard that refuses any non-loopback
endpoint. It makes two passes: one to synthesise a rationale from what the channels cited,
and one **adversarial pass that searches specifically for evidence against the candidate**.
The model may introduce no mechanism no channel cited, and every citation it emits is
validated against the channel's own reference set; unvalidated keys are dropped and
counted.

---

## 9. Scalability

**Partially demonstrated, and we are not going to overstate it.**

The Scalability claim is that the pipeline is disease-agnostic — that running it on a
different monogenic disease requires changing only the config. Nothing in L1–L4 hardcodes
MVA: the causal gene, the endpoint, the module parameters, every channel threshold and every
data release live in `config/pipeline.yaml`, and L1 refuses to run without a gene rather
than defaulting to one.

We ran a second disease — **cystic fibrosis, gene *CFTR*** — through **L1 and L2 channel B**,
changing only `l1.gene` in the config. Both layers ran unmodified.

| | MVA (`BUB1B`) | CF (`CFTR`) |
|---|---|---|
| Module seeds | 32 | 9 |
| Module size | 200 | 200 |
| Upstream / downstream split | 84 / 116 | 54 / 146 |
| Drugs scored by channel B | 2,566 | 2,566 |

*Sources: `results/l1_target/module.json`, `results/scalability_cftr/l1_target/module.json`,
`results/scalability_cftr/l2/channel_b_proximity/channel.json`.*

The interesting column is the split. A layer with MVA baked in would have produced the same
partition twice; instead the second gene seeds a different complex, expands to a different
module, and divides 54/146 rather than 84/116. The module layer is responding to the gene,
not reciting an answer.

**We did not run L3–L5 for the second disease, so the end-to-end claim is unproven.** Two
layers of six are demonstrated disease-agnostic; the rest are *written* to be. Completing
it is run time, not redesign — the remaining layers read the same artifacts in the same
shapes, and nothing in them is keyed on a gene name.

---

## 10. Governance, ethics and disclosure

This work involves a real child's genome. Several constraints were fixed before any code
was written and none was relaxed:

- **Patient data never enters the repository.** All data lives under gitignored `data/` and
  `results/`. Reference caches sit outside the working tree entirely, and the cache module
  *refuses* a path inside the repository.
- **No per-variant remote queries.** Every reference resource is downloaded as a whole
  public release and matched locally, so no coordinate of this child's ever reaches a
  third-party server. That is why annotation runs on a local snpEff rather than a REST API.
- **Clinical phenotype is patient data.** HPO terms are parsed from the data directory at
  run time and never written into config, source, tests or any committed file — a specific
  combination of features is identifying in a population this size. Tests use invented terms.
- **No recontact** with the subject, family or MVA Society, and nothing published beyond
  what the family already shares publicly through their own blog posts.
- **Licence segregation.** Non-redistributable sources (ChEMBL-derived drug-target content,
  and LINCS, whose licence we could not resolve at its canonical URL) are held in a separate
  cache and never enter a redistributed output. Only identifiers, names and derived scores
  cross that boundary, through a whitelist a test enforces.

**Generative AI disclosure.** The L3 reasoning step inside the pipeline runs a local
open-weights model (`qwen2.5-7b-instruct-q4_k_m`) over loopback; nothing patient-derived
leaves the machine. Development assistance used a commercially available assistant under a
consumer subscription with model training disabled. Both are recorded in decision D17 and
in `COMPLIANCE.md`.

**Data deletion** is committed in `COMPLIANCE.md` with a register of every custody location
in `docs/data_custody.md`. The attestation tooling (`src/purge.py`) is **not yet
implemented** — an outstanding obligation that outlives this submission.

---

## 11. Limitations

Stated plainly, because a shortlist for a child should arrive with its weaknesses attached.

1. **n = 1.** This is hypothesis generation for one patient. Nothing here is validated, and
   the pipeline cannot be validated on this disease because no therapy exists to validate
   against.
2. **No discriminating cross-channel convergence** (Section 6). The architecture's central
   bet did not pay off on this data.
3. **82 of 83 survivors have no published evidence.** Their support is a proximity score
   this pipeline computed.
4. **Two of five channels contribute nothing** — one unbuilt, one calibrated to a null.
5. **The safety triage is bounded by openFDA coverage.** 1,231 exclusions are
   "cannot establish", not "found unsafe".
6. **The benchmark rests on eight positives** and its confidence interval nearly reaches
   chance.
7. **The causal gene is a research premise, not a diagnosis.** `BUB1B` was fixed by a person
   at a recorded gate from ClinVar and gnomAD evidence; the allele pair is unphased and
   those uncertainties travel into every downstream output.
8. **End-to-end scalability is unproven.** Two layers of six are demonstrated on a second
   disease; the rest are written to be disease-agnostic but not shown to be (Section 9).
9. **Trametinib's own adversarial pass argues against shortlisting it** — and it is the one
   candidate with a citation.

---

## 12. What we would do next

**Complete the second-disease run** through L3–L5 — it is run time, not redesign, and it
converts Section 9 from two layers demonstrated to the whole pipeline.

**Try a curated aneuploidy-response gene set as channel C's proxy.** It is a genuinely
different query, derived from aneuploidy biology rather than one gene's knockdown, so it
need not lie along the generic drug-response axis that sank the knockdown proxy. The
scoring and calibration machinery already exists; it is a query swap and a rerun against the
same control-gene null. It may also come back negative.

**Build channel A** (gene-anchored knowledge-graph link prediction). It is the one
specified channel with no result at all, and it is the most likely source of a second
*discriminating* signal — the thing Section 6 shows is missing.

**Implement the deletion attestation** (`src/purge.py`), which is an obligation rather than
an improvement.

---

## Reproducing this

```bash
conda env create -f environment.lock.yml
conda activate mva-track2
snpEff download -noLog GRCh38.115
mkdir -p data results
# download the gated dataset into ./data — see DATA.md

PYTHONHASHSEED=42 python -m src.pipeline
# then open results/l5/index.html
```

`scripts/channel_c_diagnostics.py` reproduces every number in Section 5.
`mngmt/decisions.md` records each decision, what forced it, and what it obliged
downstream — including the ones that went against us.

---

## References

Full reference list with the Crossref verification log, including retraction and correction
checks and the date each ran: [`docs/references.md`](references.md). Machine-readable
source: [`docs/references.bib`](references.bib).

Key sources for the claims above:

Conlin, L. K., Thiel, B. D., Bonnemann, C. G., Medne, L., Ernst, L. M., Zackai, E. H.,
Deardorff, M. A., Krantz, I. D., Hakonarson, H., & Spinner, N. B. (2010). Mechanisms of
mosaicism, chimerism and uniparental disomy identified by single nucleotide polymorphism
array analysis. *Human Molecular Genetics, 19*(7), 1263–1275.
https://doi.org/10.1093/hmg/ddq003

Guney, E., Menche, J., Vidal, M., & Barábasi, A.-L. (2016). Network-based in silico drug
efficacy screening. *Nature Communications, 7*(1), Article 10331.
https://doi.org/10.1038/ncomms10331

Hanks, S., Coleman, K., Reid, S., Plaja, A., Firth, H., FitzPatrick, D., Kidd, A., Méhes,
K., Nash, R., Robin, N., Shannon, N., Tolmie, J., Swansbury, J., Irrthum, A., Douglas, J.,
& Rahman, N. (2004). Constitutional aneuploidy and cancer predisposition caused by biallelic
mutations in BUB1B. *Nature Genetics, 36*(11), 1159–1161. https://doi.org/10.1038/ng1449

Lamb, J., Crawford, E. D., Peck, D., Modell, J. W., Blat, I. C., Wrobel, M. J., Lerner, J.,
Brunet, J.-P., Subramanian, A., Ross, K. N., Reich, M., Hieronymus, H., Wei, G., Armstrong,
S. A., Haggarty, S. J., Clemons, P. A., Wei, R., Carr, S. A., Lander, E. S., & Golub, T. R.
(2006). The Connectivity Map: Using gene-expression signatures to connect small molecules,
genes, and disease. *Science, 313*(5795), 1929–1935. https://doi.org/10.1126/science.1132939

Loh, P.-R., Genovese, G., Handsaker, R. E., Finucane, H. K., Reshef, Y. A., Palamara, P. F.,
Birmann, B. M., Talkowski, M. E., Bakhoum, S. F., McCarroll, S. A., & Price, A. L. (2018).
Insights into clonal haematopoiesis from 8,342 mosaic chromosomal alterations. *Nature,
559*(7714), 350–355. https://doi.org/10.1038/s41586-018-0321-x

Scott, R. H., Stiller, C. A., Walker, L., & Rahman, N. (2006). Syndromes and constitutional
chromosomal abnormalities associated with Wilms tumour. *Journal of Medical Genetics,
43*(9), 705–715. https://doi.org/10.1136/jmg.2006.041723

Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., Natoli, T. E., Lu, X., Gould,
J., Davis, J. F., Tubelli, A. A., Asiedu, J. K., Lahr, D. L., Hirschman, J. E., Liu, Z.,
Donahue, M., Julian, B., Khan, M., Wadden, D., Smith, I. C., Lam, D., … Golub, T. R. (2017).
A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. *Cell,
171*(6), 1437–1452.e17. https://doi.org/10.1016/j.cell.2017.10.049

---

*Hypothesis generation only. Not medical advice, not a clinical recommendation, and not a
claim that any drug named here is safe or effective for any person.*
