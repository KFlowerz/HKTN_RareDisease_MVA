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

## D5 — LoF calls are made on MANE Select; LoF on other transcripts is kept, flagged (2026-09-11)

`annotator.transcript_policy: mane_select_tiered` in
[config/pipeline.yaml](../config/pipeline.yaml), implemented in
[src/l0_genomics/transcripts.py](../src/l0_genomics/transcripts.py).

**What it means.** L0 annotates twice. The **primary** pass uses only the MANE Select transcript
(`-tag MANE_Select`); a LoF consequence there is the call (`lof_tier = mane_select`). The
**secondary** pass uses every transcript; a variant that is LoF only there is kept as
`lof_tier = non_mane`, naming the transcript, and never becomes a call on its own.

**Why.** The review and data are in
[docs/research/transcript-policy-lof.md](../docs/research/transcript-policy-lof.md). In short:

- MANE Select is the transcript set designed as the default for clinical reporting (Morales et
  al., 2022), and all but 67 of 33,736 PubMed-supported pathogenic variants map to it (Pozo et
  al., 2022).
- LoF confined to a subset of a gene's transcripts is common and enriched for false positives
  (MacArthur et al., 2012; Singer-Berk et al., 2023), and PVS1 is not applied when the affected
  exon is missing from biologically relevant transcripts (Abou Tayoun et al., 2018).
- On this panel, 448 bp can receive a LoF call only from a non-MANE transcript
  (`results/transcript_policy/territory.tsv`, seed=42). Too much to promote to calls — but in an
  n=1 search the causal allele could sit there, so it is demoted rather than dropped.

**Why not the alternatives.**

- `-canon` — snpEff's canonical is the longest CDS, not MANE (SnpEff & SnpSift documentation,
  n.d.). On this panel it picks a different transcript for BUB1, BUB1B and CEP192, and leaves
  20 bp of MANE Select sequence in BUB1B unable to receive a LoF call (`picks.tsv`,
  `territory.tsv`). Refused by name in code.
- All transcripts as the call — promotes the partial-LoF territory above to calls.
- MANE Select alone — defensible, but silently discards the rare pathogenic allele that lies
  outside MANE, which is the costliest possible miss for one child.

**What it does not rest on.** Known pathogenic variants could not discriminate: all 118 public
ClinVar P/LP variants in the panel get the same call under every policy
(`results/transcript_policy/clinvar_summary.tsv`). The decision rests on the literature and on
the territory measurement, not on a benchmark win.

**Obliges.**

- Every L0 variant call is a `VariantCall` carrying `transcript`, `mane_release`, `lof_tier`,
  `annotation_pass` and `transcript_policy`. A `non_mane` tier can only come from the secondary
  pass and `mane_select` only from the primary — enforced in code and tested.
- `annotator.mane_release` names the MANE release the database's tags were checked against.
  Re-run `scripts/transcript_policy_check.py` on any snpEff database change.
- NMD escape — a stop in the last exon or the 3′-most 50 nt of the penultimate exon (Abou Tayoun
  et al., 2018) — is flagged, not treated as equivalent to an early truncation.
- L0 is not LoF-only. 6 of the 118 known pathogenic panel variants are not LoF under any policy;
  "no biallelic LoF found" is reported as exactly that, never as "no causal variant".
- A new policy needs a new decision here. Code refuses unknown values rather than defaulting.

---

## D6 — L0 cross-references ClinVar locally; the archive is evidence, not a classifier (2026-09-16)

`clinvar` block in [config/pipeline.yaml](../config/pipeline.yaml), implemented in
[src/l0_genomics/clinvar.py](../src/l0_genomics/clinvar.py). Output:
`results/l0_genomics/clinvar_crossref.json`, summarised onto every configuration in
`causal_gene_call.json`.

**What it means.** The whole public ClinVar GRCh38 release is downloaded once, cached outside the
repository, and matched against locally. Each panel allele is reported as an **exact** match, a
**same_position** record (a different allele at the same site — context, not identity), or
**none**. Records are matched on the minimal representation of `(contig, pos, ref, alt)` (Tan et
al., 2015), so the two files' padding conventions cannot hide a match. Pathogenic and
likely-pathogenic records in the panel regions are additionally annotated on the MANE Select
transcript, and those changing the same protein residue as one of the subject's alleles are
reported alongside it. Every significance travels with the submitter's review status.

**Why.** Gate G1 turns on a pair — a LoF allele and a protein-altering second allele — and
consequence prediction cannot weigh the second one at all. ClinVar is the public record of what
submitters have already concluded about these exact variants, with the evidence and review status
behind each submission (Landrum et al., 2018). Same-residue and same-change pathogenic records are
the observations behind ACMG/AMP criteria PM5 and PS1 (Richards et al., 2015), which is why they
are worth surfacing next to an allele.

**Why not the alternatives.**

- **A ClinVar API lookup per variant** — the obvious implementation, and refused. It would put
  this child's coordinates on an NCBI server, which [COMPLIANCE.md](../COMPLIANCE.md) prohibits
  for the same reason annotation runs on a local snpEff rather than Ensembl's VEP REST API. The
  whole-release download costs about 190 MB once and leaks nothing — not even a gene symbol.
- **Applying ACMG/AMP criteria automatically** — would make L0 a variant classifier. PVS1, PS1 and
  PM5 each carry conditions (transcript relevance, comparable amino acid change, splice effects
  excluded) that this pipeline does not evaluate, and the project is hypothesis generation, not
  clinical classification ([COMPLIANCE.md](../COMPLIANCE.md)). The observations are reported; the
  criteria are not applied.
- **A missense pathogenicity predictor (REVEL, CADD, AlphaMissense)** — would put a number where
  the missing evidence is population frequency, and an uncalibrated score reads as more certain
  than it is. Not adopted; the gap stays a stated caveat.

**What it does not settle.** Absence from ClinVar is not evidence of benignity — for a private
allele in a disease with roughly 50 known patients, absence is the expected state. Population
allele frequency is still unassessed: the frequency fields ClinVar carries (ESP, ExAC, 1000
Genomes) are legacy and absent from most records. Both remain caveats on every causal-gene call.
*(Population frequency: closed by D8, 2026-09-17.)*

**Obliges.**

- A reported significance always carries its review status and star rating. An unrecognised review
  status yields `null` stars, never `0` — a guessed zero is a claim of its own.
- The whole-release rule is asserted in `tests/test_l0_clinvar.py`, not just documented, because
  "just checking this one variant" is an easy way to break it later.
- The release's fileDate, byte size and SHA-256 are written into the artifact, so a cross-reference
  can be traced to the exact bytes. Re-run L0 after a ClinVar release change rather than reasoning
  from an older answer.
- A conflicting classification is never counted as pathogenic support — it is the one significance
  that names pathogenicity while meaning "submitters disagree".

---

## D7 — Drug–target data is used, not redistributed (2026-09-16)

Channel B scores drugs against the disease module, which needs each drug's protein targets.
Those target sets live in [src/l2_channels/enrichment.py](../src/l2_channels/enrichment.py)
and are never published; the ranking in `results/l2/channel_b_proximity/candidates.tsv`
carries identifiers, names and scores. Ledger entry: Table 5 in
[src/l4_validate/sources.md](../src/l4_validate/sources.md).

**What it means.** The target sets, action types and mechanism text stay in
`enrichment_dir` — a directory outside the repository, separate from the public reference
cache, so the line between what may be republished and what may not is physical rather
than conventional. Everything that reaches a published table passes through a whitelist,
`enrichment.PUBLISHABLE_FIELDS`: ChEMBL id, drug name, drug type, clinical stage, and the
*count* of targets. Identifiers and INN-style names are what Table 2 already permits for
performing a join.

**Why.** No comprehensive, gene-level drug–target resource is unambiguously clean for a
redistributed CC-BY-4.0 output. What the licence pages actually say (read 2026-09-16):

- **Open Targets** marks its Platform data CC0 1.0 and says downstream users may consume it
  without restriction — and, in a table on the same page, lists **ChEMBL as CC BY-SA 3.0**
  among its sources. The drug–target content is ChEMBL's.
- **Guide to PHARMACOLOGY** is ODbL plus CC BY-SA 4.0. **DrugCentral** is Creative Commons
  with a ShareAlike term. **DrugBank** is CC BY-NC. All already in Table 2.
- **RxClass / MED-RT** is genuinely clean (US-Gov, "no license is needed"), and is what
  Table 1 nominates as the DrugBank substitute — but see below.

ShareAlike is the viral one: bundling that content into a derived table would relicense
every output this project publishes ([COMPLIANCE.md](../COMPLIANCE.md)). Taking Open
Targets' CC0 mark at face value would have been faster and would have let the report name
each drug's targets everywhere. It was declined: the project's own ledger explains why
ShareAlike matters more than NonCommercial, and an output that contradicted that paragraph
would cost more in Scientific Rigor than the convenience is worth.

**Why not the clean source alone.** RxClass/MED-RT gives 781 mechanism-of-action classes
with real drug coverage, but they are *categories* — "Calcium Channel Antagonists", not
`CACNA1C`. Matching the class names against STRING's 2.5 million protein aliases resolved
**0 of 4** sampled classes, because a calcium channel is about twenty proteins rather than
one. Worse than the labour: even a perfect mapping yields a family, so amlodipine and
nifedipine would get identical target sets and Channel B would be ranking drug classes
rather than drugs. That data belongs in L4's safety triage and in per-candidate rationales,
not here.

**Obliges.**

- `enrichment.publishable` is the only sanctioned way out of the zone, and it is a
  whitelist: a field added to `DrugRecord` later stays inside until someone decides
  otherwise. `tests/test_l2_channel_b.py` asserts the written table carries no target
  identifier.
- The crosswalk from Ensembl gene id to interactome protein is built from STRING's own
  alias table (CC BY 4.0), so drug *identity* is the only thing taken from the restricted
  zone.
- A judge reproducing this downloads the release themselves — it is public and free. The
  release, its URL, byte size and SHA-256 travel in `channel.json`, marked
  `redistributable: false`.
- Where a published rationale needs to name a drug's target, it comes from the FDA label
  (openFDA, CC0), not from this source.
- Any future source of drug annotation gets a row in Table 5 and a zone before it is wired
  in. "It was already downloaded" is not a licence.

---

## D8 — L0 reports gnomAD population frequency, read by whole gene span (2026-09-17)

`gnomad` block in [config/pipeline.yaml](../config/pipeline.yaml), implemented in
[src/l0_genomics/gnomad.py](../src/l0_genomics/gnomad.py). Output:
`results/l0_genomics/gnomad_frequencies.json`, summarised onto every configuration in
`causal_gene_call.json`. Closes the gap D6 left open.

**What it means.** For every panel allele, L0 reports whether gnomAD v4.1.1 (Chen et al., 2024)
observed it in its exomes and its genomes, with allele count, allele number, frequency,
homozygote count and highest genetic-ancestry-group frequency. Where an allele was not observed,
the allele number of sites within 50 bp is reported beside the absence, so "not seen" can be told
apart from "not sampled". Only the five core CC0 fields are read.

**Why.** At G1 the open question about a protein-altering second allele was whether it is simply
common. Consequence prediction cannot answer that, and ClinVar's frequency fields are legacy and
mostly empty (D6). A population frequency is the direct measurement.

**Why not the alternatives.**

- **Download the whole release, as for ClinVar** — chromosome 15 alone is 24 GB across the two
  data sets (7.4 GB exomes, 17.0 GB genomes; server Content-Length, 2026-09-17), and the panel
  spans five chromosomes. Not proportionate to six genes.
- **Query gnomAD's API per variant** — refused for the reason D6 refused ClinVar's: it would put
  this child's coordinates on a third-party server.
- **Read by gene span (adopted)** — htslib fetches the tabix index and the compressed blocks
  covering each panel gene's whole padded span. The request is determined by the panel alone,
  which is already public in this repository, and is identical for every proband. A server log
  shows that someone read six mitotic-checkpoint genes, not which variants anyone carries.
- **A missense predictor instead** — still declined, for D6's reason: an uncalibrated score in
  place of a measurement.

**What it does not settle.** Rarity is necessary, not sufficient: most very rare variants are
benign, and no ACMG/AMP criterion (PM2 or other) is applied (Richards et al., 2015). Presence
does not imply benign either: healthy heterozygous carriers of a recessive pathogenic allele are
expected in any large population sample, so the homozygote count matters more. Phase is untouched: a rare pair on the same copy of a gene is still not
biallelic.

**Obliges.**

- `gnomad.extract` takes a gene span and nothing else, and `tests/test_l0_gnomad.py` asserts the
  queried regions do not change when the alleles do. A later "just look up this one position" has
  to break a test first.
- Only `AC`, `AN`, `AF`, `nhomalt`, `AF_grpmax` are read. SpliceAI scores in the same files are
  CC BY-NC 4.0 (Table 4b, [src/l4_validate/sources.md](../src/l4_validate/sources.md)).
- Each extract's URL, region, server ETag and Last-Modified, and SHA-256 travel in the artifact.
  A new gnomAD release means a new extract and a re-run, not reasoning from the old one.
- A failed read is fatal, like the ClinVar cross-reference: an empty extract must never read as
  "this gene has no variation", and the G1 artifact is not written without the evidence it claims
  to weigh. `gnomad.enabled: false` is the deliberate way to run without it, and the artifact then
  says frequency was not assessed.
- gnomAD's terms bind this project not to attempt to re-identify its participants.

---

## Open

- **`causal_gene`** (gate G1) — still `null`. A finding from L0, never a setting. The evidence is
  assembled: `results/l0_genomics/causal_gene_call.json` (seed=42) reports one gene with a LoF
  allele paired with a protein-altering second allele; D6's cross-reference
  (`clinvar_crossref.json`) says what ClinVar holds about each, and D8's lookup
  (`gnomad_frequencies.json`) how often gnomAD observed them. What the evidence cannot settle is
  recorded with it: the pair is unphased, and the second allele is absent from ClinVar. G1 is a
  person reading that and deciding.
- **Hackathon close date** — unknown, so the 30-day deletion deadline cannot be computed. See
  [docs/data_custody.md](../docs/data_custody.md).
- **APA title casing** in [docs/references.md](../docs/references.md) — Crossref preserves publisher
  casing; a sentence-case pass is owed before submission.
