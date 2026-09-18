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

## D9 — Gate G1: `causal_gene` is BUB1B (2026-09-17)

Gate G1. `causal_gene: BUB1B` in [config/pipeline.yaml](../config/pipeline.yaml). Decided by the
maintainer on 2026-09-17, as drafted, after reviewing the L0 evidence below. L0 proposed; a person
decided — the pipeline never set this value itself.

The evidence: `results/l0_genomics/causal_gene_call.json` (seed=42), with D6's
`clinvar_crossref.json` and D8's `gnomad_frequencies.json`. Stated here at category level on
purpose: exact alleles, record ids and counts stay in those gitignored artifacts, because in a
population of roughly 50 patients a specific variant can identify a child.

**For** (pipeline results, seed=42, ClinVar fileDate 2026-09-13, gnomAD v4.1.1):
- Of the six panel genes, only one carries any PASS loss-of-function allele (`per_gene` in
  `causal_gene_call.json`).
- That allele truncates the protein on the MANE Select transcript, outside the last exon, and
  ClinVar holds it as pathogenic with two-star review status for mosaic variegated aneuploidy
  syndrome 1 (`clinvar_crossref.json`).
- The second allele is a missense variant in the protein kinase domain, residues 766–1050 of
  the 1,050-residue protein (The UniProt Consortium, 2026; matched to the MANE Select
  translation's length). It has no ClinVar record, and gnomAD reports it extremely rare and
  never homozygous, at a position both gnomAD data sets sampled
  (`gnomad_frequencies.json`).
- The pattern is the reported one. BUB1B was established as an MVA gene through truncating and
  missense mutations (Hanks et al., 2004), and in patients with biallelic mutations "a missense
  mutation pairs with a truncating mutation", with the missense change "consistently" in or near
  the BUBR1 kinase domain (Suijkerbuijk et al., 2010). Both quotations are from the abstracts;
  the full texts were not read.

**Against, or unsettled** — carried into the report wherever the call is used:
- **Phase.** The pair is unphased. With no parental samples and no recontact, a cis arrangement
  — both on one copy, the other intact — cannot be excluded. This is the largest remaining
  uncertainty and no available data can close it.
- **The missense allele is untested.** The literature supports the *pattern*; nothing shows
  *this* change impairs BUBR1. Rarity is necessary, not sufficient.
- **Unseen second hits.** Copy number, structural, deep-intronic and mosaic variants are not
  assessed from a called VCF (`caveats` in `causal_gene_call.json`).
- **No detectable constitutional aneuploidy.** `aneuploidy_burden.json` (seed=42) reports every
  chromosome as `baseline` or `below_detection_limit`. This does not exclude low-level mosaicism
  below the method's limit, but it is not supporting evidence, and the report must say so.

**What it does and does not mean.** It fixes the gene the drug search is organised around.
It is not a diagnosis, not a classification of either allele, and not a clinical statement
([COMPLIANCE.md](../COMPLIANCE.md)).

**Why decide now rather than wait.** No data this project may use can close the largest
uncertainty: phase needs parental samples or recontact, both excluded. Leaving G1 open would not
make the call more certain; it would leave every downstream layer running on the development
fixture instead of this child's disease.

**Obliges.**

- The four uncertainties above travel with the gene into the report, the video and every
  candidate page. The call is never described as a diagnosis or as "confirmed".
- L1 and every L2 channel are re-run against BUB1B; results produced on the development fixture
  are not reported as this child's.
- `tests/test_smoke.py` asserts the configured gene is a panel gene recorded in this log, and,
  where L0's artifact exists locally, that it is one of L0's candidates — so the value cannot
  drift from the evidence without a failing test.
- Revisit if L0 is re-run and its candidates change (a new ClinVar or gnomAD release, a panel
  change, a corrected VCF). A changed candidate list reopens G1; it does not silently update it.

---

## D10 — Channel D matches the phenotype locally; OMIM stays in, Mondo names come out (2026-09-17)

`l2.channel_d` in [config/pipeline.yaml](../config/pipeline.yaml), implemented in
[src/l2_channels/channel_d_phenotype.py](../src/l2_channels/channel_d_phenotype.py),
[phenotype.py](../src/l2_channels/phenotype.py) and [phenosim.py](../src/l2_channels/phenosim.py).
Output: `results/l2/channel_d_phenotype/`. Decided by the maintainer on 2026-09-17 after the
licence and liability review summarised in Table 6 of
[src/l4_validate/sources.md](../src/l4_validate/sources.md).

**What it means.** The subject's HPO terms are read from the phenotype document at run time and
never leave the machine. The HPO ontology (release 2026-09-01) and Monarch's Mondo-keyed
disease–phenotype table (release 2026-09-02) are downloaded whole and matched locally by
one-sided semantic similarity with an empirical p-value (Köhler et al., 2009; Resnik, 1999).
Approved drugs come from Open Targets' `clinical_indication` for the closest diseases, or for HPO
terms on the lineage of the subject's own terms, and stay in the D7 enrichment zone. Every
candidate is tagged `endpoint=symptomatic`.

**Why not the alternatives.**

- **Monarch's web API** — faster to build, and refused: it would send this child's feature
  combination to a third party, whose website states it runs HotJar and Google Analytics.
  Constraint 7 of `CLAUDE.md` already treats that combination as identifying.
- **FDA labels (openFDA) for the drug step** — licence-clean, but indications are free text,
  so the join from a matched disease to a drug would itself need extraction. Kept for L4's
  rationale text, per D7.
- **Monarch's full graph or its 3.6 GB semantic-similarity database** — the 4.4 MB
  association table already carries what the method needs, keyed by Mondo.

**Licences, as read.** The HPO licence is custom, not Creative Commons: use is free on condition
of citation, a displayed version, and no alteration. It is met by caching the files unaltered
and writing the version into `channel.json`. Mondo is CC BY 4.0; Orphanet's Science datasets
are CC BY 4.0. **OMIM's use agreement could not be read** (HTTP 403), and OMIM originates most
rare-disease annotations, so OMIM content is used for matching only and nothing OMIM-authored
is output: diseases are identified and named by Mondo.

**Liability, as read.** Monarch provides its data "as is", disclaims warranties and liability,
places "total and exclusive responsibility and risk" on the user, and says its tools "should not
be used for direct diagnostic use or medical decision-making". That matches this project's
scope (hypothesis generation, no clinical recommendation) and travels as a caveat.

**What it does not settle.** Coverage is thin by construction: diseases resembling a rare
presentation are mostly rare themselves, and few carry an approved indication recorded against
that exact Mondo disease. The first run on the real data found that most selected diseases
contributed no drug (`selected_diseases_contributing_a_drug` in `channel.json`, seed=42).
Widening the join to a disease's Mondo parents would add coverage and noise together.
**Decided by the maintainer on 2026-09-17: not widened.** Indications are joined to the exact
Mondo disease only; a drug indicated for a broader parent category is not treated as indicated
for the rare disease beneath it. Thin coverage is reported as thin, not padded.

**Obliges.**

- `evidence.json` is patient-derived and never published, quoted or committed: which diseases
  resemble this child, and which of the child's features an indication matched, are facts
  about the child. `candidates.tsv` and `channel.json` carry no HPO id, no disease, no term,
  and no value computed from the subject's terms; a test asserts that, and that the log
  carries none either.
- A negated row in the phenotype document stops the run rather than being guessed at, and
  errors name row numbers, never content. HPO's own term names are exempt: "Absent speech"
  names a feature that is *present*, and the ontology's wording may not veto its own row.
- **Negation in the prose around the table needs a person, and the review is recorded outside
  the repository.** No word list can separate "the following were excluded" from "no single
  feature in isolation is diagnostic" — a sentence this dataset's own notes contain. So the
  run stops until a reviewer sets `MVA_PHENOTYPE_REVIEWED` to the document's SHA-256 in the
  local environment file (`~/.config/mva/env`, never committed: it is a handle on a patient
  file). Editing the document changes the fingerprint and asks again. Reviewed for the
  current document on 2026-09-17: the two negation-carrying paragraphs are guidance on
  interpreting the phenotype, and describe no listed feature as absent.
- Any public display of channel D results carries the HPO version and acknowledgement from
  `channel.json`'s `attribution` field, and cites Gargano et al. (2024) and Putman et al. (2024).
- Releases are pinned by date in config; `latest` is refused by a test.
- Revisit the OMIM row in Table 6 once its agreement can be read.

---

## D11 — Channel E ranks a curated prior by verified literature, with no model in the loop (2026-09-17)

`l2.channel_e` in [config/pipeline.yaml](../config/pipeline.yaml), implemented in
[src/l2_channels/channel_e_prior.py](../src/l2_channels/channel_e_prior.py) and
[literature.py](../src/l2_channels/literature.py), seeded by
[aneuploidy_prior.tsv](../src/l2_channels/aneuploidy_prior.tsv). Output:
`results/l2/channel_e_prior/`. Decided by the maintainer on 2026-09-17, to reach gate G3
(three channels producing) without making the gate depend on an unevaluated model.

**What it means.** A committed seed file lists compounds the literature has tested against an
aneuploid or chromosomally unstable state. For each, the channel resolves its curated citations
against Europe PMC, retrieves the records where the compound and the aneuploidy vocabulary
co-occur in title or abstract, grades each record from NLM indexing
(`clinical` / `in_vivo` / `in_vitro` / `ungraded`), drops the retracted, and scores what remains
as a weighted count. Compounds resolve to a ChEMBL identity through the D7 enrichment zone so L3
can aggregate them with Channels B and D. Every candidate is tagged
`endpoint=chemoprevention` (gate G2, D4).

**Why no extraction model, for now.** CLAUDE.md's design for this channel is agentic literature
mining. That is still the intent, and the extraction step is kept behind an interface. But the
choice between a local model and the Claude API has not been evaluated, and G3 is due before it
can be. Shipping the retrieval-and-verification path first means the gate does not rest on an
unmeasured component, and the evaluation can be decided on numbers afterwards rather than under
deadline. A pinned local model would arguably be *more* reproducible than a server-side API for
the Scientific Rigor criterion — which is exactly the kind of claim that should be measured, not
assumed.

**Why Europe PMC and Crossref.** Europe PMC needs no key, has documented rate limits, and returns
retraction and correction status directly in `pubTypeList` and `commentCorrectionList` — verified
on a known retracted record. Crossref resolves DOIs independently and is public-domain metadata.
Checked on 2026-09-17: Crossref's `update-to` field was **empty** for that same retraction, so
Europe PMC is the primary retraction signal and Crossref the second opinion, not the reverse.

**Europe PMC's own terms make the REST API the sanctioned route.** Its copyright page prohibits
crawlers retrieving batches of articles from the website, and names OAI, RESTful, SOAP and bulk
download as "the only services that may be used for automated downloading of articles in Europe
PMC". This channel queries the RESTful service and fetches no page from the site. Article content
stays under publisher copyright with per-article licence terms that "are not identical for all the
articles", which settles the abstract question: since the terms vary per article and the pipeline
does not read each article's licence statement, no output reproduces abstract or full text at all.
The live page is unreachable from this machine (Cloudflare, HTTP 403, to `curl` and to a browser
user agent alike, and the Internet Archive's own 2026 captures recorded the same 403), so these
clauses were read from the Archive's 2025-06-05 snapshot and Table 7 names that date. An earlier
version of this decision and of Table 7 recorded the page as simply unread; that was corrected on
2026-09-17 once the snapshot was found.

**The outbound-query exception, and its bound.** Every other external read in this pipeline is a
whole-file download, because a per-item query would put the subject's data on someone's server
(`COMPLIANCE.md`). A literature index has no release to download, so this channel queries. The
exception is bounded in code, not by convention: `literature.refuse_private` rejects any query
carrying an HPO id, an OMIM id, a genomic coordinate, an HGVS expression, a dbSNP id or an allele
change, and it runs before the socket opens. Structurally, the module does not import the
phenotype reader at all, which `tests/test_l2_channel_e.py` asserts from the parsed module.

**What the first real run showed (seed=42).** Of 13 curated compounds, 8 ranked and 5 were
excluded — 3 unapproved, 2 with no ChEMBL identity. Two findings changed the design:

- **Co-occurrence in one record is not evidence about the two together.** Bortezomib initially
  ranked first on multiple-myeloma papers that name the drug in one place and chromosomal
  instability in another. Requiring both in a single sentence dropped its score from 13.0 to
  2.25 and removed lupus-with-trisomy-X and Down-syndrome-COVID papers from other compounds.
- **The lexical filters must not be applied to curated citations.** They were, and they discarded
  the best evidence in the file: the paper that identified chloroquine as aneuploidy-selective
  never names the compound in its abstract. A curated anchor is now checked for existence and
  retraction only, and `n_curated_citations_verified` says how much of each rank rests on that.

**Corrected after a multi-agent review, 2026-09-17.** Two of the seven findings changed what
the channel does rather than how it reads:

- **The wildcards in the config were dead.** `build_query` quoted every term, and Europe PMC
  does not expand a wildcard inside quotes: measured that day, `TITLE_ABS:"aneuploid*"` returns
  9,711 hits — exactly what `TITLE_ABS:"aneuploid"` returns — while the bare
  `TITLE_ABS:aneuploid*` returns 28,075. The query looked like it had a wildcard and did not,
  and it no longer meant the same thing as the local `term_pattern` check. Wildcard terms are
  now unquoted and phrases still quoted, since an unquoted phrase reads as separate terms. A
  term with both is refused: Europe PMC can express neither reading. **Retrieval widened and
  the ranking did not move** — chloroquine 1 → 4 retrieved records, reversine 12 → 17,
  tanespimycin 6 → 8 — because the sentence co-occurrence rule absorbed the additional hits.
  That the two changes cancel is evidence the precision filter is doing what it claims.
- **A Crossref outage was cached permanently.** A failed lookup was written to the response
  cache like any other answer, so a DOI Crossref could not resolve once would have stayed
  unresolved for the life of the cache, and a retraction deposited later would never be seen.
  Failures are no longer cached. Relatedly, a read timeout (`socket.timeout`, an `OSError` but
  not a `URLError`) escaped the handler and would have killed the channel.

The other five: a public AACR-style DOI tripped the patient-data guard, because
`10.1158/0008-5472.CAN-13-1174` ends in something shaped like chromosome 13 at position 1174
(DOIs are now masked before the coordinate and HGVS checks — a guard that cries wolf on a
citation gets widened until it stops guarding); two molecules sharing a preferred name were
reported as unambiguous and resolved by Parquet row order; a record with no PMID, DOI or
Europe PMC id could be counted as support although no reader could look it up; and the
channel-level counts of withdrawn and adverse records summed over ranked compounds only,
hiding exactly the compounds whose citations were dropped.

**What it does not settle.** After filtering, most ranked compounds rest on their curated
citation alone (`compounds_supported_only_by_curation` in `channel.json`). The candidate set is
curated, so the channel verifies and ranks what a person chose; it discovers nothing. Direction
is decided by fixed phrases and cannot tell a finding from a proposal. Grade weights are a
declared choice, not a measurement. Each of these is a caveat on the channel's output, and each
is a thing the model evaluation is meant to move.

## D12 — Hackathon close is taken as submission close, 2026-10-24 (2026-09-17)

Recorded by the maintainer on 2026-09-17 from the organizers' published timeline. Tentative:
the timeline itself says "Dates are subject to change. Any updates to key milestones will be
reflected in this table posted to the Community page."

**The organizers' milestones.** Launch 2026-08-24; submissions open 2026-08-25; **submissions
close 2026-10-24, 23:59 UTC**, Track 1 leaderboard frozen; Track 1 qualitative evaluation and
Track 2 expert-panel judging 2026-10-24 to 2026-11-24; winners announced 2026-11-25.

**What it means.** Two things that had been assumptions become facts:

- **The submission deadline is 2026-10-24, not 2026-10-05.** The delivery plan's target was
  labelled an assumption in the plan itself ("the submission date is an assumption, not a
  repository fact"). It was 19 days early. The plan's W1–W5 structure now ends almost three
  weeks before the real deadline.
- **The deletion deadline is 2026-11-23 23:59 UTC**, 30 days after submission close, which
  `COMPLIANCE.md` commits this project to. `docs/data_custody.md` carries it, and the purge
  attestation's `hackathon_close_date` and `deletion_deadline_utc` are filled in.

**Why submission close and not the announced end of the event.** "Hackathon close" is not
defined by the timeline, and the two readings are a month apart: 2026-11-23 from submission
close, or 2026-12-25 from the 2026-11-25 announcement. The earlier one is taken because this is
an obligation this project made about a child's genome, and where the wording is ambiguous the
reading that deletes sooner is the one to be held to. Choosing the later date would be choosing
to hold patient data for an extra month on our own interpretation.

**What it costs, stated rather than discovered later.** Judging runs to 2026-11-24 and winners
are announced 2026-11-25, so on this reading the data is deleted *during* the judging window,
one day before it ends. If the expert panel asks something that needs the pipeline re-run
against the real data, the answer will be that the data is gone and why. That is the correct
answer and it should not be avoided by quietly extending the deadline. The way to change it is
to ask the organizers which date they mean and record their answer as a new decision.

**What it does not settle.** Whether the extra 19 days change the plan. More time is not more
scope by default: L3, L4 and L5 are all still scaffolds, and the gates they serve (G4, and the
submission sweep) were sized for a schedule that assumed less time, not more. Any replan is its
own decision.

## D13 — Gate G3 passed: three channels producing independently (2026-09-18)

Signed by the maintainer on 2026-09-18, against the L2 run that finished 2026-09-17 22:01 UTC.

**The evidence.** `results/l2/channels.json`, seed 42, `causal_gene: BUB1B`:
`complete: [phenotype, prior, proximity]`, `n_complete: 3`. Channel B wrote 2,568 ranked rows,
channel D 28, channel E 8, each with per-candidate provenance. `kg` and `signature` are recorded
`not_implemented` rather than empty, which is the distinction the runner exists to preserve — an
empty list is a scientific claim, an exception is not.

**What the signature covers.** That the three channels are independent enough for L3 to
rank-aggregate across them, and that each one's weaknesses travel into the report rather than
being smoothed away.

**What it explicitly does not cover: that the candidates are good.** They are not yet, and the
sign-off was made with these three facts in front of it:

- **Channel B's top five are not five drugs.** Bortezomib, bortezomib D-mannitol, ixazomib,
  ixazomib citrate, carfilzomib — salt and co-crystal forms competing with their own parents.
  Open Targets records `parentId = CHEMBL325041` on the D-mannitol form, so L3's identity
  harmonisation can collapse them. All five also carry `n_targets_in_module = 0`: proximity is a
  distance measure, so that is legal, but the top of B is not "hits the disease module".
- **Channel D reaches antileishmanials** (miltefosine, amphotericin B) through the disease route.
  The method matched a disease whose *annotated features* resemble the phenotype, not one the
  child has. That is D10's similarity-is-not-identity caveat firing, visibly, in the output.
- **Bortezomib sits near the top of both B and E**, which will read as convergence. It is
  cytotoxic, and channel E's own caution column already calls it unsuitable for chronic
  prevention. For a cancer-predisposed child it is close to the worst possible nomination.

**What this obliges.**

- **L4's genotoxic and cytotoxic exclusion binds hard**, and the compounds above are the first
  test of it. G4 publishes the exclusions table as a headline result, not an appendix.
- **L3 harmonises identity before aggregating**, on the RxCUI backbone, so a molecule and its
  salt form cannot occupy two ranks (D2's field contract already requires this).
- **Convergence is not counted naively.** Channel B scores 2,566 of the approved drugs it could
  reach, so a high rank there is weak evidence alone and agreement between B and anything else is
  close to uninformative. Agreement between D and E, or a drug placing well in B *and* appearing
  in both others, is the signal worth weighting. L3 must say which kind of convergence produced
  each rank rather than reporting a bare channel count.
- **The independence claim has a ceiling, and the report states it.** B, D and E all reach drug
  identity and indications through Open Targets/ChEMBL. Three channels agreeing is not three
  independent pieces of evidence when they share a substrate.

**What it does not settle.** Channels A and C remain unbuilt. G3's condition is three, not five,
and the two stubs are recorded as such — but cross-channel convergence over three channels that
share a drug substrate is a thinner claim than the architecture was designed to make, and the
report says so rather than implying five-channel consensus was achieved and trimmed.

## D14 — L3 harmonises on the parent molecule and aggregates by RRA; RxCUI is attached, not required (2026-09-18)

Implemented in [src/l3_integrate/harmonize.py](../src/l3_integrate/harmonize.py),
[aggregate.py](../src/l3_integrate/aggregate.py) and [run.py](../src/l3_integrate/run.py).
Output: `results/l3/`. Licences in Table 8 of
[src/l4_validate/sources.md](../src/l4_validate/sources.md).

**Identity: the primary key is the parent ChEMBL molecule, not the RxCUI.** `CLAUDE.md`
specifies RxCUI as the backbone, and it remains the cross-system identifier every candidate
carries where one exists. It is not the *key*, because it cannot be: of 2,575 ChEMBL ids
across the three channels, 45% resolve to no RxCUI at all from the public sources this
project can use. Making RxCUI the key would silently drop nearly half the pipeline —
mostly drugs approved outside the United States — while looking like plumbing.

The key is therefore Open Targets' `parentId`, which is complete, local, and a curated
assertion rather than our string matching. A second pass merges identities whose *whole*
RxCUI set is identical and non-empty, which catches formulation pairs `parentId` does not
record. On the real data that is not hypothetical: ixazomib and ixazomib citrate arrived as
two identities holding identical RxCUI sets and took ranks 2 and 3. 2,575 ChEMBL ids
collapse to 1,963 molecules — 612 duplicates that would otherwise have competed with their
own parents for rank and split the convergence signal.

**Crosswalk: ChEMBL → UNII → RxCUI, from two bulk files.**

- **Not RxNav**, which would answer this in one call per drug. Channel D's candidates are
  derived from the subject's phenotype, so the *set* of drugs asked about is itself a weak
  statement about the child. A bulk file is asked about nothing. This is the pipeline's
  standing rule applied to the one case where the convenient route is a per-item query.
- **Not RxNorm's full release**, which needs a UMLS licence key. A judge rebuilding from
  this repository alone could not run it.
- **UniChem has no RxNorm source** — checked 2026-09-18, 25 sources, none of them RxNorm —
  so the route goes through UNII, which both UniChem and openFDA carry.
- Resolution route is recorded per candidate: 865 by UNII, 251 by name, 854 unresolved. A
  name match is a weaker claim than a registered-substance match and is labelled as one.

**Aggregation: Robust Rank Aggregation** (Kolde et al., 2012), over ranks normalised within
each channel. A candidate absent from a channel contributes nothing rather than a penalty —
channel E ranks 8 compounds from a curated file, and the 2,560 drugs it never considered
are not evidence against them.

**Convergence is reported by kind, not as a count.** A channel that ranks more than half the
candidate pool cannot discriminate, whatever its internal score says. Channel B ranks 2,566
of 1,963 molecules' worth of the pool and is automatically flagged *broad*; agreement with
it is labelled `mixed`, agreement between two narrow channels `discriminating`. Reporting a
bare "supported by 2 channels" would let the weakest possible agreement read as the
strongest.

**The first real run produced a result that matters more than the ranking (seed 42).**

> **No candidate has `discriminating` convergence. Channels D and E do not overlap at all.**

26 rows are `mixed` — each is one narrow channel plus channel B — and 1,937 are single-
channel. The architecture's headline claim is cross-channel convergence, and at the level
that would carry weight, there is currently none. The explanation is not a bug: channel D
answers a symptomatic question and channel E a chemoprevention one, so they are ranking
against different endpoints and are not expected to agree. But it means the report cannot
present convergence as an achieved result on three channels, and must say this plainly.
It is also the strongest argument yet for building channel C, which would be the first
additional narrow channel aimed at the same endpoint as E.

**What the top of the table says, and why it is not alarming.** The highest-ranked rows are
proteasome inhibitors and cytotoxic nucleoside analogues — bortezomib, ixazomib,
carfilzomib, gemcitabine, fludarabine, clofarabine, pemetrexed — plus a cluster of HIF-PHI
anaemia drugs. Channel D contributes withdrawn anti-obesity agents, among them **lorcaserin,
withdrawn in 2020 over an increased incidence of cancer**. None of this has been through a
safety filter, because L3 ranks and does not judge. It is L4's first and most important
test, and the exclusions table is the result that belongs in the report beside the survivors.

It also exposes a limitation to encode in L4: Open Targets' `maximumClinicalStage` of
`APPROVAL` means *was approved*, not *is marketed*. Withdrawn drugs carry it.

**What this does not settle.** The Claude reasoning step is deliberately not built yet, so
D2's per-candidate field contract is **not** satisfied — `rationale`,
`contradicting_evidence`, `confidence` and the safety fields are all outstanding.
`integration.json` records the gap explicitly rather than letting L5 discover it at render
time. Keeping the reasoning step separate is the point: it costs API spend, and a failure
in a later stage must never force it to re-run.

## D15 — L4 triages on openFDA label text, fails closed, and the exclusions are the result (2026-09-18)

Implemented in [src/l4_validate/safety_triage.py](../src/l4_validate/safety_triage.py),
[labels.py](../src/l4_validate/labels.py), [benchmark.py](../src/l4_validate/benchmark.py)
and [run.py](../src/l4_validate/run.py). Output: `results/l4/`. Sources in Table 9 of
[src/l4_validate/sources.md](../src/l4_validate/sources.md).

**Three hard rules, in order: genotoxicity, paediatric use, currently marketed.** Each is an
exclusion, never a penalty, and no channel-consensus score overrides one. Two further rules
are recorded rather than scored: blood-brain-barrier penetration is `not_applicable` under
this endpoint (D4), and clinical status quotes whatever boxed warning or contraindication
the label carries.

**Genotoxicity is decided from three independent signals**, any one of which excludes: the
FDA's own Established Pharmacologic Class vocabulary, an affirmative finding in SPL section
13.1, and an increased-malignancy statement in a boxed warning. More than one source
matters because section 13.1 is absent from many labels that a pharmacologic class alone
would condemn.

**The headline result is an exclusion, not a nomination.** Of 1,963 candidates, **83
survive and 1,880 are excluded**. The largest single reason is `insufficient_evidence`
(1,231): openFDA describes drugs marketed in the United States, and the channels nominate
from a much wider pool. That is a coverage fact and is reported as one, never disguised as
a safety finding. Only 5 survivors carry support from more than one channel, and
**trametinib is the only one of those on a mechanistic axis** -- the other four reach the
list through channel D's symptomatic route.

**Selumetinib is excluded by its own label, and this is the most important thing L4 found.**
The curated prior called it the most paediatric-ready compound on the RAF/MEK/ERK axis.
Its FDA label says: *"Selumetinib did result in an increase in micronucleated immature
erythrocytes (chromosome aberrations) in mouse micronucleus studies, predominantly via an
**aneugenic** mode of action."* An aneugenic drug for a child whose syndrome is mosaic
variegated **aneuploidy** is the worst nomination this pipeline could make, and the
curated prior would have made it. The safety layer caught it from primary label text.

**Trametinib survives on the same axis**, with the same channel-E evidence behind it
(Zerbib 2024): *"Trametinib was not genotoxic in studies evaluating reverse mutations in
bacteria, chromosomal aberrations in mammalian cells…"*, and paediatric use established
from 1 year of age. So L4's output is not "the MEK axis is unsafe" but "on the MEK axis,
trametinib is the defensible candidate and selumetinib is not" — a conclusion no channel
reached and the curated prior got backwards.

**The dose caveat is carried, not resolved.** Selumetinib's aneugenicity is reported at
doses roughly 38 times the clinical Cmax. The rule excludes on the finding regardless of
that margin; a pharmacologist would weigh it. The exclusion therefore ships with the
sentence that produced it, dose context included, so a reader can disagree with the rule
rather than having to trust it.

**Negation handling is the difference between those two verdicts**, and it took two
iterations to get right. An affirmative-pattern match is discarded when a negation cue
precedes it in the same clause, with clause scope ending at sentence punctuation *and* at
contrast conjunctions — "negative in the Ames assay **but** was clastogenic" is a positive
finding. Commas are deliberately not boundaries, because "not carcinogenic, mutagenic, or
clastogenic" is one negation governing three terms. The first version had no negation
handling at all and excluded selumetinib on the sentence saying it was *not* clastogenic —
the right answer for the wrong reason, which is worse than a wrong answer because it looks
correct.

**The benchmark is blinded by excluding the channel that cannot be blinded.** The positive
set *is* channel E's seed file, so channel E's recovery of it is 100% by construction and
measures nothing; it is reported as `circular` and kept out of the headline. So is the
combined ranking, which inherits channel E's contribution.

> **Channel B — network proximity, which never saw the positive set — recovers 7 of 18
> known aneuploidy-selective compound names with AUROC 0.80 (95% bootstrap CI
> 0.63–0.94).**

That is the project's first independent validation signal: a method that knows nothing
about the aneuploidy literature ranks its compounds above chance. The interval is wide
because the set has 18 names covering 13 compounds, and the point estimate is not a
measurement on its own.
Channel D recovers none, which is expected — it ranks against a symptomatic endpoint.

**What this obliges.**

- The exclusions table ships with the results and is a headline result in the report, not
  an appendix. What the pipeline refused to propose for a cancer-predisposed child is as
  informative as what it nominated, and here it is more so.
- No rule may be relaxed to increase the survivor count. If 41 is too few, the answer is
  more channels or better identity coverage, not a softer gate.
- Nothing in the benchmark may be used to re-tune a channel without re-freezing the
  positive set and disclosing the reuse.

**Corrected after a multi-agent review, 2026-09-18.** The review found eight issues in the
committed L3 and L4 code, and the first of them meant this decision's original numbers
described a filter that was not doing what it claimed.

- **The hard genotoxic gate was failing open.** The negation handling cancelled an
  affirmative finding whenever *any* negation cue appeared earlier in the clause, so
  ordinary label prose cleared real findings: "Patients with **non**-Hodgkin lymphoma had
  an increased risk of secondary malignancies", "Although **no** increase in tumors was
  seen at low dose, the drug was clastogenic", "In **non**-clinical studies the compound
  was clastogenic", "Inventib, which has **no** effect on fertility, was carcinogenic in
  mice". Cue-scanning is replaced by explicit negated *spans* that must **cover** the
  matched term. Failing open is the one direction this rule must never fail, and it was.
- **The paediatric rule excluded drugs that are approved for children.** Almost every
  paediatric approval states its own lower bound as a denial -- trametinib is approved from
  1 year of age and the same section says "have not been established ... less than 1 year
  old" -- and checking denials first let that floor veto the approval. Establishment is now
  checked first, with denial spans masking any age band inside them so a denial cannot
  clear a drug either.

  This changes what a passing paediatric verdict *means*, and the caveats now say so: it
  means the label establishes use in **some** paediatric band, which is quoted in the
  verdict. Whether that band covers this child is an age question this layer cannot answer,
  because the proband's age is patient data and never enters the pipeline.

The counts above are from the corrected run. For the record, the three states were 41
survivors (both bugs present), 14 (genotoxic gate fixed), and 83 (both fixed). The first
figure was the one originally recorded here, and it was wrong in both directions at once.

The remaining six findings did not change a verdict but each would have: merging every
matching label instead of an arbitrary first one, AUROC emitting a bare `NaN` that made
`recovery.json` unparseable, bootstrap intervals falling outside [0, 1] from resampling
ranks as though distinct, the positive set ignoring aliases and under-reporting recovery,
a clearing verdict quoting the section's opening rather than the sentence that cleared it,
and blank drug names when a parent id was absent from the molecule release.

**What this says about the process.** Both behaviour-changing bugs were in code I had read
line by line and believed correct, and one of them was introduced *by* my own earlier fix
to the same function. The review is not optional for this layer.

**What it does not settle.** Drug interactions are not assessed at all: the NLM RxNav
interaction API was discontinued on 2024-01-02 and no structured, licence-clean source with
severity exists, so this layer makes no interaction claim rather than an unsourced one.
The rules are lexical and cannot read a label the way a pharmacologist does. And the
survivors are still unreasoned — D2's `rationale`, `contradicting_evidence` and
`confidence` fields await L3's Claude step.

## Open

- **Which date the organizers mean by "Hackathon close"** — taken as submission close (D12), which is
  the earlier and safer reading. Worth asking them, because the alternative is a month later and the
  chosen deadline falls one day before judging ends. See [docs/data_custody.md](../docs/data_custody.md).
- **Whether the extra 19 days change the plan** — the delivery plan still ends 2026-10-05 and the
  remaining 2.7 weeks are unallocated on purpose (D12).
- **APA title casing** in [docs/references.md](../docs/references.md) — Crossref preserves publisher
  casing; a sentence-case pass is owed before submission.
