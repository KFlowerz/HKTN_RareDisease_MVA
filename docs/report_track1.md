# Identifying the causal variant pair in a child with mosaic variegated aneuploidy

**Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1 (Variant Prediction)**

Seed 42 · GRCh38 · snpEff GRCh38.115 · ClinVar release 2026-09-13 · gnomAD v4.1.1
Every number below names the artifact under `results/` that produced it.
Every term, abbreviation and acronym is defined in [`glossary.md`](glossary.md).

> **A note on what this document does not contain.** No variant coordinate, HGVS expression
> or ClinVar accession appears here. The submission CSV carries the identifiers, as the
> format requires; this report describes the method and the evidence at category level. In a
> disease with roughly fifty patients worldwide, a specific variant can identify the person
> it came from, and the readable document is the one most likely to be forwarded. The
> reasoning is recorded as decision D21 in the repository.

---

## Summary

I submit a **compound-heterozygous pair in *BUB1B*** as the primary finding: one truncating
allele that ClinVar holds as pathogenic at two-star review status for mosaic variegated
aneuploidy syndrome 1, paired with a rare missense variant in the BUBR1 kinase domain.

The call was produced by the L0 layer of a pipeline built for Track 2, where the causal gene
is the premise the drug search is organised around. Identifying it was not a side quest — it
is the gate the rest of the work depends on, and it was decided by a person reviewing the
evidence, not by the pipeline setting its own premise.

**The pair is unphased and I say so in the submission itself.** With no parental samples
and no recontact permitted, a *cis* arrangement cannot be excluded. That is the largest
uncertainty in this call and no data available to me closes it.

**This is a computational prediction, not a diagnosis**, and not a classification of either
allele. Section 8 states what that costs.

---

## 1. The problem, and what makes it hard

MVA is caused by biallelic loss of function in the genes that build the mitotic
spindle-assembly checkpoint — most often *BUB1B* (Hanks et al., 2004). Chromosome
missegregation produces mosaic aneuploidy across tissues, with growth restriction,
microcephaly, developmental delay and a substantial cancer predisposition (Scott et al.,
2006). Roughly fifty patients are known worldwide.

**The patient is a child.** That is stated once, here, because it bears on the work: it is
why Track 2's safety triage is paediatric, and part of why the privacy constraints in
Section 7 are as strict as they are. Everywhere else this report says *patient*.

Four properties of this case shaped every decision that follows:

**No parental samples, and no recontact.** Phase — which parental copy each variant sits on
— is the single piece of evidence that would settle a recessive call, and it is the one
piece that cannot be obtained. No analysis recovers it from a single sample's called
variants.

**No alignments.** The dataset ships raw reads and a called VCF, with no BAM, verified at
revision `59e322d2…` on 2026-09-02. Anything a called VCF cannot express — copy number,
structural rearrangement, deep-intronic change, low-level mosaicism — is invisible to this
analysis as run.

**The disease is recessive.** A single damaged copy is not an answer. The unit of evidence
is a *configuration* — two alleles considered together — not a variant.

**The population is about fifty people.** A specific variant, or a specific pattern of
affected chromosomes, is close to an identifier. That constrains what this report can say,
and Section 6 records one result withheld on exactly those grounds.

---

## 2. Method

### Annotation runs locally, by design

Variants are annotated with a **local snpEff** (`GRCh38.115`, pinned; Cingolani et al.,
2012) rather than a REST service. This is a privacy constraint, not a performance one: a
per-variant lookup against a remote annotation API would put the patient's coordinates on a
third-party server. The same rule governs every reference resource in the pipeline — each is
downloaded as a **whole public release** and matched locally, never queried per variant.

Every snpEff invocation passes `-noLog` (it otherwise reports usage statistics to its
server), `-nodownload` and `-noStats`.

### Transcript policy is declared, not defaulted

The call is made on the **MANE Select** transcript — the single transcript NCBI and EMBL-EBI
agree on per gene, and the one clinical reporting is being standardised around (Morales et
al., 2022; Wright et al., 2023). An all-transcript pass runs alongside it, and loss of
function found *only* on a non-MANE transcript is recorded as `lof_tier = non_mane` — kept
visible, never promoted to a call.

snpEff's `-canon` flag is explicitly refused: it selects the longest CDS, not MANE (SnpEff &
SnpSift documentation, n.d.), and on this panel it picks a different transcript for three of
the six genes (`results/transcript_policy/picks.tsv`, seed=42). The MANE release the
database's tags were cross-checked against is recorded in the config (`MANE.GRCh38.v1.5`),
and a script verifies the panel genes against it.

### A six-gene panel, not a genome-wide sweep

Candidates are restricted to a six-gene spindle-assembly-checkpoint panel, because the
phenotype is specific and the prior is strong. The six are not equally supported, and the
panel does not pretend they are:

| Gene | Why it is on the panel |
|---|---|
| ***BUB1B*** | Biallelic mutations established as a cause of MVA (Hanks et al., 2004), with the truncating-plus-missense architecture characterised since (Suijkerbuijk et al., 2010) |
| ***CEP57*** | Biallelic mutations reported to cause MVA (Snape et al., 2011) |
| ***TRIP13*** | Biallelic loss of function reported in MVA with Wilms tumour (Yost et al., 2017) |
| ***CEP192*** | Biallelic variants reported in MVA with tetraploidy (Guo et al., 2024) |
| ***BUB1*** | Biallelic mutations reported to cause microcephaly and developmental delay with *variable* effects on chromosome segregation (Carvalhal et al., 2022) — a checkpoint-gene disorder, not an established MVA gene |
| ***BUB3*** | **No biallelic MVA report.** Carried as a core checkpoint component on the same mechanism, so that a variant in it would be seen rather than filtered out unexamined |

This is a real limitation and I state it rather than presenting the result as genome-wide:
**a causal variant in a gene outside this panel would not be found.** For a patient whose
presentation is characteristic, the panel is the higher-precision choice; for one whose is
not, it is the wrong one.

### Evidence assembled per configuration

For each panel gene the pipeline enumerates biallelic configurations and attaches:

- **Effect class and LoF tier** per allele, on the MANE Select transcript, with exon rank and
  a last-exon flag. A truncation in the last exon escapes nonsense-mediated decay and is not
  equivalent to one that does not — the reason the ACMG/AMP PVS1 criterion is downgraded for
  such variants rather than applied flat (Abou Tayoun et al., 2018). Two flags are recorded
  separately and neither is used to filter: `nmd_escape_last_exon` from the exon rank, and
  `snpeff_predicts_no_nmd` from the annotator's own prediction.
- **ClinVar cross-reference** — submitted interpretations and **review status**, from the
  whole GRCh38 release matched locally (Landrum et al., 2018; National Center for
  Biotechnology Information, 2026). The pipeline reports what ClinVar holds; it classifies
  nothing itself.
- **gnomAD v4.1.1 frequencies** — allele count, number, frequency, homozygote count and
  highest-group frequency, read over the whole panel-gene span (Chen et al., 2024). Only CC0
  count and frequency fields are read; the SpliceAI scores shipped in the same files are
  **CC BY-NC and are not read** (Genome Aggregation Database, 2026), which a test enforces.
- **Phase**, where the VCF carries it, and an explicit `unphased_pair` flag where it does not.

---

## 3. The finding

Of the six panel genes, **exactly one carries any PASS loss-of-function allele**: *BUB1B*.
The other five yield no biallelic configuration at all.

*Source: `results/l0_genomics/causal_gene_call.json` and `variant_calls.json` (seed 42),
produced by `python -m src.pipeline --only l0_genomics` on the config recorded in
`results/_manifest.json`.*

The *BUB1B* configuration is classed `lof_plus_protein_altering`:

| Allele | Evidence |
|---|---|
| **Truncating** | Loss of function on MANE Select, **outside the last exon** so nonsense-mediated decay applies. ClinVar holds it as **pathogenic at two-star review status** for mosaic variegated aneuploidy syndrome 1 |
| **Missense** | In the **protein kinase domain** — residues 766–1050 of the 1,050-residue protein, matched to the MANE Select translation's length (The UniProt Consortium, 2026). **No ClinVar record.** gnomAD reports it extremely rare and **never homozygous**, at a position both gnomAD data sets sampled |

### Why this configuration and not another

**The pattern is the reported one.** *BUB1B* was established as an MVA gene through
truncating and missense mutations (Hanks et al., 2004). In patients with biallelic mutations,
"a missense mutation pairs with a truncating mutation," with the missense change
"consistently" in or near the BUBR1 kinase domain (Suijkerbuijk et al., 2010). Both
quotations are from the abstracts; the full texts were not read, and I mark that rather than
implying otherwise.

So the configuration is not merely *a* biallelic pair — it is the specific architecture the
literature reports for this gene, in the only panel gene carrying a qualifying allele, with
the truncating half independently classified pathogenic for this exact syndrome.

---

## 4. The submission

Two rows, ranked by my own estimated probability of causal relationship.

| Rank | Row | EPCR | Why |
|---|---|---|---|
| 1 | The compound-heterozygous **pair** | 0.70 | The best-supported configuration. Not higher, because the pair is unphased and the missense allele is untested |
| 2 | The **truncating allele alone** | 0.45 | Hedge against the second hit being something a called VCF cannot surface |

**The second row is scoring-neutral and I know it.** I checked the submission against the
organizers' published `evaluation.py` — the parser and scorer only, never the answer key —
and found that partial credit is already awarded from the pair row, because the scorer tests
for intersection with the true set. The row is kept because its note states a real failure
mode for the human reviewer, not because it helps the metric.

I also measured the alternative ordering. Leading with the truncating allele alone would
raise F-max from 0.500 to 0.667 *if* the second allele is wrong, and cost 50 rank points if
the pair is right. Rank points run 0–100 and F-max 0–1; that trade is bad and I did not take
it.

Chromosome names are `chr`-prefixed to match the template — the dataset's VCF uses bare
contig names, and the builder converts and asserts it, because an unprefixed chromosome
scores zero while looking entirely correct in the file.

**The limitations in Section 8 are not confined to this report.** Three of the four travel
with the call into the submission's own `notes` field, so a reader who never opens this
document still gets them.

---

## 5. Scientific rigor

### Evidence discipline

Every claim in this project is one of two admissible kinds: a **pipeline data result**
(artifact path + seed + config, regenerable) or a **cited source** in APA 7th with a
resolvable DOI. Every DOI was resolved against Crossref and checked for retraction and
correction notices; the verification log is in `docs/references.md` with the date each check
ran. An LLM-produced citation is treated as a hypothesis about the literature until
resolved — and dropped if it does not resolve.

The pipeline reports what its sources say and classifies nothing itself. ClinVar's verdicts
are ClinVar's; gnomAD's frequencies are gnomAD's; snpEff's consequences are snpEff's. Where
those disagree with each other, the disagreement is recorded rather than resolved silently.

### Reproducibility

Re-running a layer on unchanged inputs reproduces it **byte for byte**; this is verified in
the repository for the downstream layers, and every dependency, data release and seed is
pinned.

| Pinned | How |
|---|---|
| Libraries | `environment.lock.yml` (exact linux-64 solve) and `environment.yml` (portable, fully pinned) |
| ClinVar | a **dated** weekly release, not the rolling `clinvar.vcf.gz`, with an archive fallback |
| snpEff database, MANE, gnomAD | version or release date in `config/pipeline.yaml` |
| Every cached download | SHA-256 recorded beside the result |

`results/_manifest.json` records the seed, the SHA-256 of the config file, a digest of the
resolved settings with machine-specific paths excluded, and the versions of the libraries
that decide the numbers. ClinVar is pinned to a **dated** weekly release because the rolling
file is replaced weekly: two runs a fortnight apart would cross-reference the patient's
alleles against different archives, and no output would say so.

---

## 6. Innovation: aneuploidy burden without alignments

MVA's phenotype *is* mosaic aneuploidy, so quantifying it should corroborate a call. **The
dataset ships no BAM** — raw reads plus a called VCF — so depth-based methods have no input.

I estimate burden from **per-chromosome B-allele frequency in the VCF itself**: mean
|BAF − 0.5| at heterozygous biallelic SNVs, in 10 Mb windows, at minimum depth 10, against a
baseline of the sample's own autosomal windows (median + scaled MAD). An allele ratio at a
single locus is self-normalising, so it needs none of the GC-bias and mappability correction
that depth would (Conlin et al., 2010; Loh et al., 2018).

The raw statistic sits at a noise floor when the true shift is zero and responds
quadratically to small shifts, so the mosaic-fraction estimator inverts the folded-normal
mean for the true BAF shift before converting — without that step, small fractions are
systematically overstated.

**Resolving power on this sample:** smallest detectable mosaic fraction **≈ 0.098**, with a
conservative bound of **≈ 0.245**, at z = 5, from the sample's own median depth (44×) and
window-to-window scatter rather than from a configured constant.

*Source: `results/l0_genomics/aneuploidy_burden.json`, `sensitivity` (seed 42).* The two
figures are quoted as a pair on purpose: the artifact records the estimate as
**model-dependent** — the optimistic figure assumes the baseline excess and a mosaic shift
add in shift space, and a spike-in simulation at this depth put the crossing nearer 0.25.
Quoting only the lower number would report the most favourable model as if it were the
measurement.

### What I deliberately do not show

**The per-chromosome result is withheld** — in every form, including aggregate counts. MVA
case series publish karyotypes; in a population of roughly fifty, which chromosomes are
involved and at what fraction is close to a fingerprint. A count of affected chromosomes is
not a safe aggregate either: in a population that small it narrows the field by itself.

I publish the **method and its resolving power**, because that is what the innovation claim
is about — the claim is *"allele ratios alone resolve mosaic aneuploidy without alignments,
down to ~10%"* — and not this patient's result from it. The withholding is **named, not
silent**: omitting it quietly would be indistinguishable from never having built it
(repository decision D20).

---

## 7. Governance, ethics and disclosure

This work involves a real child's genome. Several constraints were fixed before any code was
written and none was relaxed:

- **Patient data never enters the repository.** All data lives under gitignored directories,
  and reference caches sit outside the working tree entirely — the cache module *refuses* a
  path inside the repository.
- **No per-variant remote query, anywhere in the pipeline.** Every reference resource is
  downloaded as a whole public release and matched locally, so no coordinate of this
  patient's ever reaches a third-party server.
- **Clinical phenotype is treated as patient data.** HPO terms are parsed from the data
  directory at run time and never written into config, source, tests or any committed file —
  a specific combination of features is identifying in a population this size.
- **No recontact** with the patient, the family or MVA Society, and nothing published beyond
  what the family already shares publicly through their own blog posts.
- **The submission CSV is never committed.** It is written under `results/`, which is
  gitignored and inside the project's deletion scope, and uploaded from there.
- **Deletion of all held data is committed**, with a register of every custody location —
  including one row for the copy this submission itself creates, which the organizers hold
  and this project cannot purge.

**Generative AI disclosure.** The pipeline's own reasoning step (used in Track 2, not in this
call) runs a **local open-weights model** — `qwen2.5-7b-instruct-q4_k_m` — served over
loopback, with a guard that refuses any non-loopback endpoint. Nothing patient-derived leaves
the machine. Development assistance used a commercially available assistant under a consumer
subscription with model training on inputs and outputs disabled. Both are recorded in
decision D17 and in `COMPLIANCE.md`.

**No AI system made the causal-gene call.** L0 assembled the evidence; a person reviewed it
and decided, and the decision is recorded with its date, its evidence and its dissenting
points (repository decision D9).

---

## 8. Limitations

Stated plainly, because a variant call for a patient should arrive with its weaknesses
attached. They are the reason the primary row's EPCR is 0.70 rather than 0.95, and three of
the four travel with the call into the submission's own `notes` field.

1. **Phase is unknown — the largest uncertainty.** No parental samples, no recontact
   permitted. If both variants sit on the same copy, the other copy is intact and the
   configuration is not causal. Nothing available to me can settle it. A confident EPCR here
   would be a statement about my data that my data does not support.
2. **The missense allele is untested.** The literature supports the *pattern*; no functional
   study shows that *this* change impairs BUBR1. Rarity is necessary, not sufficient — under
   the ACMG/AMP framework low frequency is supporting evidence, never on its own a
   classification (Richards et al., 2015) — and absence from ClinVar is absence of a
   submission, not evidence of benignity.
3. **Unseen second hits.** Copy-number, structural, deep-intronic and low-level mosaic
   variants are **not assessable from a called VCF**. The missense allele is the best second
   candidate *in the called set*, which is a narrower claim than the best that exists. The
   submission's second row hedges exactly this case.
4. **No detectable constitutional aneuploidy.** The pipeline estimates per-chromosome
   aneuploidy burden from B-allele frequency (Section 6) and every chromosome came back at
   baseline or below the detection limit. This does **not** exclude low-level mosaicism
   beneath the method's floor, but it is honestly *not* supporting evidence, and I decline to
   present it as such. This is the one caveat that does not reach the `notes` field, because
   it is an absence rather than a property of the call.
5. **A six-gene panel, not a genome-wide sweep** (Section 2). A causal variant outside the
   panel would not be found.
6. **The causal gene is a research premise, not a diagnosis.** `BUB1B` was fixed by a person
   at a recorded gate from ClinVar and gnomAD evidence; the allele pair is unphased and those
   uncertainties travel into every downstream output, including all of Track 2.

---

## 9. What I would do next

**Align the supplied reads and attempt read-backed phasing.** The dataset ships raw reads;
the pipeline does not use them. If the two variants fall within the span of a single fragment
or read pair, phase can be read directly off the alignments — which would close, or overturn,
the largest uncertainty in this call. If they do not, the attempt returns nothing, and that
is still worth knowing rather than assuming. This is run time plus alignment, not new method.

**Call copy-number and structural variants from those same alignments.** Limitation 3 exists
because a called VCF cannot express them. It is the only route to the second allele if the
missense is not it.

**Functional characterisation of the missense allele** is what would settle limitation 2, and
it is **out of scope for this project** — it is wet-lab work, and this is a computational
nomination. Naming it is not proposing it.

---

## Reproducing this

```bash
conda env create -f environment.lock.yml     # exact linux-64 solve
conda activate mva-track2
snpEff download -noLog GRCh38.115
mkdir -p data results
# download the gated dataset into ./data — see DATA.md

PYTHONHASHSEED=42 python -m src.pipeline --only l0_genomics
python scripts/build_track1_submission.py    # -> results/submissions/track1/ (gitignored)
python scripts/package_submissions.py        # adds the report and an upload checklist
```

`scripts/transcript_policy_check.py` verifies the panel's MANE Select transcripts against the
recorded release. `mngmt/decisions.md` records each decision, what forced it, and what it
obliged downstream — including the ones that went against me.

---

## References

Every entry below is cited in the text above, and every citation in the text above appears
below. DOIs were resolved against Crossref and checked for retraction and correction
notices; the log with the date each check ran is in [`docs/references.md`](references.md).
Machine-readable source: [`docs/references.bib`](references.bib).

Abou Tayoun, A. N., Pesaran, T., DiStefano, M. T., Oza, A., Rehm, H. L., Biesecker, L. G., &
Harrison, S. M. (2018). Recommendations for interpreting the loss of function PVS1 ACMG/AMP
variant criterion. *Human Mutation, 39*(11), 1517–1524. https://doi.org/10.1002/humu.23626

Carvalhal, S., Bader, I., Rooimans, M. A., Oostra, A. B., Balk, J. A., Feichtinger, R. G.,
Beichler, C., Speicher, M. R., van Hagen, J. M., Waisfisz, Q., van Haelst, M., Bruijn, M.,
Tavares, A., Mayr, J. A., Wolthuis, R. M. F., Oliveira, R. A., & de Lange, J. (2022).
Biallelic BUB1 mutations cause microcephaly, developmental delay, and variable effects on
cohesion and chromosome segregation. *Science Advances, 8*(3), Article eabk0114.
https://doi.org/10.1126/sciadv.abk0114

Chen, S., Francioli, L. C., Goodrich, J. K., Collins, R. L., Kanai, M., Wang, Q., Alföldi, J.,
Watts, N. A., Vittal, C., Gauthier, L. D., Poterba, T., Wilson, M. W., Tarasova, Y., Phu, W.,
Grant, R., Yohannes, M. T., Koenig, Z., Farjoun, Y., Banks, E., … Karczewski, K. J. (2024). A
genomic mutational constraint map using variation in 76,156 human genomes. *Nature,
625*(7993), 92–100. https://doi.org/10.1038/s41586-023-06045-0

Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., Land, S. J., Lu, X., &
Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide
polymorphisms, SnpEff: SNPs in the genome of *Drosophila melanogaster* strain w1118; iso-2;
iso-3. *Fly, 6*(2), 80–92. https://doi.org/10.4161/fly.19695

Conlin, L. K., Thiel, B. D., Bonnemann, C. G., Medne, L., Ernst, L. M., Zackai, E. H.,
Deardorff, M. A., Krantz, I. D., Hakonarson, H., & Spinner, N. B. (2010). Mechanisms of
mosaicism, chimerism and uniparental disomy identified by single nucleotide polymorphism array
analysis. *Human Molecular Genetics, 19*(7), 1263–1275. https://doi.org/10.1093/hmg/ddq003

Genome Aggregation Database. (2026). *Policies*. Retrieved September 17, 2026, from
https://gnomad.broadinstitute.org/policies

Guo, J., He, W.-B., Dai, L., Tian, F., Luo, Z., Shen, F., Tu, M., Zheng, Y., Zhao, L., Tan, C.,
Guo, Y., Meng, L.-L., Liu, W., Deng, M., Wu, X., Peng, Y., Zhang, S., Lu, G.-X., Lin, G., …
Yang, Y. (2024). Mosaic variegated aneuploidy syndrome with tetraploid, and predisposition to
male infertility triggered by mutant CEP192. *Human Genetics and Genomics Advances, 5*(1),
Article 100256. https://doi.org/10.1016/j.xhgg.2023.100256

Hanks, S., Coleman, K., Reid, S., Plaja, A., Firth, H., FitzPatrick, D., Kidd, A., Méhes, K.,
Nash, R., Robin, N., Shannon, N., Tolmie, J., Swansbury, J., Irrthum, A., Douglas, J., &
Rahman, N. (2004). Constitutional aneuploidy and cancer predisposition caused by biallelic
mutations in BUB1B. *Nature Genetics, 36*(11), 1159–1161. https://doi.org/10.1038/ng1449

Landrum, M. J., Lee, J. M., Benson, M., Brown, G. R., Chao, C., Chitipiralla, S., Gu, B.,
Hart, J., Hoffman, D., Jang, W., Karapetyan, K., Katz, K., Liu, C., Maddipatla, Z., Malheiro,
A., McDaniel, K., Ovetsky, M., Riley, G., Zhou, G., … Maglott, D. R. (2018). ClinVar: Improving
access to variant interpretations and supporting evidence. *Nucleic Acids Research, 46*(D1),
D1062–D1067. https://doi.org/10.1093/nar/gkx1153

Loh, P.-R., Genovese, G., Handsaker, R. E., Finucane, H. K., Reshef, Y. A., Palamara, P. F.,
Birmann, B. M., Talkowski, M. E., Bakhoum, S. F., McCarroll, S. A., & Price, A. L. (2018).
Insights into clonal haematopoiesis from 8,342 mosaic chromosomal alterations. *Nature,
559*(7714), 350–355. https://doi.org/10.1038/s41586-018-0321-x

Morales, J., Pujar, S., Loveland, J. E., Astashyn, A., Bennett, R., Berry, A., Cox, E.,
Davidson, C., Ermolaeva, O., Farrell, C. M., Fatima, R., Gil, L., Goldfarb, T., Gonzalez,
J. M., Haddad, D., Hardy, M., Hunt, T., Jackson, J., Joardar, V. S., … Murphy, T. D. (2022). A
joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature, 604*(7905),
310–315. https://doi.org/10.1038/s41586-022-04558-8

National Center for Biotechnology Information. (2026). *ClinVar review status*. Retrieved
September 16, 2026, from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde,
M., Lyon, E., Spector, E., Voelkerding, K., & Rehm, H. L. (2015). Standards and guidelines for
the interpretation of sequence variants: A joint consensus recommendation of the American
College of Medical Genetics and Genomics and the Association for Molecular Pathology.
*Genetics in Medicine, 17*(5), 405–424. https://doi.org/10.1038/gim.2015.30

Scott, R. H., Stiller, C. A., Walker, L., & Rahman, N. (2006). Syndromes and constitutional
chromosomal abnormalities associated with Wilms tumour. *Journal of Medical Genetics, 43*(9),
705–715. https://doi.org/10.1136/jmg.2006.041723

Snape, K., Hanks, S., Ruark, E., Barros-Núñez, P., Elliott, A., Murray, A., Lane, A. H.,
Shannon, N., Callier, P., Chitayat, D., Clayton-Smith, J., FitzPatrick, D. R., Gisselsson, D.,
Jacquemont, S., Asakura-Hay, K., Micale, M. A., Tolmie, J., Turnpenny, P. D., Wright, M., …
Rahman, N. (2011). Mutations in CEP57 cause mosaic variegated aneuploidy syndrome. *Nature
Genetics, 43*(6), 527–529. https://doi.org/10.1038/ng.822

SnpEff & SnpSift documentation. (n.d.). *Commands & command line options*. Retrieved September
11, 2026, from https://pcingola.github.io/SnpEff/snpeff/commandline/

Suijkerbuijk, S. J. E., van Osch, M. H. J., Bos, F. L., Hanks, S., Rahman, N., & Kops, G. J. P.
L. (2010). Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome
mosaic variegated aneuploidy. *Cancer Research, 70*(12), 4891–4900.
https://doi.org/10.1158/0008-5472.can-09-4319

The UniProt Consortium. (2026). *O60566 · BUB1B_HUMAN: Mitotic checkpoint
serine/threonine-protein kinase BUB1 beta* (UniProtKB/Swiss-Prot release 2026_03) [Data set].
Retrieved September 17, 2026, from https://rest.uniprot.org/uniprotkb/O60566.json

Wright, C. F., FitzPatrick, D. R., Ware, J. S., Rehm, H. L., & Firth, H. V. (2023). Importance
of adopting standardized MANE transcripts in clinical reporting. *Genetics in Medicine,
25*(2), Article 100331. https://doi.org/10.1016/j.gim.2022.10.013

Yost, S., de Wolf, B., Hanks, S., Zachariou, A., Marcozzi, C., Clarke, M., de Voer, R. M.,
Etemad, B., Uijttewaal, E., Ramsay, E., Wylie, H., Elliott, A., Picton, S., Smith, A.,
Smithson, S., Seal, S., Ruark, E., Houge, G., Pines, J., … Rahman, N. (2017). Biallelic TRIP13
mutations predispose to Wilms tumor and chromosome missegregation. *Nature Genetics, 49*(7),
1148–1151. https://doi.org/10.1038/ng.3883

---

*This is a computational prediction, not a diagnosis and not a classification of either
allele.*
