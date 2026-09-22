# Identifying the causal variant pair in a proband with mosaic variegated aneuploidy

**Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1 (Variant Prediction)**

Seed 42 · GRCh38 · snpEff GRCh38.115 · ClinVar release 2026-09-13 · gnomAD v4.1.1

> **A note on what this document does not contain.** No variant coordinate, HGVS expression
> or ClinVar accession appears here. The submission CSV carries the identifiers, as the
> format requires; this report describes the method and the evidence at category level. In a
> disease with roughly fifty patients worldwide, a specific variant can identify a child, and
> the readable document is the one most likely to be forwarded. The reasoning is recorded as
> decision D21 in the repository.

---

## Summary

We submit a **compound-heterozygous pair in *BUB1B*** as the primary finding: one truncating
allele that ClinVar holds as pathogenic at two-star review status for mosaic variegated
aneuploidy syndrome 1, paired with a rare missense variant in the BUBR1 kinase domain.

The call was produced by the L0 layer of a pipeline built for Track 2, where the causal gene
is the premise the drug search is organised around. Identifying it was not a side quest — it
is the gate the rest of the work depends on, and it was decided by a person reviewing the
evidence, not by the pipeline setting its own premise.

**The pair is unphased and we say so in the submission itself.** With no parental samples
and no recontact permitted, a *cis* arrangement cannot be excluded. That is the largest
uncertainty in this call and no data available to us closes it.

---

## 1. Method

### Annotation runs locally, by design

Variants are annotated with a **local snpEff** (`GRCh38.115`, pinned) rather than a REST
service. This is a privacy constraint, not a performance one: a per-variant lookup against
a remote annotation API would put this child's coordinates on a third-party server. The same
rule governs every reference resource in the pipeline — each is downloaded as a **whole
public release** and matched locally, never queried per variant.

Every snpEff invocation passes `-noLog` (it otherwise reports usage statistics to its
server), `-nodownload` and `-noStats`.

### Transcript policy is declared, not defaulted

The call is made on the **MANE Select** transcript. An all-transcript pass runs alongside it,
and loss of function found *only* on a non-MANE transcript is recorded as `lof_tier =
non_mane` — kept visible, never promoted to a call.

snpEff's `-canon` flag is explicitly refused: it selects the longest CDS, which is not MANE
and differs from it on several genes in this panel. The MANE release the database's tags were
cross-checked against is recorded in the config (`MANE.GRCh38.v1.5`), and a script verifies
the panel genes against it.

### A six-gene panel, not a genome-wide sweep

Candidates are restricted to the established MVA genes — ***BUB1B*, *CEP57*, *TRIP13*,
*BUB1*, *BUB3*, *CEP192*** — because the phenotype is specific and the prior is strong
(Hanks et al., 2004; Suijkerbuijk et al., 2010).

This is a real limitation and we state it rather than presenting the result as genome-wide:
**a causal variant in a gene outside this panel would not be found.** For a proband whose
presentation is characteristic, the panel is the higher-precision choice; for one whose is
not, it is the wrong one.

### Evidence assembled per configuration

For each panel gene the pipeline enumerates biallelic configurations and attaches:

- **Effect class and LoF tier** per allele, on the MANE Select transcript, with exon rank and
  a last-exon flag (a truncation in the last exon escapes nonsense-mediated decay and is not
  equivalent).
- **ClinVar cross-reference** — submitted interpretations and **review status**, from the
  whole GRCh38 release matched locally. The pipeline reports what ClinVar holds; it
  classifies nothing itself.
- **gnomAD v4.1.1 frequencies** — allele count, number, frequency, homozygote count and
  highest-group frequency, read over the whole panel-gene span. Only CC0 count and frequency
  fields are read; the SpliceAI scores shipped in the same files are **CC BY-NC and are not
  read**, which a test enforces.
- **Phase**, where the VCF carries it, and an explicit `unphased_pair` flag where it does not.

---

## 2. The finding

Of the six panel genes, **exactly one carries any PASS loss-of-function allele**: *BUB1B*.
The other five yield no biallelic configuration at all.

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
quotations are from the abstracts; the full texts were not read, and we mark that rather than
implying otherwise.

So the configuration is not merely *a* biallelic pair — it is the specific architecture the
literature reports for this gene, in the only panel gene carrying a qualifying allele, with
the truncating half independently classified pathogenic for this exact syndrome.

---

## 3. What we are not claiming

Four caveats travel with this call into the submission's own `notes` field. They are the
reason the primary row's EPCR is 0.70 rather than 0.95.

**1. Phase is unknown — the largest uncertainty.** No parental samples, no recontact
permitted. If both variants sit on the same copy, the other copy is intact and the
configuration is not causal. Nothing available to us can settle it. A confident EPCR here
would be a statement about our data that our data does not support.

**2. The missense allele is untested.** The literature supports the *pattern*; no functional
study shows that *this* change impairs BUBR1. Rarity is necessary, not sufficient, and
absence from ClinVar is absence of a submission, not evidence of benignity.

**3. Unseen second hits.** Copy-number, structural, deep-intronic and low-level mosaic
variants are **not assessable from a called VCF**. The missense allele is the best second
candidate *in the called set*, which is a narrower claim than the best that exists. The
submission's second row hedges exactly this case.

**4. No detectable constitutional aneuploidy.** The pipeline estimates per-chromosome
aneuploidy burden from B-allele frequency (see below) and every chromosome came back at
baseline or below the detection limit. This does **not** exclude low-level mosaicism beneath
the method's floor, but it is honestly *not* supporting evidence, and we decline to present
it as such.

---

## 4. Aneuploidy burden without alignments

MVA's phenotype *is* mosaic aneuploidy, so quantifying it should corroborate a call. **The
dataset ships no BAM** — raw reads plus a called VCF — so depth-based methods have no input.

We estimate burden from **per-chromosome B-allele frequency in the VCF itself**: mean
|BAF − 0.5| at heterozygous biallelic SNVs, in 10 Mb windows, at minimum depth 10, against a
baseline of the sample's own autosomal windows (median + scaled MAD). An allele ratio at a
single locus is self-normalising, so it needs none of the GC-bias and mappability correction
that depth would (Conlin et al., 2010; Loh et al., 2018).

The raw statistic sits at a noise floor when the true shift is zero and responds
quadratically to small shifts, so the mosaic-fraction estimator inverts the folded-normal
mean for the true BAF shift before converting — without that step, small fractions are
systematically overstated.

**Resolving power on this sample:** smallest detectable mosaic fraction **≈ 0.098**, with a
conservative bound of **≈ 0.245**, at z = 5.

**The per-chromosome result is withheld** — in every form, including aggregate counts. MVA
case series publish karyotypes; in a population of roughly fifty, which chromosomes are
involved and at what fraction is close to a fingerprint. We publish the method and what it
can resolve, because that is what the innovation claim is about. The withholding is named
rather than silent (repository decision D20).

---

## 5. The submission

Two rows, ranked by our own estimated probability of causal relationship.

| Rank | Row | EPCR | Why |
|---|---|---|---|
| 1 | The compound-heterozygous **pair** | 0.70 | The best-supported configuration. Not higher, because the pair is unphased and the missense allele is untested |
| 2 | The **truncating allele alone** | 0.45 | Hedge against the second hit being something a called VCF cannot surface |

**The second row is scoring-neutral and we know it.** We checked the submission against the
organizers' published `evaluation.py` — the parser and scorer only, never the answer key —
and found that partial credit is already awarded from the pair row, because the scorer tests
for intersection with the true set. The row is kept because its note states a real failure
mode for the human reviewer, not because it helps the metric.

We also measured the alternative ordering. Leading with the truncating allele alone would
raise F-max from 0.500 to 0.667 *if* the second allele is wrong, and cost 50 rank points if
the pair is right. Rank points run 0–100 and F-max 0–1; that trade is bad and we did not take
it.

Chromosome names are `chr`-prefixed to match the template — the dataset's VCF uses bare
contig names, and the builder converts and asserts it, because an unprefixed chromosome
scores zero while looking entirely correct in the file.

---

## 6. Reproducibility

```bash
conda env create -f environment.lock.yml     # exact linux-64 solve
conda activate mva-track2
snpEff download -noLog GRCh38.115
PYTHONHASHSEED=42 python -m src.pipeline --only l0_genomics
python scripts/build_track1_submission.py    # -> results/submissions/track1/ (gitignored)
python scripts/package_submissions.py        # adds the report and an upload checklist
```

Re-running a layer on unchanged inputs reproduces it byte for byte; this is verified in the
repository for the downstream layers, and every dependency, data release and seed is pinned.
`results/_manifest.json` records the seed, the SHA-256 of the config, and the versions of the
libraries that decide the numbers. ClinVar is pinned to a **dated** weekly release rather than
the rolling `clinvar.vcf.gz`, whose contents change weekly.

**The submission CSV is never committed.** It is written under `results/`, which is gitignored
and inside the project's deletion scope, and uploaded from there.

---

## 7. Generative AI disclosure

The pipeline's own reasoning step (used in Track 2, not in this call) runs a **local
open-weights model** — `qwen2.5-7b-instruct-q4_k_m` — served over loopback, with a guard that
refuses any non-loopback endpoint. Nothing patient-derived leaves the machine.

**Development assistance: Anthropic, Claude, consumer subscription (Pro), with model training
on inputs and outputs disabled in account settings.**

No AI system made the causal-gene call. L0 assembled the evidence; a person reviewed it and
decided, and the decision is recorded with its date, its evidence and its dissenting points
(repository decision D9).

---

## 8. Ethics and data handling

- Patient data never enters the repository. All data lives under gitignored directories, and
  reference caches sit outside the working tree entirely — the cache module *refuses* a path
  inside the repository.
- No per-variant remote query, anywhere in the pipeline.
- Clinical phenotype is treated as patient data: HPO terms are parsed at run time and never
  written into config, source, tests or any committed file.
- No recontact with the subject, family or MVA Society.
- Deletion of all held data is committed with a register of every custody location, including
  one row for the copy this submission itself creates — which the organizers hold and this
  project cannot purge.

---

## References

Conlin, L. K., Thiel, B. D., Bonnemann, C. G., Medne, L., Ernst, L. M., Zackai, E. H.,
Deardorff, M. A., Krantz, I. D., Hakonarson, H., & Spinner, N. B. (2010). Mechanisms of
mosaicism, chimerism and uniparental disomy identified by single nucleotide polymorphism array
analysis. *Human Molecular Genetics, 19*(7), 1263–1275. https://doi.org/10.1093/hmg/ddq003

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

Suijkerbuijk, S. J. E., van Osch, M. H. J., Bos, F. L., Hanks, S., Rahman, N., & Kops, G. J. P.
L. (2010). Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome
mosaic variegated aneuploidy. *Cancer Research, 70*(12), 4891–4900.
https://doi.org/10.1158/0008-5472.can-09-4319

Full reference list with the Crossref verification log, including retraction and correction
checks: [`docs/references.md`](references.md).

---

*This is a computational prediction, not a diagnosis and not a classification of either
allele.*
