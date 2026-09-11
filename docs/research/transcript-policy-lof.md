# Which Transcripts Count for a Loss-of-Function Call

> **Adopted as decision D5 on 2026-09-11** — see [`../../mngmt/decisions.md`](../../mngmt/decisions.md).
> Reviewed and citation-checked 2026-09-11. This document is the evidence behind the decision; the
> policy itself lives in [`../../config/pipeline.yaml`](../../config/pipeline.yaml)
> (`annotator.transcript_policy`) and
> [`../../src/l0_genomics/transcripts.py`](../../src/l0_genomics/transcripts.py).

## Status and how this may be used

**Citations: verified.** Every journal source was resolved against Crossref and retraction-checked
on 2026-09-11 (`scripts/verify_refs.py`). None is retracted. Two carry notices that do not touch
the results used here — Cummings et al. (2020) an Author Correction adding a consortium member,
Karczewski et al. (2020) the same plus an Addendum extending its disease-gene analysis. Keys and
the verification log are in [`../references.bib`](../references.bib) and
[`../references.md`](../references.md).

**Data results: reproducible, and patient-free.** Every number under
[What this repository's data show](#what-this-repositorys-data-show) is an evidence-type-A result
produced by [`../../scripts/transcript_policy_check.py`](../../scripts/transcript_policy_check.py)
into `results/transcript_policy/` (`picks.tsv`, `clinvar_summary.tsv`, `clinvar_variants.tsv`,
`territory.tsv`, `manifest.json`; `seed=42`, `config/pipeline.yaml`; the script uses no RNG). Its
inputs are public — the MANE v1.5 summary table, ClinVar's GRCh38 VCF (fileDate 2026-09-05), the
snpEff GRCh38.115 database — plus invented substitutions. **The proband's variants were not
annotated for this review.**

**Scope.** Which transcripts L0 should consult when it classifies a variant as loss-of-function
(LoF). Out of scope: missense interpretation, the causal-gene call itself, and any clinical
classification — this pipeline generates hypotheses and makes no clinical variant calls.

---

## TL;DR

- **Recommended: MANE Select as the primary transcript, with an all-transcript pass kept as
  flagged, lower-tier evidence.** Never snpEff's `-canon`.
- The literature is consistent: transcript choice changes LoF calls substantially (McCarthy et al.,
  2014); LoF on a subset of a gene's transcripts is common and enriched for false positives
  (MacArthur et al., 2012; Singer-Berk et al., 2023; Cummings et al., 2020); clinical LoF criteria
  are assessed on biologically relevant transcripts (Abou Tayoun et al., 2018); and MANE Select is
  the transcript set built to be that default (Morales et al., 2022), capturing all but a sliver of
  known pathogenic variants (Pozo et al., 2022).
- On this panel, **snpEff `-canon` is not MANE for BUB1, BUB1B and CEP192**, and in BUB1B it would
  leave 20 bp of MANE Select coding/splice sequence unable to receive a LoF call. `-tag MANE_Select`
  reproduces NCBI's MANE Select transcript for all six genes.
- Known pathogenic variants cannot choose between policies — all 118 public ClinVar P/LP variants
  in the panel get the same call under all three. The choice matters for **novel** variants, and
  the panel has 448 bp where only a non-MANE transcript would permit a LoF call.

---

## The question

L0 must identify biallelic LoF in spindle-assembly-checkpoint genes — MVA is caused by biallelic
mutation (Hanks et al., 2004). snpEff annotates every transcript in the database, so a single
variant receives a different consequence per transcript: the install smoke test gave one invented
BUB1B variant a 5′-UTR consequence on one transcript and an upstream consequence on two others.
Whether that variant "is LoF" therefore depends on which transcripts are allowed to vote.

| Option | snpEff | Rule |
|---|---|---|
| **A. All transcripts** | default | LoF if LoF on *any* transcript |
| **B. snpEff canonical** | `-canon` | LoF if LoF on snpEff's canonical transcript |
| **C. MANE Select** | `-tag MANE_Select` | LoF if LoF on the MANE Select transcript |
| **D. Tiered** | C, then A | C is the call; LoF found *only* under A is kept, flagged, one tier down |

---

## What the literature says

**1. The choice of transcripts changes LoF calls more than any other annotation decision.**
Annotating the same whole genomes against two transcript sets, McCarthy et al. (2014) "found only
44% agreement in annotations for putative loss-of-function variants when using the RefSeq and
Ensembl transcript sets"; with the transcript set fixed, two annotators still agreed on only 65%
of LoF variants. Software matters too — which is why this review measures snpEff itself rather than
borrowing results obtained with VEP.

**2. LoF confined to some of a gene's transcripts is common, and it is where false positives
concentrate.** In 185 genomes, "415 (32.3%) of our high-confidence LoF variants are partial LoF
variants, affecting only a subset of the known transcripts from the affected gene" (MacArthur et
al., 2012), and annotation/reference errors accounted for 26.8% of examined candidates. gnomAD's
LOFTEE removes "terminal truncation variants, as well as rescued splice variants, that are predicted
to escape nonsense-mediated decay", and of 443,769 high-confidence pLoF variants, 413,097 "fall on
the canonical transcripts of 16,694 genes" (Karczewski et al., 2020). Applied to 22 recessive
disease genes, a curation framework found predicted LoF evasion or artifacts in 27.3% of
high-confidence pLoF variants — the leading reasons being last-exon location, homopolymers, low
proportion-expressed-across-transcripts (pext) scores, and cryptic in-frame splice rescue
(Singer-Berk et al., 2023). Expression is the underlying issue: de novo pLoF variants in regions
with little evidence of expression "are as equally distributed in cases versus controls as de novo
synonymous variants" (Cummings et al., 2020).

*Evidence grade: population-genetic and curation studies. They bear on annotation accuracy, not on
SAC-gene biology.*

**3. Clinical LoF criteria are applied on biologically relevant transcripts.** The ACMG/AMP
framework's very-strong LoF criterion, PVS1 (Richards et al., 2015), was refined by ClinGen: "PVS1
at any strength level should not be applied if the putative loss of function variant affects
exon(s) which is/are missing from alternate biologically relevant transcript(s)", and a premature
stop's impact "depends on the location of the new termination codon within the most biologically
relevant transcript(s)" (Abou Tayoun et al., 2018). Relevance is judged "based on functional and/or
expression evidence", and the presence of pathogenic variants in an exon supports it.

*Evidence grade: expert consensus recommendations.*

**4. MANE Select is the transcript set designed to be the default.** "The MANE Select set
identifies a representative transcript for each human protein-coding gene, whereas the MANE Plus
Clinical set provides additional transcripts at loci where the Select transcripts alone are not
sufficient to report all currently known clinical variants" — 55 genes at publication (Morales et
al., 2022). Ensembl's own canonical transcript *is* the MANE Select transcript wherever curators at
both groups have approved it (Ensembl, 2026), and Wright et al. (2023) argue for adopting MANE in
clinical reporting. On coverage: "Prior to manual curation, all but 67 of the 33,736 pathogenic
variants with PubMed support mapped to MANE Select transcripts rather than alternative
transcripts", and "very few 'Pathogenic' variants mapped to alternative exons (0.37%)" (Pozo et
al., 2022).

**5. "Canonical" in snpEff means something else.** snpEff defines canonical transcripts as "the
longest CDS of amongst the protein coding transcripts in a gene", and warns there is "no warranties
that what SnpEff considers a canonical transcript will match exactly what UCSC or ENSEMBL consider
a canonical transcript" (SnpEff & SnpSift documentation, n.d.; Cingolani et al., 2012). Pozo et al.
(2022) note that "although the longest transcript has traditionally been chosen as the reference,
APPRIS principal and MANE Select transcripts, biologically supported reference sequences, are now
available". Length is a heuristic; MANE is curated.

---

## What this repository's data show

All from `results/transcript_policy/` (seed=42, `config/pipeline.yaml`, snpEff 5.4c, database
GRCh38.115). No patient data.

**1. What each policy keeps** (`picks.tsv`). MANE Select per NCBI's MANE v1.5 summary (National
Center for Biotechnology Information & EMBL-EBI, n.d.) against what each snpEff policy keeps:

| Gene | MANE Select (v1.5) | `-tag MANE_Select` | `-canon` |
|---|---|---|---|
| BUB1 | ENST00000302759.11 | same | **ENST00000922602.1** |
| BUB1B | ENST00000287598.11 | same | **ENST00000918306.1** |
| BUB3 | ENST00000368865.9 | same | same |
| CEP57 | ENST00000325542.10 | same | same |
| CEP192 | ENST00000506447.5 | same | **ENST00000912813.1** |
| TRIP13 | ENST00000166345.8 | same | same |

`-canon` departs from MANE in half the panel, including BUB1B, and each time selects a
version-1 transcript model. No panel gene has a MANE Plus Clinical transcript in v1.5.

**2. Known pathogenic variants** (`clinvar_summary.tsv`, `clinvar_variants.tsv`). The 118 public
ClinVar P/LP variants in the panel (National Center for Biotechnology Information, 2026) — BUB1B 84,
CEP57 19, TRIP13 9, BUB1 6, none in BUB3 or CEP192:

| | All transcripts | `-canon` | MANE Select |
|---|---|---|---|
| HIGH-impact (LoF) calls | 112 | 112 | 112 |
| Variants whose call differs between policies | — | 0 | 0 |

The policies are indistinguishable on known variants. That is consistent with ascertainment rather
than equivalence: ClinVar variants are described on the transcripts submitters use — "RefSeq
transcripts are typically used for variant submissions to ClinVar" (Morales et al., 2022) — so they
sit in exons every sensible policy shares. **This benchmark cannot choose a policy.** It does show
something else: 6 of the 118 are not LoF under any policy — 5 missense (MODERATE impact; one also
splice-region) and 1 intronic splice-region variant (LOW) — so a LoF-only filter will miss
pathogenic alleles in these genes.

**3. Where the choice would matter for a novel variant** (`territory.tsv`). An invented substitution
at every position of each gene, counting base pairs where a policy would permit a LoF call
(coding sequence or canonical splice site, ±1/±2):

| Gene | All | `-canon` | MANE | All, not MANE | `-canon`, not MANE | MANE, not `-canon` |
|---|---|---|---|---|---|---|
| BUB1 | 3,472 | 3,442 | 3,354 | 118 | 88 | 0 |
| BUB1B | 3,442 | 3,347 | 3,241 | 201 | 126 | **20** |
| BUB3 | 1,036 | 1,013 | 1,013 | 23 | 0 | 0 |
| CEP57 | 1,567 | 1,543 | 1,543 | 24 | 0 | 0 |
| CEP192 | 7,789 | 7,789 | 7,786 | 3 | 3 | 0 |
| TRIP13 | 1,426 | 1,347 | 1,347 | 79 | 0 | 0 |
| **Total** | **18,732** | **18,481** | **18,284** | **448** | **217** | **20** |

448 bp across the panel can receive a LoF call only from a non-MANE transcript — the "partial LoF"
territory the literature flags as false-positive-prone. `-canon` admits 217 bp of it, and in BUB1B
it also **loses 20 bp of MANE Select sequence**: a truncating variant there would not be called LoF
under `-canon` at all.

---

## Assessment

| | Known P/LP recovered | False-positive exposure | False-negative exposure | Literature | Verdict |
|---|---|---|---|---|---|
| **A. All** | 112/118 | 448 bp of non-MANE territory, where partial LoF concentrates | none within annotated transcripts | against, as a primary call | secondary pass only |
| **B. `-canon`** | 112/118 | 217 bp non-MANE | **20 bp of MANE in BUB1B** | length heuristic, not curated | **reject** |
| **C. MANE** | 112/118 | none beyond MANE | variants in exons absent from MANE — rare (0.37% of pathogenic), not zero | designed default | primary |
| **D. Tiered** | 112/118 | contained — non-MANE LoF is flagged, never the call | none within annotated transcripts | matches PVS1's relevance test | **recommended** |

C alone would be defensible. D costs one extra annotation pass and removes C's one weakness,
which matters more here than usual: this is n=1 in a population of roughly 50, and silently
discarding the one variant that explains the child's disease is the worst failure available.
Keeping it — visibly, one tier down — is cheaper than missing it.

---

## Recommendation — adopted 2026-09-11 as decision D5

1. **Primary LoF call on MANE Select** (`-tag MANE_Select`). In GRCh38.115 this reproduces NCBI's
   MANE v1.5 Select transcript for all six genes; re-check that on any database update, and record
   the MANE release beside every call.
2. **Secondary all-transcript pass.** A variant that is LoF only on non-MANE transcripts is kept as
   `lof_tier = non_mane`, naming the transcript. It never becomes a causal call on its own: it needs
   independent evidence the transcript is biologically relevant (Abou Tayoun et al., 2018), for
   which pext is the natural refinement (Cummings et al., 2020).
3. **Never `-canon`.**
4. **Carry the NMD caveat.** A premature stop in the last exon or the 3′-most 50 nt of the
   penultimate exon is not expected to trigger NMD (Abou Tayoun et al., 2018) — flag it rather than
   treat it as equivalent to an early truncation. snpEff adds LOF and NMD tags with `-lof` (per
   `snpEff ann -h`, v5.4c); its NMD rule has not been compared with the PVS1 rule here.
5. **Do not let L0 be LoF-only.** 6 of 118 known pathogenic variants in the panel are not LoF under
   any policy. "No biallelic LoF found" must be reported as exactly that — never as "no causal
   variant".
6. **Scalability.** MANE Plus Clinical is empty for this panel; for a second disease, include it
   (confirm the tag name in the snpEff database first — it has not been tested here).

**Implemented with the decision:** `annotator.transcript_policy: mane_select_tiered` and
`annotator.mane_release` in the config; the `VariantCall` record in
[`../../src/l0_genomics/transcripts.py`](../../src/l0_genomics/transcripts.py), carrying
`transcript`, `mane_release`, `lof_tier`, `annotation_pass` and `transcript_policy`; and
[`../../tests/test_l0_transcripts.py`](../../tests/test_l0_transcripts.py), which checks that the two
passes stay separate — a `non_mane` tier can only come from the secondary pass — and that
`-canon` and a missing policy are refused.

---

## Limits

- **The ClinVar benchmark cannot discriminate** between policies, for the ascertainment reason
  above. The territory scan is the discriminating measurement, and it counts eligible base pairs —
  it says nothing about how likely a variant at any of them is to be pathogenic.
- **snpEff is not VEP + LOFTEE.** Consequence rules differ between annotators (McCarthy et al.,
  2014), and none of LOFTEE's filters (Karczewski et al., 2020) are applied here.
- **MANE and Ensembl are versioned.** The match between `-tag MANE_Select` and NCBI's table was
  checked for release 115 against v1.5 only.
- **Six genes.** The panel is read from `src/l0_genomics/run.py`, so the script extends with it.

---

## References

Full APA entries and the verification log are in [`../references.md`](../references.md).

| Citation | Key |
|---|---|
| Abou Tayoun et al. (2018) | `aboutayoun2018` |
| Cingolani et al. (2012) | `cingolani2012` |
| Cummings et al. (2020) | `cummings2020` |
| Ensembl (2026) | `ensemblcanonical` |
| Hanks et al. (2004) | `hanks2004` |
| Karczewski et al. (2020) | `karczewski2020` |
| MacArthur et al. (2012) | `macarthur2012` |
| McCarthy et al. (2014) | `mccarthy2014` |
| Morales et al. (2022) | `morales2022` |
| National Center for Biotechnology Information (2026) | `clinvar20260905` |
| National Center for Biotechnology Information & EMBL-EBI (n.d.) | `mane_v15` |
| Pozo et al. (2022) | `pozo2022` |
| Richards et al. (2015) | `richards2015` |
| Singer-Berk et al. (2023) | `singerberk2023` |
| SnpEff & SnpSift documentation (n.d.) | `snpeffdocs` |
| Wright et al. (2023) | `wright2023` |
