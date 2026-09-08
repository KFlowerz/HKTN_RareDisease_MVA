# Aneuploidy-Selective Compounds and Supporting Evidence

> **Research input, not a pipeline decision.** Reviewed and citation-checked 2026-09-08.
> Read this header before using anything below.

## Status and how this may be used

**Citations: verified.** All 21 distinct journal citations in this document were resolved against
Crossref and retraction-checked on 2026-09-08 (`scripts/verify_refs.py`). None are retracted. They
are now in [`../references.bib`](../references.bib) with keys, and rendered in
[`../references.md`](../references.md). Cite by key from there, never by the loose
"Author Year, Journal Vol:Page" strings in the prose below.

**Two errors found and corrected in that sweep:**

- **Torres et al. 2007 is *Science* 317:916–924, not *Cell*** as stated in Section 1A/1B. Key
  `torres2007`, doi:10.1126/science.1142210. (A separate Torres et al. 2010 paper *is* in *Cell*;
  the two appear to have been conflated.)
- The **Oncogene 2012** trisomy-7/AICAR reference is **Ly et al.**, *Oncogene* 32:3139–3146, key
  `ly2012` — first author unnamed in the prose below.

**Scope limits — this document presupposes two decisions the pipeline exists to make.**
Its original title named "BUB1B & TRIP13 Genotypes" and its framing assumes
`therapeutic_endpoint = chemoprevention`. Both are `null` in
[`../../config/pipeline.yaml`](../../config/pipeline.yaml) by hard constraint 5, and `causal_gene`
comes from L0. Read every genotype-specific passage as **conditional** — *if* L0 calls BUB1B or
TRIP13, *then* the following applies. Nothing here may set either value.

**Channel scoping — required to protect the design.** This document may seed the **Channel E**
literature prior, which is that channel's job. It must **not** additionally select Channel C's
signatures or weight Channel C's targets: L2's guardrail requires channels to be independent
because L3's rank aggregation treats cross-channel convergence as evidence, and two channels drawing
on one source would manufacture that convergence.

**Benchmark contamination — do not do what Stage 0 suggests.** Recommendations below propose seeding
the Channel E prior *and* drawing the L4 benchmark positive-control set from the same compound list.
That would make the blinded benchmark measure whether the pipeline recovers what it was told, which
is not a validation. There is no external ground truth for Track 2, so the benchmark is the only
honest evidence the pipeline works. **Positive controls must have provenance disjoint from the
channel priors, and each compound's source must be recorded.**

**Still unverified:** the GEO/ArrayExpress accessions (GSE47830, GSE39768, GSE247267, GSE92742,
GSE70138, GSE106127), the clinical-trial identifiers, the patent number, and the regulatory claims
(approval dates, pediatric indications, BBB penetration). Verify before any of them drives a
pipeline step or reaches the report. The document's own Caveats section flags two of these honestly.

**Coverage:** serves Channels C and E. Contributes nothing to Channel A (knowledge graph),
Channel B (network proximity), or Channel D (phenotype) — including nothing on the features
that actually drive Channel D.

**Quoted third-party text** (Orphanet, NEJM, and others) appears below. Quotation is fine for
internal reasoning; check licensing before any of it reaches a CC-BY-4.0 output.

---


## TL;DR
- The strongest, most-cited aneuploidy-selective vulnerabilities are proteotoxic/chaperone stress (HSP90 inhibitor 17-AAG; proteasome inhibitors), energy stress (AICAR/AMPK), autophagy dependence (chloroquine), the RAF/MEK/ERK–CRAF axis, the mitotic kinesin KIF18A dependency, and SAC/TTK inhibition (reversine) — but nearly all evidence comes from transformed, frequently p53-null aneuploid cancer cells, so translation to a constitutionally-aneuploid, non-transformed child requires caution; only a few (AICAR as a clone-selection proof-of-concept, hydroxychloroquine, and pediatric-approved selumetinib) are plausibly repurposable for chronic pediatric chemoprevention.
- For Channel C (signature reversion) there is **no MVA patient RNA-seq**, but usable public proxy "disease-state" aneuploidy signatures exist with confirmed accessions: Dürrbaum 2014 (GSE47830 / E-GEOD-47830), Stingele 2012 (GSE39768 / E-GEOD-39768), Sheltzer 2013 (TRI70/CIN70/HET70 gene lists), the Ben-David RPE1-hTERT isogenic aneuploid clones (GSE247267), plus LINCS L1000 (GSE92742/GSE70138/GSE106127) for reversine/SAC-perturbation reference signatures.
- For the chemoprevention endpoint, BUB1B- and TRIP13-mutant MVA carry the high embryonal-tumor risk (Wilms tumor, rhabdomyosarcoma, leukemia; cancer usually before age 2–3), whereas CEP57/type-2 MVA lacks increased cancer risk — making SAC-severity a genotype stratifier and TRIP13 itself both a disease gene and (in the opposite direction) a sporadic-cancer drug target.

## Key Findings

1. **Aneuploidy imposes recurrent, druggable stresses.** Independent of which chromosome is gained, aneuploid cells share proteotoxic stress (impaired HSF1/HSP90 folding), increased autophagy/lysosomal load, energy/metabolic stress, replication stress and DNA-damage-response activation, elevated RAF/MEK/ERK signaling, and heightened dependence on faithful mitosis (SAC, KIF18A). Each is a candidate vulnerability.

2. **The best-validated selective compounds** are 17-AAG (HSP90), AICAR (AMPK/energy stress), and chloroquine (autophagy) from the Amon lab (Tang, Williams, Siegel & Amon 2011, Cell 144:499–512: "We have identified the energy stress-inducing agent AICAR, the protein folding inhibitor 17-AAG, and the autophagy inhibitor chloroquine as exhibiting this property. AICAR induces p53-mediated apoptosis in primary mouse embryonic fibroblasts (MEFs) trisomic for chromosome 1, 13, 16, or 19"); SAC/TTK inhibitors (reversine, MPS1-IN-1, MPI-0479605) and KIF18A inhibitors from the Ben-David/DepMap work; CRAF/RAF-MEK-ERK inhibitors from Ben-David/Santaguida; proteasome inhibitors; and trisomy-1q-specific nucleotide analogs (3-deazauridine, RX-3117) from the Sheltzer ReDACT work.

3. **Most evidence is in transformed cells.** The Amon trisomic-MEF work, the Pallister-Killian AICAR data, and the non-transformed Ben-David RPE1 clones are the main exceptions; the SAC/KIF18A/CRAF vulnerabilities are largely documented in p53-mutant cancer lines. This is the central translational caveat for MVA chemoprevention.

4. **Genotype matters for the cancer endpoint.** BUB1B and TRIP13 (severe SAC impairment) confer high embryonal-tumor risk; CEP57 (type 2) does not (Orphanet: "No individual with a CEP57 mutation has thus far been diagnosed with cancer"). This maps onto the pipeline's decision to group BUB1B with TRIP13.

5. **Proxy signatures are available and specific.** Confirmed GEO/ArrayExpress accessions let the team build a "disease-state" transcriptomic signature to reverse even without patient RNA.

## Details

### SECTION 1 — ANEUPLOIDY-SELECTIVE COMPOUNDS

**A. Proteotoxic / chaperone (HSP90) stress**
- **17-AAG (tanespimycin; 17-allylamino-17-demethoxygeldanamycin)**, geldanamycin, radicicol. Target: HSP90. Vulnerability: aneuploid cells have impaired HSF1/HSP90-dependent folding and proteotoxic stress; extra chromosomes burden protein-quality control. Systems: disomic budding yeast sensitive to radicicol/geldanamycin (Torres et al. 2007, Cell); trisomic MEFs (chr 1, 13, 16, 19) more sensitive to 17-AAG than wild-type (Tang, Williams, Siegel & Amon 2011, Cell 144:499–512); aneuploid human cancer lines (17-AAG, especially combined with AICAR). Mechanistic support: Donnelly/Dürrbaum/Storchová 2014, EMBO J 33:2374 (HSF1 deficiency + impaired HSP90 folding are hallmarks of aneuploid human cells).
- **Proteasome inhibitors (bortezomib, carfilzomib; MG132 as tool).** Target: 26S proteasome β5 chymotryptic site. Vulnerability: aneuploidy raises misfolded-protein/UPR load. Systems: disomic yeast sensitive to MG132 (Torres 2007); Cohen-Sharir/Ben-David flagged "increased sensitivity to proteasome inhibition" as a candidate vulnerability in the highly-aneuploid cancer-line analysis (Nature 2020, Supplementary Table 3). Selectivity in constitutional (non-transformed) aneuploidy is not established.

**B. Energy / metabolic stress**
- **AICAR (acadesine).** Target: AMPK activation (via the AMP analog ZMP), mimicking energy stress; also drives EGFR degradation in trisomy-7 colonic cells. Vulnerability: aneuploid cells have altered metabolism/energy stress. Systems: p53-mediated apoptosis in trisomic MEFs (Tang & Amon 2011); trisomy-7 human colonic epithelial cells via EGFR degradation (Oncogene 2012, doi 10.1038/onc.2012.339); selective reduction of isochromosome-12p mosaic cells in Pallister-Killian syndrome patient-derived fibroblasts (bioRxiv 2022) — notable as a *constitutional-mosaicism* clone-normalization proof of concept directly relevant to MVA. AICAR/17-AAG combination is synergistic against high-grade aneuploid cells.

**C. Autophagy dependence**
- **Chloroquine / hydroxychloroquine.** Target: lysosomal autophagy inhibition. Vulnerability: aneuploid cells upregulate p62-dependent autophagy and lysosomal degradation and rely on it to clear protein aggregates (Stingele 2012/2013; Santaguida, Vasile, White & Amon 2015, Genes Dev 29:2010 — aneuploidy-induced stresses limit autophagic degradation). Systems: apoptosis in trisomic MEFs (Tang & Amon 2011).

**D. RAF/MEK/ERK — CRAF**
- **CRAF (RAF1) inhibition; pan-RAF (e.g., LY3009120, PLX7904), MEK inhibitors (trametinib, selumetinib).** Vulnerability: aneuploid cells show elevated RAF/MEK/ERK activity and depend particularly on CRAF to buffer increased DNA damage. Systems: non-transformed, p53-WT RPE1-hTERT isogenic aneuploid clones + validation in cancer lines (Zerbib et al. 2024, Nat Commun 15:7772). CRAF/MEK inhibition sensitizes aneuploid cells to DNA-damaging chemotherapy and PARP inhibitors; olaparib resistance correlated with high RAF/MEK/ERK signaling in highly-aneuploid tumors.

**E. Mitotic dependencies — SAC/TTK and KIF18A**
- **SAC/TTK(MPS1) inhibitors: reversine, MPS1-IN-1, MPI-0479605.** Vulnerability and paradox: highly-aneuploid cancer cells depend on *genetic* perturbation of core SAC genes yet are *less* sensitive to short-term *chemical* SAC inhibition (Cohen-Sharir et al. 2021, Nature 590:486–491: "aneuploid cancer cells show increased sensitivity to genetic perturbation of core components of the spindle assembly checkpoint (SAC)… we also found that aneuploid cancer cells were less sensitive than diploid cells to short-term exposure to multiple SAC inhibitors"). Systems: aneuploidy landscape of ~1,000 cancer lines + Achilles/PRISM (repurposing 19Q4) + isogenic validation. *Caveat for MVA:* BUB1B/TRIP13 patients already have SAC impairment, so further SAC inhibition is mechanistically counterproductive for baseline cells.
- **KIF18A inhibitors: AMG 650 (Amgen), sovilnesib, VLS-1488 (Volastra), ATX020.** Target: mitotic kinesin KIF18A. Vulnerability: CIN-high/aneuploid/whole-genome-doubled cancer cells depend on KIF18A to avoid lethal multipolar mitosis, while KIF18A is dispensable in normal euploid division (Cohen-Sharir 2021: "Aneuploid cancer cells were particularly vulnerable to depletion of KIF18A"). Systems: DepMap mining + CIN-high cancer lines, in vivo xenografts; AMG 650 phase 1 (NCT04293094) in TP53-mutant/CIN solid tumors; VLS-1488 in phase 1/2; reported minimal effect on normal bone marrow mononuclear cells.

**F. Karyotype-specific "aneuploidy addiction" (ReDACT / Sheltzer)**
- **3-deazauridine and RX-3117 (fluorocyclopentenyl-cytosine).** Target: nucleotide analogs activated by UCK2 (pyrimidine salvage kinase encoded on chr 1q). Vulnerability: 1q-trisomy raises UCK2, creating collateral sensitivity (Girish et al. 2023, Science 381:eadg4521 — chr1q UCK2 phosphorylates these compounds to cytotoxic derivatives; 1q gain also "increases the expression of MDM4 and suppresses p53 signaling," and "TP53 mutations are mutually exclusive with 1q aneuploidy"). Systems: ReDACT isogenic 1q-trisomy vs disomy in A2780/MCF10A; MPS1-inhibitor–derived DLD1 trisomy-1; NCI-60 UCK2 correlation. Conceptually the closest to a "hold the expanding clone in a low-malignant state" chemoprevention strategy.
- **AHR-directed (CGS-15943)** for chr7p (AHR) gains — noted in the associated patent (US20250102493A1) as a complementary karyotype-specific approach.

**G. Sphingolipid / ceramide metabolism**
- Aneuploid yeast are sensitive to perturbation of sphingolipid metabolism (cited as yeast precedent within Cohen-Sharir 2021). **Myriocin** (serine palmitoyltransferase inhibitor) and **FTY720/fingolimod** (sphingosine analog; also inhibits ceramide synthases, Ki ≈ 2.15 µM for CerS2, and reactivates PP2A) are the tool/clinical agents. Direct selectivity in mammalian *constitutional* aneuploidy is not established; FTY720's anticancer activity (e.g., t(8;21) AML) proceeds via PP2A/ceramide.

### SECTION 2 — DRUGGABILITY, APPROVAL, PEDIATRIC SAFETY

- **17-AAG (tanespimycin):** investigational; multiple adult oncology trials, never approved; hepatotoxicity; limited BBB penetration; **unsuitable for chronic pediatric chemoprevention** (cytotoxic, tolerability). Retain as positive control.
- **AICAR (acadesine):** investigational; studied in CLL and as a cardioprotectant in cardiac surgery; poor oral bioavailability. Best value here is as a *proof-of-concept clone-selection* agent (Pallister-Killian data) rather than a finished pediatric preventive.
- **Chloroquine / hydroxychloroquine:** FDA/EMA-approved (malaria, autoimmune); extensive **pediatric safety record**; crosses the BBB; chronic dosing feasible with retinal/cardiac monitoring. **Among the most plausibly repurposable** agents for chronic use in a child on the proteostasis/autophagy axis.
- **Bortezomib / carfilzomib:** approved (myeloma/MCL); pediatric use exists within ALL protocols; peripheral neuropathy (bortezomib), cardiotoxicity (carfilzomib); poor CNS penetration; cytotoxic — **unsuitable for chronic prevention**.
- **MEK inhibitors — selumetinib (KOSELUGO):** **FDA-approved in children** — originally April 10, 2020 for pediatric patients ≥2 years with symptomatic, inoperable NF1 plexiform neurofibromas, with the label subsequently broadened (FDA, September 10, 2025) to "pediatric patients 1 year of age and older"; dose 25 mg/m² twice daily; warnings include cardiomyopathy and ocular toxicity. This makes selumetinib the **most pediatric-ready RAF/MEK/ERK-axis candidate**. Trametinib is also approved (adult/pediatric contexts). CRAF-selective agents remain investigational.
- **KIF18A inhibitors (AMG 650, sovilnesib, VLS-1488):** investigational, adult phase 1/2 only; no pediatric data; oral; designed to spare normal cells — mechanistically attractive "watch" class.
- **Reversine / MPS1 inhibitors:** tool/investigational; mechanistically contraindicated for baseline SAC-deficient MVA cells.
- **3-deazauridine / RX-3117:** investigational cytotoxic nucleotide analogs (RX-3117 reached adult oncology trials); not chronic-prevention agents, but the conceptual model for a clone-selective "anti-trisomy" approach.
- **FTY720 / fingolimod:** approved for relapsing MS **including a pediatric MS indication (≥10 years)**; crosses BBB; chronic oral dosing; bradycardia/immunosuppression/macular-edema risks. Pediatric-ready pharmacology but weak aneuploidy-selectivity evidence.
- **Myriocin:** tool compound only.

**Central caveat:** the selectivity data derive overwhelmingly from transformed, often p53-null, chromosomally-unstable *cancer* cells. Constitutional MVA cells are non-transformed and (in BUB1B/TRIP13) have functional-but-impaired checkpoints. The evidence bearing on *premalignant/expanding-clone* selectivity vs baseline aneuploid tissue is (a) the Sheltzer "aneuploidy addiction" concept — an anti-trisomy compound could keep a population "in a low-malignant state" — and (b) the AICAR Pallister-Killian clone-normalization data. p53 status is a key effect-modifier: AICAR's pro-apoptotic effect is p53-mediated in primary MEFs, yet AICAR/17-AAG still act in p53-null tumor cells.

### SECTION 3 — PROXY DISEASE-STATE SIGNATURES (Channel C)

- **Dürrbaum et al. 2014, BMC Genomics 15:139** — "uniform/model aneuploidy response" in human cells (trisomic/tetrasomic HCT116 & RPE-1-derived lines). **Accession: GSE47830 (E-GEOD-47830).** Best single "common aneuploidy response" human signature.
- **Stingele et al. 2012, Mol Syst Biol 8:608** — genome/transcriptome/proteome response to aneuploidy (HCT116 & RPE-1 trisomies/tetrasomies). **Accession: GSE39768 (E-GEOD-39768).**
- **Sheltzer 2013, Cancer Res 73:6401–6412** — transcriptional/metabolic signature of primary aneuploidy present in CIN cancers; defines **TRI70/HET70** and uses **CIN70** (Carter 2006). Gene lists are in the supplements (no dedicated GEO deposit was located). Good for prognostic overlay and as a reversal target set.
- **Ben-David RPE1-hTERT isogenic aneuploid clones (Zerbib et al. 2024, Nat Commun 15:7772) — RNA-seq accession GSE247267.** Non-transformed, p53-WT, karyotype-defined near-diploid→aneuploid clones — the **most MVA-relevant human model** and the ideal validation substrate.
- **LINCS/L1000: GSE92742 (Phase I, chemical+genetic), GSE70138 (Phase II), GSE106127 (genetic shRNA/CRISPR subset).** Query these for reversine / MPS1-TTK perturbation signatures to build the SAC-perturbation reference (reversine = Mps1 inhibitor; Santaguida et al. 2010, J Cell Biol 190:73). The presence of the specific reversine perturbagen must be confirmed in the compound metadata (clue.io/iLINCS).
- **Additional models to consider:** Santaguida et al. 2017 (Dev Cell 41:638) Mps1i-treated RPE1 arrested-aneuploid-cell signatures; trisomic MEF arrays (Williams et al. 2008). BUB1B/TRIP13 loss-of-function *transcriptomic* signatures do not exist off-the-shelf (Yost 2017 and Hanks 2004 are functional/cytogenetic); derive a TRIP13/BUB1B-knockdown signature de novo from GSE106127 shRNA data or generate it in the RPE1 system.

### SECTION 4 — MVA GENOTYPE-SPECIFIC CANCER EPIDEMIOLOGY

- **Overall malignancy proportion:** a widely used figure is **~37% of reported MVA patients develop cancer** — rhabdomyosarcoma, Wilms tumor, and leukemia, mostly within the first ~3 years of life, sometimes in utero (traceable to the Kops-lab review literature, Suijkerbuijk et al. 2010, Cancer Res 70:4891). **Important disambiguation:** this proportion must not be confused with the separate statement in Rio Frio et al. 2010, NEJM 363:2628 that "Worldwide, 37 cases of the mosaic variegated aneuploidy syndrome have been reported" — that "37" is a *case count*, not a percentage. Treat the ~37% cancer-proportion as approximate given ultra-small samples.
- **BUB1B-mutant (MVA1):** high cancer incidence — Orphanet (reviewers Hanks, Rahman, Snape): "Individuals with BUB1B mutations have a high incidence of cancer (approximately 75%)." Tumor spectrum: embryonal rhabdomyosarcoma, Wilms tumor (nephroblastoma), leukemia (ALL), myelodysplasia; cancer usually diagnosed **before age 2** (NEJM 2010). Original gene identification: Hanks et al. 2004, Nat Genet 36:1159 (biallelic BUB1B in five families, two with embryonal RMS — "the first to relate germline mutations in a spindle checkpoint gene with a human disorder"). Adult-onset extension: Rio Frio et al. 2010, NEJM — homozygous intronic BUB1B mutation with GI adenocarcinomas (ampulla of Vater, then colon and stomach) in a 68-year-old, expanding the phenotype to common adult cancers.
- **TRIP13-mutant (MVA3):** biallelic loss-of-function — Yost et al. 2017, Nat Genet 49:1148–1151: "we identified six individuals with biallelic loss-of-function mutations in TRIP13. All six developed Wilms tumor… patient cells have no detectable TRIP13 and have substantial impairment of the spindle assembly checkpoint (SAC)." MVA features (mosaic aneuploidy, microcephaly, developmental delay, seizures) were variably present; a sibling of one proband died at age 4 of a pelvic Sertoli–Leydig cell tumor.
- **CEP57-mutant (MVA2):** centrosomal/microtubule-stabilizing protein; **no increased cancer risk documented** (Orphanet: "No individual with a CEP57 mutation has thus far been diagnosed with cancer"; NEJM 2010: in children "who have neither BUB1B mutations nor premature chromatid separation, cancer is rarely, if ever, present"). CEP57 cells show minimal SAC deficiency (Yost 2017) — the mechanistic basis for the tumor-risk split and validation that SAC-severity, not aneuploidy per se, tracks embryonal-tumor risk.
- **Other/rarer genes:** BUB1, BUB3, MAD2L1BP/p31^comet (associated with juvenile granulosa cell tumors), CEP192, SMC5 — expanding the genotype list.
- **Secondary downstream target set (associated pediatric tumors):**
  - *Wilms tumor* (Scott/Rahman stratification, 120 tumors): 11p15 IGF2/H19 imprinting defects (~69%; H19 epimutation ~37%, paternal UPD ~32%), WTX/AMER1 (~32%), CTNNB1/β-catenin (~15%), WT1 (~12%), TP53 (~5%). Actionable nodes: WNT/β-catenin and IGF2/IGF1R axes.
  - *Rhabdomyosarcoma:* alveolar/fusion-positive driven by PAX3-FOXO1 or PAX7-FOXO1 with transcriptional targets MET, MYCN, IGF1R, FGFR4, ALK; embryonal/fusion-negative driven by RAS-pathway mutations and TP53 with low mutational burden. MVA cases are predominantly *embryonal* RMS, so RAS/RTK and IGF axes are most relevant.
  - *Leukemia / MDS:* monosomy-7 MDS and ALL reported; aneuploidy-driven with no single recurrent driver.

### SECTION 5 — TRIP13-SPECIFIC EVIDENCE

- **Role:** TRIP13 is an AAA+ ATPase that, with the adapter p31^comet, disassembles MAD2-containing mitotic-checkpoint complexes to silence the SAC and permit mitotic exit; it also functions in meiotic/mitotic DNA double-strand-break repair and HORMAD regulation. It regulates **both activation and inactivation** of the SAC (Ma & Poon 2016, Cell Rep 14:1086; structural mechanism in Alfieri, Chang & Barford 2018, Nature 559:274).
- **In MVA:** biallelic loss-of-function → no detectable TRIP13 → severe SAC impairment and a high rate of chromosome missegregation; restoring TRIP13 function rescues SAC proficiency and accurate segregation (Yost 2017). This groups TRIP13 mechanistically with BUB1B — both cause severe SAC impairment and both confer embryonal-tumor predisposition — which is the rationale for the pipeline's BUB1B+TRIP13 pairing.
- **As a sporadic-cancer drug target (opposite direction — TRIP13 is often *over*expressed and oncogenic in sporadic tumors):** small-molecule inhibitors **DCZ0415** (first-in-class, structure-based; direct binding confirmed by NMR/SPR; active in multiple myeloma and in colorectal cancer regardless of p53/KRAS/BRAF/EGFR/MSI status; induces G2/M arrest and apoptosis; reduces xenograft growth and metastasis), **TI17**, **DCZ5417**, and **DCZ5418**. TRIP13 is oncogenic in multiple myeloma, brain tumors, lung, renal, and colorectal cancers.
- **Nuance for MVA chemoprevention:** in constitutional TRIP13-null MVA, baseline cells already *lack* TRIP13, so a TRIP13 *inhibitor* is **not** a chemoprevention strategy for those cells. The relevant strategies are (a) TRIP13-function *restoration*/gene therapy (Yost demonstrated phenotypic rescue in patient cells and oocytes) and (b) general aneuploidy-selective agents exploiting the resulting CIN. TRIP13 inhibitors belong to the *downstream-tumor* target set only if a specific tumor over-expresses TRIP13.

### COMPACT CANDIDATE COMPOUND TABLE

| Compound (class/target) | Vulnerability exploited | Key evidence system | Approval / pediatric-safety status | Chemoprevention suitability |
|---|---|---|---|---|
| **17-AAG / tanespimycin** (HSP90 inhibitor) | Proteotoxic/chaperone stress | Trisomic MEFs; disomic yeast; aneuploid cancer lines (Tang & Amon 2011) | Investigational; hepatotoxic; poor BBB | Low — cytotoxic; positive control only |
| **AICAR / acadesine** (AMPK activator) | Energy/metabolic stress | Trisomic MEFs; trisomy-7 colon cells; Pallister-Killian mosaic fibroblasts | Investigational; poor oral bioavailability | Moderate as clone-selection proof-of-concept; not a finished drug |
| **Chloroquine / hydroxychloroquine** (autophagy inhibitor) | Autophagy dependence | Trisomic MEFs (Tang & Amon 2011); Storchová/Santaguida mechanistic | **Approved; pediatric-safe; CNS-penetrant** | **High** — repurposable for chronic use w/ monitoring |
| **Bortezomib / carfilzomib** (proteasome) | Proteostasis/UPR load | Disomic yeast (MG132); DepMap candidate | Approved (myeloma/MCL); neuropathy/cardiotox; poor CNS | Low — cytotoxic |
| **Selumetinib** (MEK1/2) | RAF/MEK/ERK–CRAF dependence | RPE1 aneuploid clones + cancer lines (Zerbib 2024) | **FDA-approved in children (NF1, ≥1 yr)**; cardiac/ocular warnings | **Highest pediatric readiness** on this axis |
| **CRAF-selective / pan-RAF inhibitors** | RAF/MEK/ERK–CRAF dependence | RPE1 aneuploid clones (Zerbib 2024) | Investigational | Moderate — mechanism strong, agents immature |
| **KIF18A inhibitors** (AMG 650, sovilnesib, VLS-1488, ATX020) | Mitotic/multipolar-spindle dependence in CIN-high cells | DepMap + CIN-high lines; adult phase 1/2 | Investigational; no pediatric data; spares normal cells | Watch class — promising, pediatric-immature |
| **Reversine / MPS1(TTK) inhibitors** (MPS1-IN-1, MPI-0479605) | SAC dependence (genetic) | ~1,000 cancer lines + PRISM/Achilles (Cohen-Sharir 2021) | Tool/investigational | Contraindicated for baseline SAC-deficient MVA cells |
| **3-deazauridine / RX-3117** (UCK2-activated nucleotide analogs) | 1q-trisomy "aneuploidy addiction" via UCK2 | ReDACT isogenic 1q clones (Girish 2023) | Investigational cytotoxic antimetabolites | Low as chronic drug; key conceptual anti-trisomy model |
| **FTY720 / fingolimod, myriocin** (sphingolipid/ceramide) | Sphingolipid-metabolism dependence | Yeast precedent; FTY720 anticancer via PP2A/ceramide | Fingolimod approved (MS, pediatric ≥10 yr); myriocin tool only | Low-moderate; selectivity evidence weak |
| **TRIP13 inhibitors** (DCZ0415, TI17, DCZ5417/5418) | TRIP13-overexpressing sporadic tumors | Myeloma/CRC/renal cell lines + xenografts | Investigational | Not for baseline TRIP13-null MVA; downstream-tumor target only |

## Recommendations

**Stage 0 — Build the priors and signatures now (weeks).**
- Seed the Channel E "aneuploidy-stress prior" with the mechanism classes weighted by strength/breadth of evidence: highest weight to proteotoxic/HSP90, autophagy, energy/AMPK, RAF-MEK-ERK/CRAF, KIF18A, SAC; secondary weight to sphingolipid and karyotype-specific (1q/UCK2, 7p/AHR) axes.
- Build the Channel C proxy disease signature by combining **GSE47830 (Dürrbaum)**, **GSE39768 (Stingele)**, and **GSE247267 (Ben-David RPE1 clones)**, intersecting to a chromosome-agnostic "common aneuploidy response" core; retain **TRI70/CIN70** as an orthogonal prognostic overlay; derive a **TRIP13/BUB1B-loss signature** de novo from GSE106127 shRNA data or generate it in the RPE1 system.
- Positive-control set for pipeline validation: 17-AAG, AICAR, chloroquine (Amon triad), reversine (SAC), a KIF18A inhibitor, a MEK inhibitor (selumetinib/trametinib), and 3-deazauridine (1q-specific). Negative controls: compounds showing no aneuploidy selectivity in Tang 2011.

**Stage 1 — Prioritize pediatric-viable chemoprevention candidates.** Rank by (chronic pediatric safety × BBB penetration × strength of aneuploidy-selectivity evidence). Front-runners for a chronic preventive concept: **hydroxychloroquine** (approved, pediatric-safe, CNS-penetrant, autophagy axis) and **selumetinib** (pediatric-approved ≥1 yr, MEK axis). Treat KIF18A inhibitors as a "watch" class pending pediatric data. Explicitly exclude cytotoxic-only agents (17-AAG, proteasome inhibitors, nucleotide analogs, reversine) from chronic-prevention consideration; retain them as positive controls or acute-intervention concepts only.

**Stage 2 — De-risk the transformation caveat.** Before any translational claim, require validation in a *non-transformed, p53-WT, constitutionally-aneuploid* model (Ben-David RPE1 clones; ideally BUB1B/TRIP13-hypomorph patient lymphoblasts/fibroblasts) demonstrating selective growth-suppression of aneuploid/expanding clones without harming euploid cells — replicating the AICAR/Pallister-Killian clone-normalization result.

**Benchmarks that change the plan:**
- If a candidate fails to show clone-selectivity in a p53-WT non-transformed model → drop from chemoprevention; keep only as an anti-tumor agent for the downstream target set.
- If reversine/SAC or a TRIP13 inhibitor scores highly in Channel C reversion → down-weight for MVA baseline cells (mechanistically contraindicated given pre-existing SAC deficiency) but retain for the *downstream tumor* target set (e.g., TRIP13-high sporadic tumor).
- If a genotype-stratified readout separates BUB1B/TRIP13 (SAC-severe) from CEP57 (SAC-mild) signatures → prioritize the SAC-severe arm for the cancer endpoint, consistent with the observed tumor-risk split.

## Caveats
- The evidence base is dominated by transformed, frequently p53-null, CIN cancer cells; extrapolation to a non-transformed child is uncertain. The non-transformed anchors (trisomic MEFs; Pallister-Killian AICAR clone-selection; Ben-David p53-WT RPE1 clones) are the exceptions to lean on.
- MVA is ultra-rare; epidemiological percentages (~37% overall developing cancer; ~75% for BUB1B) come from case aggregations, not prospective cohorts, and carry wide uncertainty. The ~37% cancer *proportion* (Suijkerbuijk 2010) must not be conflated with the "37 cases reported worldwide" *count* (NEJM 2010).
- Aneuploidy-targeting therapy carries a theoretical risk of *promoting* transformation of non-aneuploid cells (by exacerbating CIN) or of selecting resistant near-euploid/less-dependent clones (flagged in the Royal Society Open Biology review, 2020) — a real concern for a lifelong preventive.
- SAC- and TRIP13-directed inhibition is mechanistically counterproductive for baseline BUB1B/TRIP13-deficient cells (already impaired SAC); these belong to the downstream-tumor target set, not baseline chemoprevention.
- No LINCS reversine signature was individually confirmed; the master series (GSE92742/GSE70138) are the correct place to query but presence of the specific perturbagen needs verification.
- No dedicated GEO deposit was located for Sheltzer 2013; use its published TRI70/CIN70/HET70 gene lists.
- BBB-penetration and pediatric-safety statements vary by agent and dose and must be re-verified against current regulatory labels before any clinical step.