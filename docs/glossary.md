# Glossary — terminology, abbreviations and acronyms

**Rare Disease, Real Kid: MVA Hackathon 2026**

Every term, abbreviation and acronym this project uses in its reports, its code and its
decision log. Each entry says what the term means **in this project**, which is not always
the broadest meaning it has elsewhere.

Nothing here is patient data: no coordinates, no variant identifiers, no phenotype terms.

**Citations** are given in-text in APA 7th. The full entries, with the Crossref verification
log and the date each source was checked, are in [`references.md`](references.md) — the
canonical reference list for the whole project.

**Version numbers** shown are the ones pinned at the time of writing;
[`config/pipeline.yaml`](../config/pipeline.yaml) is authoritative, and
`results/_manifest.json` records what actually ran.

---

## 1. The disease

| Term | Stands for | What it means here |
|---|---|---|
| **MVA** | Mosaic variegated aneuploidy | The patient's disease. OMIM 257300. Biallelic loss of function in a spindle-assembly-checkpoint gene lets cells divide with the wrong number of chromosomes, differently in different cells (Hanks et al., 2004) |
| **Aneuploidy** | — | A cell having the wrong number of chromosomes — not the normal 46 in a human cell |
| **Mosaic** | — | Present in some cells and not others, rather than in every cell of the body |
| **Variegated** | — | Different cells are affected in *different* ways — one cell gains a chromosome, another loses a different one |
| **SAC** | Spindle assembly checkpoint | The mechanism that halts cell division until the chromosomes are correctly lined up. MVA is what happens when it fails |
| **BUB1B** | — | The gene this project treats as causal. One of the genes that builds the SAC |
| **BUBR1** | — | The protein *BUB1B* codes for. Gene names are italicised, protein names are not |
| **Panel genes** | — | The six SAC genes L0 searches: *BUB1B*, *CEP57*, *TRIP13*, *BUB1*, *BUB3*, *CEP192*. Their supporting evidence differs and is tabulated in [`report_track1.md`](report_track1.md) (Snape et al., 2011; Yost et al., 2017; Carvalhal et al., 2022; Guo et al., 2024) |
| **Cancer predisposition** | — | An inherited raised risk of developing cancer. In MVA it is substantial, and it is the reason safety here is a gate rather than a score (Scott et al., 2006) |
| **RAF/MEK/ERK** | — | A chain of proteins relaying growth signals inside a cell. Cells with the wrong chromosome count lean on it to cope with the DNA damage they accumulate, which is the published rationale behind the one literature-supported candidate (Zerbib et al., 2024) |

---

## 2. Genetics and reading a genome

| Term | Stands for | What it means here |
|---|---|---|
| **WGS** | Whole-genome sequencing | Sequencing the entire genome rather than selected genes. What the dataset contains |
| **VCF** | Variant Call Format | The text file listing how a genome differs from the reference. **The pipeline's input** |
| **BAM** | Binary Alignment Map | Sequencing reads aligned to the genome. **The dataset ships none**, which is why depth-based methods were unavailable |
| **CSV** | Comma-separated values | A plain table. The format the Track 1 predictions are submitted in |
| **PASS** | — | The VCF `FILTER` value meaning the caller had no reservations about a variant. Alleles that are not `PASS` are flagged, not silently dropped |
| **CDS** | Coding sequence | The part of a transcript that becomes protein. snpEff's `-canon` picks the *longest* CDS, which is **not** the same as MANE Select — the reason that flag is refused |
| **GRCh38** | Genome Reference Consortium human build 38 | The reference genome everything is expressed against |
| **Variant** | — | One difference between this genome and the reference |
| **Allele** | — | One of the two copies of a position, one inherited from each parent |
| **SNV** | Single-nucleotide variant | A one-letter change |
| **CNV** | Copy-number variant | A duplicated or deleted stretch. **Not callable from a VCF alone** |
| **SV** | Structural variant | A large rearrangement. Also not callable from a VCF alone |
| **LoF** | Loss of function | A change that stops the gene producing working protein |
| **Truncating** | — | A change that cuts the protein short |
| **Missense** | — | A change that swaps one building block for another |
| **NMD** | Nonsense-mediated decay | The cell's disposal of messages carrying a premature stop. A truncation in the **last exon** escapes it and so does less damage — which is why the PVS1 criterion is downgraded there (Abou Tayoun et al., 2018) |
| **Exon / intron** | — | The parts of a gene that do and do not end up in the protein recipe |
| **Biallelic** | — | Both copies of the gene affected. MVA is recessive, so this is the unit of evidence |
| **Compound heterozygous** | — | Biallelic by way of **two different** changes, one on each copy — as opposed to the same change twice |
| **Phase** | — | Which parental copy each variant sits on. **Unresolvable here**: it needs parental samples |
| **In *cis* / in *trans*** | — | Both variants on the same copy (*cis*) or one on each (*trans*). Only *trans* is causal for a recessive disease |
| **Unphased** | — | Phase unknown. The largest stated uncertainty in the Track 1 call |
| **Karyotype** | — | The description of a cell's chromosome complement. **Deliberately not published** for this patient (decision D20) |
| **Proband** | — | The affected individual through whom a family comes to medical attention. The reports say *patient* instead |
| **HGVS** | Human Genome Variation Society | The standard notation for naming a variant precisely. Treated as an identifier and never published |
| **MANE** | Matched Annotation from NCBI and EMBL-EBI | One agreed transcript per gene across the two major archives. **MANE Select** is the one the calls are made on (Morales et al., 2022) |
| **Transcript** | — | One of several possible readings of a gene. Which one you pick changes what a variant appears to do |
| **BAF** | B-allele frequency | The proportion of reads supporting the alternate allele. Used to estimate aneuploidy burden without alignments (Conlin et al., 2010; Loh et al., 2018) |
| **FORMAT/AD** | Allelic depth | The VCF field holding per-allele read counts, from which BAF is computed |
| **Mosaic fraction** | — | The proportion of cells carrying an abnormality |
| **ACMG/AMP** | American College of Medical Genetics and Genomics / Association for Molecular Pathology | The standard framework for classifying variants. This pipeline **does not classify** — it reports what others have (Richards et al., 2015) |
| **PVS1** | Pathogenic, very strong, criterion 1 | The ACMG/AMP criterion for loss-of-function variants (Abou Tayoun et al., 2018) |

---

## 3. This project's pipeline

| Term | Stands for | What it means here |
|---|---|---|
| **L0 – L5** | Layer 0 to layer 5 | The six pipeline stages: genomics, target, channels, integrate, validate, report. Each writes files the next one reads |
| **Channel A – E** | — | The five independent candidate generators inside L2: knowledge graph, network proximity, signature reversion, phenotype, literature prior |
| **Artifact** | — | A file the pipeline wrote under `results/`. Every number in every report names the artifact it came from |
| **Seed** | — | The fixed random-number starting point, 42, so a re-run reproduces the same output |
| **`PYTHONHASHSEED`** | — | An environment variable fixed at launch so Python's hashing is deterministic too |
| **Config-driven** | — | Behaviour lives in `config/pipeline.yaml`, not in the code, so a reader can change a threshold and see what moves |
| **Fail closed** | — | When evidence is missing, exclude rather than allow. Why "no safety data found" is an exclusion |
| **Publish guard** | — | The check that refuses to write anything identifying into the published dossier |
| **Dossier** | — | The static HTML output of L5 — one page per surviving candidate, plus the exclusions |
| **Proxy signature** | — | A stand-in for the disease's effect on cells, used because no patient RNA exists. Channel C's, from a gene knockdown, proved to carry no gene-specific signal |
| **Adversarial pass** | — | The reasoning step that deliberately searches for evidence *against* each candidate |
| **Convergence** | — | Two or more independent channels agreeing on the same drug. A **discriminating** channel is one that ranks only a narrow slice; a **broad** channel ranks nearly everything and its agreement alone is not counted |
| **CFTR** | — | The gene behind cystic fibrosis. Not part of this patient's disease — it is the **second disease** the pipeline was re-run on, by changing one config line, to test whether the method is disease-agnostic |

---

## 4. Methods and statistics

| Term | Stands for | What it means here |
|---|---|---|
| **Interactome** | — | The network of which proteins physically or functionally interact |
| **Disease module** | — | The neighbourhood of the interactome around the causal gene — the region L1 defines as the target space |
| **RWR** | Random walk with restart | How L1 grows a module outward from seed proteins, with a restart probability keeping it close to the seeds |
| **Network proximity** | — | How close a drug's targets sit to the disease module, scored against a degree-matched null (Guney et al., 2016) |
| **Degree-matched null** | — | The comparison set: random proteins with the *same* number of connections, so a hub is not mistaken for a finding |
| **RRA** | Robust rank aggregation | How L3 combines several channels' rankings into one (Kolde et al., 2012) |
| **Signature reversion** | — | Looking for drugs whose effect on cells is the opposite of the disease's — the Connectivity Map idea (Lamb et al., 2006) |
| **L1000** | — | The LINCS assay measuring about a thousand genes as a proxy for the whole transcriptome (Subramanian et al., 2017) |
| **shRNA knockdown** | Short hairpin RNA | Switching a gene off in cultured cells, to see what changes. Channel C's proxy query |
| **WTCS / NCS** | Weighted connectivity score / normalised connectivity score | The signature-matching statistic and its normalised form (Subramanian et al., 2017), built on weighted running-sum enrichment (Subramanian et al., 2005) |
| **GSEA** | Gene set enrichment analysis | The running-sum enrichment method underneath the connectivity score (Subramanian et al., 2005) |
| **Null** | — | What the statistic looks like when there is nothing to find, used as the comparison a real result has to beat. **A null that cannot fail is not a null** — three were built and rejected before one could discriminate |
| **Permutation null** | — | A null built by shuffling the data many times |
| **Control-gene null** | — | The null that settled channel C: run the identical procedure on unrelated genes and ask whether the real gene stands out |
| **Benchmark** | — | Testing whether the pipeline recovers a frozen set of already-published compounds. There is no external ground truth for this disease |
| **AUROC** | Area under the receiver operating characteristic curve | A ranking-quality score. 0.5 is chance, 1.0 is perfect |
| **Bootstrap CI** | Confidence interval by resampling | The uncertainty around a score. Reported here because the point estimate alone would imply precision that does not exist |
| **z-score** | — | How many standard deviations a value sits from its baseline |
| **MAD** | Median absolute deviation | A spread measure that a few extreme values cannot distort, used for the aneuploidy baseline |
| **Folded-normal correction** | — | The step that stops small mosaic fractions being systematically overstated, since the raw statistic cannot go below zero |
| **IC / MICA** | Information content / most informative common ancestor | How phenotype-term similarity is measured in channel D |
| **EPCR** | Estimated probability of causal relationship | The confidence column in the Track 1 submission format |
| **F-max** | — | The best F-score across all thresholds, part of the Track 1 scoring |

---

## 5. Data sources

| Source | What it is | Used by |
|---|---|---|
| **snpEff** | Variant annotation software, run **locally** so no coordinate leaves the machine (Cingolani et al., 2012) | L0 |
| **ClinVar** | Public archive of laboratories' verdicts on variants (Landrum et al., 2018) | L0 |
| **gnomAD** | Population allele frequencies from very large cohorts (Chen et al., 2024) | L0 |
| **UniProt** | Protein sequence and domain records (The UniProt Consortium, 2025) | L0, L1 |
| **STRING** | Protein–protein association network (Szklarczyk et al., 2022) | L1, L2 |
| **Reactome** | Curated biological pathways (Milacic et al., 2023) | L1 |
| **Open Targets** | Drug–target, molecule identity and indication data (Ochoa et al., 2023) | L2, L3 |
| **ChEMBL** | Bioactivity and molecule database (Zdrazil et al., 2024) | L3, via Open Targets |
| **HPO** | Human Phenotype Ontology — standard coded terms for clinical features (Gargano et al., 2024) | L2 channel D |
| **Monarch** | Knowledge graph of disease–phenotype associations (Putman et al., 2024) | L2 channel D |
| **Mondo** | Disease ontology giving diseases stable identifiers | L2 channel D |
| **LINCS** | Perturbation signature library; **GEO** accessions GSE106127 and GSE70138 (Broad Institute & NIH LINCS Program, 2017a, 2017b) | L2 channel C |
| **GEO** | Gene Expression Omnibus — NCBI's archive where the LINCS data lives | L2 channel C |
| **Europe PMC** | Literature search service used to verify each cited claim (Ferguson et al., 2021) | L2 channel E |
| **Crossref** | DOI resolution and retraction checking (Hendricks et al., 2020) | Reference verification |
| **openFDA** | US drug labels and marketing status, US-Government public domain (U.S. Food and Drug Administration, 2026) | L3, L4 |

---

## 6. Identifiers and crosswalks

| Term | Stands for | What it means here |
|---|---|---|
| **RxCUI** | RxNorm Concept Unique Identifier | The identity backbone drugs are harmonised onto, so two names for one drug become one row |
| **RxNorm** | — | The US normalised naming system for clinical drugs |
| **UNII** | Unique Ingredient Identifier | FDA's identifier for a substance |
| **NDC** | National Drug Code | The identifier for a marketed drug product |
| **ATC** | Anatomical Therapeutic Chemical | WHO's drug classification hierarchy |
| **UniChem** | — | EMBL-EBI's crosswalk service between chemical identifier systems |
| **DOI** | Digital Object Identifier | The permanent handle for a paper. Every one cited here was resolved against a real record |
| **PMID / PMCID** | PubMed identifier / PubMed Central identifier | Literature identifiers |
| **SHA-256** | — | The checksum recorded beside every downloaded file, so a re-run can prove it used the same bytes |

---

## 7. Licences

| Term | Stands for | What it means here |
|---|---|---|
| **CC0** | Creative Commons Zero | No rights reserved. Freely redistributable |
| **CC BY** | Creative Commons Attribution | Redistributable with credit |
| **CC BY-SA** | — Attribution-ShareAlike | Derivatives must carry the same licence. **Segregated**, never bundled into this project's output |
| **CC BY-NC** | — Attribution-NonCommercial | No commercial use. Also segregated |
| **CC BY-NC-ND** | — NonCommercial-NoDerivatives | No commercial use and no derivatives |
| **Public domain / US-Government work** | — | Not copyrightable. Why openFDA label sentences can be quoted in the exclusions |
| **Segregation** | — | Keeping licence-restricted data in a separate cache that never reaches a redistributed output. Only identifiers and derived scores cross that boundary |

---

## 8. Drug safety and regulatory

| Term | Stands for | What it means here |
|---|---|---|
| **Repurposing** | — | Finding a new use for an already-approved drug, rather than developing a new one |
| **Chemoprevention** | — | Using a drug to reduce cancer risk. Here it means **secondary prevention** specifically |
| **Secondary prevention** | — | Reducing recurrence or a second cancer in someone whose risk is already raised. **Not** preventing a first cancer, and **not** treating an existing one. This project's fixed endpoint (gate G2) |
| **Genotoxic** | — | Damages DNA. The opposite of the goal in a cancer-predisposition syndrome |
| **Mutagenic** | — | Causes changes in DNA sequence |
| **Clastogenic** | — | Breaks chromosomes |
| **Aneugenic** | — | Causes the wrong number of chromosomes — **the very thing this disease already does** |
| **Carcinogenicity study** | — | The animal study on a label reporting whether a drug increased tumours |
| **SPL** | Structured Product Labeling | The standard format of a US drug label, divided into named sections |
| **Boxed warning** | — | The most serious warning a US label carries |
| **Paediatric use not established** | — | The label does not support use in children. An exclusion here, not a down-rank |
| **Hard gate** | — | A rule that removes a candidate outright. Safety is never a weight to be traded off |

---

## 9. Governance and project shorthand

| Term | Stands for | What it means here |
|---|---|---|
| **D1 – D21** | Decision *n* | Numbered entries in [`mngmt/decisions.md`](../mngmt/decisions.md), each recording what was decided, what forced it, and what it obliged downstream |
| **G1 / G2 / G3** | Gate *n* | Decision points a person had to pass before the work could continue: G1 fixed the causal gene, G2 the therapeutic endpoint, G3 confirmed enough channels were producing |
| **Custody location** | — | A place a copy of the data exists, tracked in [`data_custody.md`](data_custody.md) so deletion can be attested against a list |
| **Purge** | — | The committed deletion of all held data, with an attestation |
| **Provenance** | — | The record of where a result came from: the artifact, the seed, the config and the source versions |
| **Reproducibility** | — | Re-running a stage on unchanged inputs producing byte-identical output |
| **Hypothesis generation** | — | What this project produces. **Not** a diagnosis, a recommendation, or evidence any drug works |
| **n = 1** | — | The study population is one patient. Nothing here is statistically generalisable |
| **LLM** | Large language model | The reasoning step's local, open-weights model, run over loopback so nothing patient-derived leaves the machine |
| **Loopback** | — | A network address that never leaves the machine. Enforced by a guard that refuses any other endpoint |
| **MCP** | Model Context Protocol | A way of connecting external tools to an AI assistant. Deliberately **not** connected to this repository |

---

## References

Full APA 7th entries for every citation above, with the Crossref verification log and the
date each source was checked: [`references.md`](references.md). Machine-readable source:
[`references.bib`](references.bib).

---

*Hypothesis generation only. Not medical advice, not a clinical recommendation, and not a
claim that any drug named here is safe or effective for any person.*
