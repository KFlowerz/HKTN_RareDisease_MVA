# Drug annotation sources — catalog and licensing

How every candidate gets annotated in L4, and which fields may appear in a
**CC-BY-4.0 redistributed output** versus which must stay in the non-redistributed
`enrichment/` zone.

## The rule

1. Resolve every candidate to an **RxCUI** (RxNorm) — the identity backbone.
2. Cross-map outward to DrugBank / ChEMBL / UNII / ATC / PubChem / NDC via **UniChem**.
3. Build the redistributed layer **only** from license-clean sources (table 1).
4. Everything else lives in `enrichment/`, joined by RxCUI, and is **never** bundled
   into a CC-BY-4.0 output (table 2).

Cross-mapped *identifiers* from an NC/SA source may be used to perform the join. The
*annotation content* from those sources may not be redistributed.

---

## Table 1 — CC-BY-4.0-clean core (use freely)

| Field | Vocabulary | Source | License | CC-BY-4.0 compatible? |
|---|---|---|---|---|
| Drug identity (ingredient, brand, salt form) | RxNorm RxCUI, TTY | RxNorm | US-Gov public domain | ✅ Yes |
| Therapeutic class | RxClass (ATC-US, VA, MeSH) | RxClass | US-Gov public domain | ✅ Yes |
| Mechanism of action (MoA) | MED-RT | MED-RT | US-Gov public domain | ✅ Yes |
| Physiologic effect / PK class | MED-RT | MED-RT | US-Gov public domain | ✅ Yes |
| Pharmacologic class (EPC) | MeSH / NDF-RT descendants | MeSH | US-Gov public domain | ✅ Yes |
| Indications | *prose* — Section 1 | openFDA drug label | CC0 | ✅ Yes |
| Contraindications | *prose* — Section 4 | openFDA drug label | CC0 | ✅ Yes |
| Boxed warning | *prose* — Section 5 | openFDA drug label | CC0 | ✅ Yes |
| Warnings & precautions | *prose* — Section 5/6 | openFDA drug label | CC0 | ✅ Yes |
| Adverse reactions | *prose* — Section 6 | openFDA drug label | CC0 | ✅ Yes |
| Carcinogenesis / mutagenesis | *prose* — Section 13.1 | openFDA drug label | CC0 | ✅ Yes |
| Use in specific populations (incl. pediatric, pregnancy) | *prose* — Section 8 | openFDA drug label | CC0 | ✅ Yes |
| Dosage & administration | *prose* — Section 2 | openFDA drug label | CC0 | ✅ Yes |
| Dosage forms & strengths | *prose* — Section 3 | openFDA drug label | CC0 | ✅ Yes |
| Approval status, application number | — | openFDA drugsfda | CC0 | ✅ Yes |
| Post-market adverse event signals | MedDRA PT (as reported) | openFDA FAERS | CC0 | ✅ Yes |
| Active moiety identity | UNII | UNII / GSRS | US-Gov public domain | ✅ Yes |
| Marketed product / packaging | NDC | NDC directory | US-Gov public domain | ✅ Yes |

**Most of this is prose, not structured fields.** Indications, contraindications,
warnings, pediatric use, and the carcinogenesis section arrive as free text that has to
be *extracted* before L4 can gate on it — that extraction is the Claude layer's job
([../l3_integrate/claude_reasoning.py](../l3_integrate/claude_reasoning.py)). Preserve a
quotable snippet plus its section number alongside every extracted value, so each safety
verdict is traceable back to label text.

---

## Table 2 — Must segregate (NC / ShareAlike)

Keep in `enrichment/`, joined by RxCUI. **Never** bundle into a redistributed output.

| Field | Vocabulary | Source | License | CC-BY-4.0 compatible? |
|---|---|---|---|---|
| Curated side effects | MedDRA / UMLS | SIDER | CC BY-NC-SA 4.0 | ❌ No — NC **and** SA |
| Drug targets, pathways, pharmacology | DrugBank ID | DrugBank | CC BY-NC 4.0 | ❌ No — NC |
| Drug–drug interactions with severity | DDInter | DDInter | CC BY-NC-SA 4.0 | ❌ No — NC **and** SA |
| Bioactivity, assays, target affinities | ChEMBL ID | ChEMBL | CC BY-SA 3.0 | ❌ No — SA (viral) |
| Pharmacogenomics, clinical annotations | PharmGKB | PharmGKB | CC BY-SA 4.0 | ❌ No — SA (viral) |
| Full ATC hierarchy (bulk) | ATC | WHO ATC/DDD | Proprietary / licensed | ❌ No |

**Why ShareAlike matters as much as NonCommercial here.** SA licenses are viral: mixing
ChEMBL or PharmGKB content into a derived table would force the *entire* output under
CC BY-SA, breaking the CC-BY-4.0 commitment in [../../COMPLIANCE.md](../../COMPLIANCE.md).
NC is the more obvious constraint; SA is the one that quietly relicenses everything.

**Substitutes in the clean core:** for drug targets use RxClass MoA + MED-RT instead of
DrugBank; for the ATC hierarchy use RxClass's ATC-US view rather than the WHO bulk file;
for adverse events use FAERS (CC0) rather than SIDER.

---

## Currency gotchas — encode these, don't rediscover them

**1. FDA pregnancy letter categories are retired.**
The A/B/C/D/X system was phased out under the Pregnancy and Lactation Labeling Rule
(PLLR). Current labels carry a **narrative Section 8** (Use in Specific Populations)
instead. Any code or prompt looking up a letter grade will find nothing on a modern
label, and any letter found on an old label is not authoritative. Parse the narrative.

**2. The RxNav drug-interaction API was discontinued 2024-01-02.**
NLM retired both the RxNorm and NDF-RT interaction endpoints. **There is no structured,
CC-BY-clean, DDI-with-severity source.** The options are:

- license a commercial source (First Databank, Cerner Multum, Micromedex), or
- run NLP over SPL Section 7 (Drug Interactions) prose from openFDA — clean, but
  unstructured and without graded severity, or
- use DDInter / DrugBank in the `enrichment/` zone only, never redistributed.

Whichever is chosen, **label DDI output with its provenance and its severity semantics**.
A DDI table without a severity scale attached is not actionable, and one whose scale is
undocumented is worse than none.

---

## Practical notes

- **UniChem is the cross-mapper**, and it maps *identifiers*, not licenses. An identifier
  obtained through UniChem does not make the source's content redistributable.
- **Salt forms, racemates, and combination products are distinct entities.** Normalize to
  ingredient level (TTY `IN`/`PIN`) for consensus counting, but keep the salt form
  visible — dose and formulation depend on it.
- **Record a retrieval date and source version for every annotation.** openFDA labels
  change; a safety verdict is only reproducible if the label snapshot behind it is
  identified.
- **openFDA is not exhaustive.** Absence of a label field is not evidence of safety —
  see the fail-closed rule in [safety_triage.py](safety_triage.py).

---

## Table 3 — network and pathway sources (L1)

The disease module is built from two sources, both clean for a redistributed CC-BY-4.0
output. Licences verified 2026-09-16 from each project's own licence page.

| Role | Source | License | CC-BY-4.0 compatible? | Obligation |
|---|---|---|---|---|
| Protein–protein interactions, confidence-scored | STRING v12.0 | CC BY 4.0 | ✅ Yes | Attribution, and say if the data were changed |
| Physical-subnetwork interactions (complex proxy) | STRING v12.0 | CC BY 4.0 | ✅ Yes | As above |
| Protein ↔ symbol ↔ UniProt mapping | STRING v12.0 | CC BY 4.0 | ✅ Yes | As above |
| Pathway membership (upstream/downstream split) | Reactome | CC0 1.0 | ✅ Yes | None; attribution encouraged |

STRING asks that users "provide appropriate credit — and inform users of any changes or
additions that you might have made to the data". That attribution, plus the retrieval
date, byte size and SHA-256 of every file, travels in
[module.json](../../results/l1_target/module.json) (gitignored, regenerated by L1), so a
module can be traced to the exact bytes behind it.

**Deliberately not used: CORUM.** The L1 scaffold suggested it for complex membership, but
its terms are not clean for redistribution, and a complex database is not worth the
licensing exposure when STRING's physical subnetwork approximates the same thing under
CC BY 4.0. L1's seed set is that approximation, and says so in its own provenance.

**The rule that applies here too:** an identifier obtained from a source does not carry
that source's license to its content. If a later layer needs curated complex membership,
it belongs in the non-redistributed enrichment zone, joined by identifier — never bundled
into an output.

---

## Table 4 — variant interpretation source (L0)

L0 cross-references the subject's panel alleles against ClinVar, so gate G1 is decided
against what has already been submitted about those variants and not against consequence
prediction alone. Terms read 2026-09-16 from NCBI's
[data usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) and ClinVar's
[maintenance and use](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/) page.

| Role | Source | License | CC-BY-4.0 compatible? | Obligation |
|---|---|---|---|---|
| Submitted clinical interpretations, review status, conditions | ClinVar (GRCh38 VCF) | US-Gov work; no restriction asserted by NCBI | ✅ Yes | Attribution requested — cite `landrum2018` |

**What NCBI actually says, including the part that is not a blanket grant.** "NCBI itself
places no restrictions on the use or distribution of the data contained therein. Nor do we
accept data when the submitter has requested restrictions on reuse or redistribution." It
then adds that individual submitters may claim rights in their own submitted data, that no
rights transfer to NCBI, and that NCBI therefore "cannot provide comment or unrestricted
permission" on redistribution. ClinVar's own page asks that redistributors "provide
attribution to ClinVar as a data source", citing a ClinVar publication (PMID 29165669 =
`landrum2018`).

Practically, that is clean for this project: what L0 redistributes is a classification
label, a review status and a variation id per allele — facts about a public record, carried
with their attribution — not a submitter's curated prose. Nothing here is bulk-republished.

**The whole release is downloaded, never queried per variant.** A ClinVar API lookup would
send this child's coordinates to an NCBI server, which [COMPLIANCE.md](../../COMPLIANCE.md)
prohibits — the same reason annotation runs on a local snpEff rather than the VEP REST API.
The release's URL, fileDate, byte size and SHA-256 travel in
[clinvar_crossref.json](../../results/l0_genomics/clinvar_crossref.json) (gitignored,
regenerated by L0).

**What comes back is not redistributable as a classification.** ClinVar aggregates
submissions of varying quality; the archive's own review status travels with every record
L0 reports, and nothing in the pipeline converts a match into a pathogenicity call. An
output that quoted a ClinVar significance without its review status would misrepresent a
one-submitter opinion as the archive's position.

### Table 4b — population frequency source (L0), added 2026-09-17

Whether a protein-altering second hit is simply common was the largest remaining gap in L0's
evidence; decision D8 closes it with gnomAD. Terms read 2026-09-17 from gnomAD's
[policies](https://gnomad.broadinstitute.org/policies) page, verbatim from its source,
[`browser/about/policies/terms.md`](https://github.com/broadinstitute/gnomad-browser/blob/main/browser/about/policies/terms.md)
(the rendered page is a script-only app).

| Role | Source | License | CC-BY-4.0 compatible? | Obligation |
|---|---|---|---|---|
| Allele count, allele number, frequency, homozygote count, highest group frequency | gnomAD v4.1.1 exomes + genomes sites VCFs (`AC`, `AN`, `AF`, `nhomalt`, `AF_grpmax` only) | CC0 1.0 | ✅ Yes | None legally; attribution requested — cite `chen2024` and `gnomad_v411` |
| SpliceAI scores shipped in the same files | gnomAD (computed by Illumina) | **CC BY-NC 4.0** | ❌ No | **Not read.** `src/l0_genomics/gnomad.py` names every field it parses, and a test forbids the rest |

**What gnomAD says.** "The primary data from the gnomAD exomes and genomes are available free
of restrictions under the Creative Commons Zero Public Domain Dedication"; SpliceAI annotations
are "provided with permission under a CC BY NC 4.0 license for academic and non-commercial
use". A CC0 file can carry non-CC0 columns, which is why the licence is enforced per field
rather than per file. The citation policy asks that users cite "the gnomAD flagship paper",
doi:10.1038/s41586-023-06045-0 (`chen2024`), and the terms state that "all users of gnomAD
data agree to not attempt to reidentify participants" — an obligation this project takes on
by reading the data at all.

**Read by gene span, not downloaded whole.** Chromosome 15 alone is 24 GB across the two
data sets, so the ClinVar rule (whole release, matched locally) is replaced by the nearest
thing that leaks nothing: htslib reads the byte ranges covering each panel gene's whole padded
span. The request is the same for every proband run through the panel and carries no subject
coordinate. See decision D8 and the test that asserts it.

**What comes back is a frequency, not a verdict.** Rarity is necessary but not sufficient for
pathogenicity; nothing in the pipeline applies ACMG/AMP PM2 or converts a frequency into a
classification. The caveats travel in
[gnomad_frequencies.json](../../results/l0_genomics/gnomad_frequencies.json) (gitignored,
regenerated by L0).

---

## Table 5 — drug–target sets for network proximity (L2 channel B)

Channel B needs, per drug, the human proteins it acts on, as nodes in the interactome.
Terms read 2026-09-16 from each project's own licence page.

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Drug → target protein, action type, mechanism text | Open Targets Platform (`drug_mechanism_of_action`) | Platform marked CC0 1.0; this content derives from **ChEMBL, CC BY-SA 3.0** | ❌ No — treated as SA | `enrichment_dir`, never redistributed |
| Drug identity: ChEMBL id, name, type, clinical stage | Open Targets Platform (`drug_molecule`) | as above | Identifiers and names only | published |
| Ensembl gene id → interactome protein | STRING v12.0 aliases | CC BY 4.0 | ✅ Yes | `reference_dir` |

**Why this is Table 5 and not Table 1.** Open Targets marks its Platform data CC0 and says
downstream users may consume it without restriction. The same licence page lists ChEMBL as
CC BY-SA 3.0 among its sources, and the drug–target content is ChEMBL's. This project does
not rely on one aggregator's relicensing of another's ShareAlike data — see
[decision D7](../../mngmt/decisions.md) — so the content is treated as ShareAlike and kept
in the segregated zone described in Table 2.

**What that means in practice.** The target sets never leave
[src/l2_channels/enrichment.py](enrichment.py). What reaches
`results/l2/channel_b_proximity/candidates.tsv` is the whitelist in
`enrichment.PUBLISHABLE_FIELDS` — ChEMBL id, drug name, drug type, clinical stage, and the
*count* of targets — plus the scores this pipeline computed. Identifiers and INN-style
names are exactly what Table 2 already permits for performing a join;
`tests/test_l2_channel_b.py` asserts the written table carries no target identifier.

**The crosswalk is clean by construction.** Ensembl gene ids are mapped into the
interactome using STRING's own alias table (CC BY 4.0), not the restricted source, so drug
*identity* is the only thing taken from the enrichment zone.

**If the rationale needs to name a target.** Take it from the drug's FDA label
(openFDA, CC0, already in Table 1) rather than from this source. The label states the
molecular target for most approved drugs, and it is quotable in a redistributed output.
