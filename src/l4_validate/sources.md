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
[src/l2_channels/enrichment.py](../l2_channels/enrichment.py). What reaches
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

---

## Table 6 — phenotype matching and indications (L2 channel D)

Channel D matches the subject's HPO-coded phenotype against diseases' annotated features,
then takes the approved drugs indicated for the closest diseases or for the subject's
features. Terms read 2026-09-17, each from the project's own source (web pages that are
script-only apps were read from their source repositories).

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Phenotype vocabulary and hierarchy | HPO `hp.obo`, release 2026-09-01 | **HPO licence (custom)**: free to use on condition of citation, showing the file's date/version wherever displayed publicly, and no alteration of content or logical relationships; services must acknowledge HPO use | ❌ Not CC — but only ids, counts and an acknowledgement are output | `reference_dir`, unaltered |
| Disease → feature annotations | Monarch KG 2026-09-02, `disease_phenotype.all.tsv.gz` (HPO annotations; primary sources include OMIM and Orphanet) | Annotations under the HPO licence; Monarch's own data recommendation is CC0/CC BY with per-source terms prevailing | ❌ Not redistributed | `reference_dir`; matched locally |
| Disease names and ids | Mondo, as labelled in Monarch's table | CC BY 4.0 | ✅ Yes | patient-derived `evidence.json` only (see below) |
| Orphanet-sourced annotations | Orphadata **Science** datasets | CC BY 4.0 — the rest of Orphadata needs written consent from INSERM | ✅ Yes (Science datasets only) | inside the Monarch table |
| OMIM-sourced annotations | OMIM | **Unread**: OMIM's use agreement returned HTTP 403 on 2026-09-17 | ❓ Treated as restricted | used locally only; no OMIM id, name or text is output |
| Drug → approved indication | Open Targets Platform 26.06 `clinical_indication` | as Table 5 — ChEMBL-derived | ❌ No — treated as SA | `enrichment_dir`, never redistributed |

**The HPO licence is the binding one, and it is not a Creative Commons licence.** Its three
conditions are met by construction: the ontology is cached and parsed unaltered; its
version is written into `channel.json`'s `attribution` field and must accompany any public
display of HPO-derived results; and the report cites Gargano et al. (2024) and acknowledges
HPO. What leaves the pipeline is ranked drug identities and counts, not HPO content.

**OMIM stays unread, so its content stays in.** Most rare-disease annotations originate
from OMIM (8,478 of the 12,882 diseases in the HPO annotation file). Its terms could not be
read, so nothing OMIM-authored is output: diseases are keyed and named by Mondo. Revisit
this row when the agreement can be read.

**What reaches `candidates.tsv`.** ChEMBL id, drug name, type, clinical stage (the Table 5
whitelist), the pipeline's own score, the route, support *counts*, and
`endpoint=symptomatic`. No HPO id, no disease, no indication pair.

**`evidence.json` is patient-derived, not just licence-restricted.** Which diseases resemble
this child, and which of the child's features an indication matched, are statements about
the child (`COMPLIANCE.md`). The file is marked `patient_derived: true` and
`redistributable: false` and stays under the gitignored `results_dir`.

**No warranty, no clinical use.** Monarch's terms of use state its tools "should not be
used for direct diagnostic use or medical decision-making", provide the data "as is"
without warranty, and place "total and exclusive responsibility and risk" on the user.
Channel D's caveats carry this, and the report must.

**Why local, not Monarch's API.** Querying Monarch's service would send the subject's
feature combination to a third party, and Monarch's website states it runs HotJar and
Google Analytics. The files are downloaded whole and matched on this machine; see
[decision D10](../../mngmt/decisions.md).

## Table 7 — literature verification (L2 channel E), added 2026-09-17

Channel E ranks a curated compound list by the literature Europe PMC returns for each, and
resolves every citation before reporting it. These are the pipeline's **only** outbound
queries that are not whole-file downloads; what makes that acceptable is that the query
carries public vocabulary and nothing derived from the subject, enforced in
`src/l2_channels/literature.py` (`refuse_private`) and asserted in
`tests/test_l2_literature.py`. Terms read 2026-09-17.

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Record retrieval, identifiers, publication types, MeSH, retraction status | Europe PMC RESTful Web Service | **Europe PMC copyright terms**: the RESTful service is one of the four routes Europe PMC names as "the only services that may be used for automated downloading"; scraping the website is prohibited. Plus **EMBL-EBI terms of use** (EMBL-EBI hosts the service): free, no key, attribution expected, provided "AS IS" without warranties | ✅ For bibliographic metadata | responses cached under `reference_dir/literature/` |
| Article abstracts (matched locally against fixed patterns) | Europe PMC | Publisher or author copyright; Europe PMC redistributes under each publisher's own terms | ❌ No | **cache only** — never written to an output |
| DOI resolution and a second retraction check (`update-to`) | Crossref REST API | "Almost all of the metadata we hold is reusable without restriction… considered to be 'facts' which are not copyrightable and are thus in the public domain (CC0)"; Crossref-generated data released as public domain | ✅ Yes | responses cached under `reference_dir/literature/` |
| Compound identity (ChEMBL id, clinical stage) | Open Targets Platform 26.06 `drug_molecule` | as Table 5 — ChEMBL-derived | ❌ No — treated as SA | `enrichment_dir`, never redistributed |

**The REST API is the sanctioned route, and scraping the site is not.** Europe PMC's
copyright page is explicit: "Crawlers and other automated processes may **NOT** be used to
systematically retrieve batches of articles from the Europe PMC web site. Bulk downloading
of articles from the main Europe PMC site, in any way, is prohibited because of copyright
restrictions." It then names OAI, RESTful, SOAP and bulk download as the sanctioned
alternatives: "These are the only services that may be used for automated downloading of
articles in Europe PMC." This channel uses the RESTful service and never fetches a page
from the website, which is what that paragraph requires.

**Article content is publisher copyright, and per-article licences differ.** "All of the
material available through the Europe PMC site is provided by the respective publishers or
authors. Almost all of it is protected by U.K. and/or foreign copyright laws, even though
Europe PMC provides free access to it." Material reached through the sanctioned services is
"still protected by copyright, but… distributed under a Creative Commons or similar
license… The license terms are not identical for all the articles." Since the terms vary
per article and this pipeline does not read each article's licence statement, no output
reproduces abstract or full text at all. `tests/test_l2_channel_e.py` asserts that abstract
text reaches no file the channel writes, and the first real run scanned 56 cached abstracts
against every output and found none.

**How that page was read, because the live one is unreachable.** Every request to
`europepmc.org/Copyright` from this machine is answered by a Cloudflare challenge (HTTP 403,
to `curl` and to a browser user agent alike), and the Internet Archive's own 2026 captures
recorded the same 403. The wording above is quoted from the Archive's snapshot of
**2025-06-05**, the most recent successful capture:
`https://web.archive.org/web/20250605155112/https://europepmc.org/Copyright`. An archived
page can lag the live one, so this row names its date; re-read it before submission if the
live page becomes reachable.

**What reaches `candidates.tsv`.** The Table 5 identity whitelist (ChEMBL id, drug name,
type, clinical stage), the channel's own score, the curated `axis`, `target`, `direction`
and `caution` columns from `src/l2_channels/aneuploidy_prior.tsv` (this project's own
text), record *counts* per evidence grade, and `endpoint=chemoprevention`.

**What reaches `references.json`.** Bibliographic metadata only — PMID, DOI, title,
journal, year, publication types, the assigned grade, and retraction status. This is the
redistributable layer, and it is what lets a reader check every claim the channel makes.

**Rate limits, observed rather than assumed.** Crossref's anonymous pool returned
`x-rate-limit-limit: 5` per `1s`. One throttle covers both endpoints and is set from the
tighter of the two limits, a little below it.

**No warranty, no clinical use.** EMBL-EBI provides its resources "AS IS" without
warranties of any kind and does not guarantee the accuracy of the data nor its suitability
for any purpose. Nothing this channel produces is a clinical recommendation; see
[decision D11](../../mngmt/decisions.md).

## Table 8 — drug identity harmonisation (L3), added 2026-09-18

L3 must give every candidate one identity before anything is aggregated, or a molecule and
its own salt form compete for rank and the cross-channel convergence signal is destroyed.
Two steps, and only the second produces anything publishable. Terms read 2026-09-18.

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Which molecules are formulations of which (`parentId`) | Open Targets Platform 26.06 `drug_molecule` | as Table 5 — ChEMBL-derived | ❌ No — treated as SA | `enrichment_dir`, never redistributed |
| ChEMBL id → UNII | UniChem whole-source mapping `src1src14` (EMBL-EBI) | **No licence statement readable** — see below | ❌ Treated as ChEMBL-derived | `enrichment_dir`, never redistributed |
| UNII → RxCUI, plus pharmacologic classes | openFDA NDC Directory + Drugs@FDA bulk | US-Government work, **public domain** | ✅ Yes | `reference_dir` |

**Why not RxNav, which would answer this in one call per drug.** Channel D's candidates are
derived from the subject's phenotype, so the *set* of drugs asked about is itself a weak
statement about the child. A bulk file is asked about nothing. This is the same rule the
rest of the pipeline follows, applied to the one case where the convenient route is a
per-item query — and it is why the crosswalk is assembled from two bulk files instead.

**Why not RxNorm's own full release.** It requires a UMLS licence key. A judge rebuilding
from this repository alone could not run it, which the reproducibility requirement in
`CLAUDE.md` rules out.

**UniChem states no licence this project could read**, on 2026-09-18: its FAQ is a
script-only page and the FTP `README` carries none. The mapping is therefore treated as
ChEMBL-derived — UniChem is hosted alongside ChEMBL at EMBL-EBI — and kept in the
non-redistributed zone with everything else under decision D7. This is the conservative
reading, recorded rather than assumed, the same way the OMIM row in Table 6 and the
Europe PMC row in Table 7 are.

**What reaches a redistributed output.** The Table 2 identity whitelist (ChEMBL id, drug
name, type, clinical stage) plus `rxcui`, `unii`, `pharm_class_epc` and `pharm_class_moa`,
which are openFDA's and public domain. The ChEMBL→UNII mapping itself and the `parentId`
relation never leave; they are used to compute identity and then discarded.

**Coverage is partial, and it is reported per row.** Of 2,575 ChEMBL ids across the three
channels, 865 resolved to an RxCUI through a registered substance (UNII), 251 through a
name match, and 854 not at all — roughly 45% unresolved. The gap is not an error: openFDA
covers products marketed in the United States, and many ChEMBL-approved drugs are approved
elsewhere. `rxcui_resolution` records which route produced each row, so a name match is
never mistaken for a registry match, and an unresolved candidate is kept and labelled
rather than dropped. Dropping them would silently restrict the pipeline to the US market
while looking like a technical detail.

**openFDA's own warning travels with the data.** It states its data is not for clinical use
and may be incomplete or inaccurate. Nothing in this pipeline is a clinical recommendation;
see [decision D14](../../mngmt/decisions.md).

## Table 9 — safety triage and benchmark (L4), added 2026-09-18

Everything L4 decides on is US-Government public domain, so every verdict can quote the
sentence it rested on and that quotation is redistributable. Terms as Table 1; the bulk
releases were read 2026-09-18.

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Carcinogenesis and mutagenesis text (SPL 13.1, with 13 as fallback) | openFDA `drug/label` bulk, export 2026-09-18 (262,883 records, 14 parts) | US-Gov public domain | ✅ Yes | `reference_dir` |
| Paediatric use (SPL 8.4), boxed warning, warnings, contraindications | openFDA `drug/label` bulk | US-Gov public domain | ✅ Yes | `reference_dir` |
| Current marketing status per product | openFDA `drug/drugsfda` bulk, export 2026-09-18 (29,335 records) | US-Gov public domain | ✅ Yes | `reference_dir` |
| Established Pharmacologic Class / mechanism of action | openFDA `openfda` block, attached at L3 | US-Gov public domain | ✅ Yes | carried on the L3 row |
| Benchmark positive set | `src/l2_channels/aneuploidy_prior.tsv` | This project's own curation | ✅ Yes | committed |

**Field names are verified against the release, not assumed.** The carcinogenesis section
is `carcinogenesis_and_mutagenesis_and_impairment_of_fertility`. An earlier guess at
`carcinogenesis_and_mutagenesis_of_fertility` matched nothing in 262,883 records and
silently excluded **every** candidate as `insufficient_evidence` — a wrong field name that
looked exactly like a strict safety filter working. Any new field read here is confirmed
against the corpus first.

**Snippets are quotable because the source is public domain.** Each verdict carries the
sentence it rested on, trimmed but never paraphrased, plus the SPL section number a reader
would cite. That is what makes the exclusions table auditable rather than a list of
assertions, and it is only possible because openFDA is a US-Government work.

**Absence of a field is never evidence of safety.** 1,231 of 1,880 exclusions are
`insufficient_evidence`, overwhelmingly because openFDA carries no label for that
candidate — it describes drugs marketed in the United States, and the channels nominate
from a wider pool. Reported as a coverage fact, not as a safety finding.

**openFDA's own warning travels with every verdict:** its data is not for clinical use and
may be incomplete or inaccurate. Nothing this layer produces is a clinical recommendation.

**No interaction claim is made.** The NLM RxNav drug-interaction API was discontinued on
2024-01-02 and no structured, licence-clean DDI source with severity exists (see *Currency
gotchas* above). L4 therefore makes no interaction claim rather than an unsourced one, and
says so in its caveats.

**The benchmark's positive set is this project's own file**, which is why channel E cannot
be scored on it — see [decision D15](../../mngmt/decisions.md). The headline is channel B's
recovery, because channel B never saw the set.

---

## Table 10 — transcriptional signatures (L2 channel C), added 2026-09-21

Channel C scores compounds for *reversing* a proxy disease signature. Both releases are
NCBI GEO supplementary files from the NIH LINCS Program, downloaded whole and matched
locally — never queried per compound. See [decision D19](../../mngmt/decisions.md) for why
the channel ships disabled.

| Role | Source | License | CC-BY-4.0 compatible? | Where it lives |
|---|---|---|---|---|
| Proxy disease signature: consensus shRNA knockdown per gene, 33,839 × 978 | GEO `GSE106127` (`CGS` matrix + `CGS_meta`), retrieved 2026-09-21 | See below | ❌ Not redistributed | `reference_dir/lincs` |
| Compound signatures to score: Level 5 MODZ, 118,050 × 12,328 | GEO `GSE70138` (LINCS Phase II, build 2017-03-06) | See below | ❌ Not redistributed | `reference_dir/lincs` |
| Landmark gene identity (`pr_is_lm`), Entrez ids | GEO `GSE70138` and `GSE106127` `gene_info` | See below | Identifiers only | `reference_dir/lincs` |
| Drug identity for the LINCS `pert_iname` → ChEMBL crosswalk | Open Targets `drug_molecule` | as Table 5 — ChEMBL-derived | ❌ No — treated as SA | `enrichment_dir` |

**The licence could not be resolved, so it is treated as restricted.** On 2026-09-21
`lincsproject.org` returned HTTP 404 at both its root and its data-release-policy path,
and `clue.io/terms` served a script-rendered glossary carrying no licence statement. LINCS
data is widely described as CC BY 4.0, but this project does not record a licence it could
not read. The conservative handling costs nothing here: **nothing from these files is
redistributed**. What channel C would publish is a derived score per drug plus identifiers
and the caveat string — no expression value, no signature, no matrix row. The same posture
as Table 6's HPO and Monarch rows.

**Cite `subramanian2017` for the method and the data**, `lamb2006` for the connectivity
concept, and `subramanian2005` for the running-sum enrichment score underneath it.
`lincs_gse106127` and `lincs_gse70138` carry the dataset records.

**These are large, and they are cached rather than committed.** The expanded `GSE70138`
matrix is 5.8 GB. [src/l2_channels/lincs.py](../l2_channels/lincs.py) downloads the GEO
`.gz`, expands it once, records the compressed file's SHA-256 and then removes it —
keeping both would double a 6 GB cache for nothing. Checksums of the expanded files are
memoised in a sidecar, so a 5.8 GB file is not rehashed on every run. As with every other
cached release these live outside the repository and outside the custody root, and must
never enter the deletion register in
[docs/data_custody.md](../../docs/data_custody.md).

**The GCTX axes are transposed relative to their names.** `0/META/ROW/id` names the genes
and `0/META/COL/id` names the signatures, but the HDF5 matrix is stored column-major, so
`matrix[i, j]` reads as `[signature i, gene j]`. Reading it the other way yields a matrix
of the right dtype and the wrong meaning, which nothing downstream would catch;
`lincs.axes()` asserts the orientation on every open.
