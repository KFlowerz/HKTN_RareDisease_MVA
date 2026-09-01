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
