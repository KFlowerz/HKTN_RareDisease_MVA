# References

APA 7th edition. Rendered from [references.bib](references.bib), which is the machine-readable
source of truth — add entries there, then regenerate this file.

**Workflow:** Mendeley → export BibTeX by hand → `references.bib` → regenerate this file → verify →
commit. The repo is deliberately **not** connected to Mendeley (no MCP, no API credentials); see
[CLAUDE.md § Reference manager](../CLAUDE.md#reference-manager-manual-export).
Verification needs no account — Crossref (`https://api.crossref.org/works/<DOI>`) resolves DOIs and
reports retractions.

Every entry below has had its DOI or PMID **resolved against a real record** and its retraction
status checked. See [CLAUDE.md § Evidence and citations](../CLAUDE.md#evidence-and-citations) for the
policy, including the two admissible evidence types (pipeline data result, or cited source) and the
in-code citation format.

---

## Reference list

Conlin, L. K., Thiel, B. D., Bonnemann, C. G., Medne, L., Ernst, L. M., Zackai, E. H., Deardorff,
M. A., Krantz, I. D., Hakonarson, H., & Spinner, N. B. (2010). Mechanisms of mosaicism, chimerism
and uniparental disomy identified by single nucleotide polymorphism array analysis. *Human
Molecular Genetics, 19*(7), 1263–1275. https://doi.org/10.1093/hmg/ddq003

Hanks, S., Coleman, K., Reid, S., Plaja, A., Firth, H., FitzPatrick, D., Kidd, A., Méhes, K., Nash,
R., Robin, N., Shannon, N., Tolmie, J., Swansbury, J., Irrthum, A., Douglas, J., & Rahman, N.
(2004). Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B.
*Nature Genetics, 36*(11), 1159–1161. https://doi.org/10.1038/ng1449

Loh, P.-R., Genovese, G., Handsaker, R. E., Finucane, H. K., Reshef, Y. A., Palamara, P. F.,
Birmann, B. M., Talkowski, M. E., Bakhoum, S. F., McCarroll, S. A., & Price, A. L. (2018). Insights
into clonal haematopoiesis from 8,342 mosaic chromosomal alterations. *Nature, 559*(7714), 350–355.
https://doi.org/10.1038/s41586-018-0321-x

---

## Verification log

| Key | Identifier | Resolved | Retraction/correction | Checked |
|---|---|---|---|---|
| `hanks2004` | doi:10.1038/ng1449 | ✅ Crossref | None on record | 2026-09-01 |
| `conlin2010` | doi:10.1093/hmg/ddq003 | ✅ Crossref | None on record | 2026-09-02 |
| `loh2018` | doi:10.1038/s41586-018-0321-x | ✅ Crossref | None on record | 2026-09-02 |

Re-run the retraction check before submission — status can change between drafting and publication.
