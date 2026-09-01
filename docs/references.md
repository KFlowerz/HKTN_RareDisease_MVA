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

Hanks, S., Coleman, K., Reid, S., Plaja, A., Firth, H., FitzPatrick, D., Kidd, A., Méhes, K., Nash,
R., Robin, N., Shannon, N., Tolmie, J., Swansbury, J., Irrthum, A., Douglas, J., & Rahman, N.
(2004). Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B.
*Nature Genetics, 36*(11), 1159–1161. https://doi.org/10.1038/ng1449

---

## Verification log

| Key | Identifier | Resolved | Retraction/correction | Checked |
|---|---|---|---|---|
| `hanks2004` | doi:10.1038/ng1449 | ✅ Crossref | None on record | 2026-09-01 |

Re-run the retraction check before submission — status can change between drafting and publication.
