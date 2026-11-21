# Data destruction attestation

**Generated:** 2026-11-21T17:21:15Z  
**Dataset:** `SageBio/mva-hackathon-2026-data` at revision `59e322d27f399006b398d366d33e703e48a29914`  
**Deletion deadline:** 2026-11-23T23:59:00Z

Identified by Hub revision only. No filenames or per-file digests are recorded, by design.

## Locations

| ID | Method | Before | After |
|---|---|---:|---:|
| `C1` | rm -rf, then fstrim on the mount | 13 files, 317,500,000 B | absent |
| `C2` | hf cache rm + hf cache prune, then rm -rf of the root | 4 files, 453,301 B | absent |
| `C3` | rm -rf contents, keep .gitkeep | 242 files, 13,105,072 B | 0 files, 0 B |
| `C4` | Strip outputs; verify no genomic content in tracked .ipynb | 1 files, 0 B | 0 files, 0 B |
| `C7` | Delete the local store at the path above (this tool); delete the conversations in the account, which removes them from history immediately and from back-end storage within 30 days (manual — not performed by this tool) | 1 files, 621,046 B (earlier run: 1 files, 711,986 B) (earlier run: 13 files, 63,625,301 B) | absent |

## Checks

| Check | Result | Detail |
|---|---|---|
| custody_register_complete | **pass** | 9 of 9 register rows reached a terminal state |
| repo_scan_genomic_extensions | **pass** | tests/test_smoke.py::test_no_patient_data_files_present |
| hf_cache_empty | **pass** |  |
| incomplete_blobs_removed | **pass** |  |
| freed_blocks_discarded | **pass** | fstrim -av completed |

## Attestation

I verified each registered location was emptied on the date recorded above.

— Karen Flores, Maintainer

*No filenames, per-file digests, variant coordinates or phenotype content appear in this record, by construction. See `docs/data_custody.md` for the register and the honest limits of any destruction record.*
