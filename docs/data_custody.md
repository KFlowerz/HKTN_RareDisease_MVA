# Data custody and destruction evidence

[COMPLIANCE.md](../COMPLIANCE.md) commits this project to deleting all data within 30 days of
Hackathon close and emailing confirmation to the organizers. This document is how that confirmation
is made *evidenced* rather than merely asserted.

## Why this file exists before the data arrives

A destruction record assembled at deletion time is reconstructed from memory, and memory omits
things — the cache directory nobody set, the intermediate a layer wrote three weeks ago, the
notebook cell that rendered a variant table. **Custody is therefore instrumented at creation time.**
Every location that will hold patient-derived bytes is registered here the moment it is created, and
the purge tool refuses to attest to anything the register does not list.

The register is the contract. If a location is missing from it, the attestation is false — which is
why adding a row is part of creating a location, not a cleanup step afterwards.

---

## What counts as patient-derived

Not just the download. Anything computed from it that retains subject-specific signal:

- variant coordinates, genotypes, and zygosity calls (L0)
- per-chromosome read depth, bin-level coverage, and the aneuploidy burden vector (L0)
- any phenotype narrative or identifier from the gated dataset
- any figure rendered at a resolution that permits re-identification

In a worldwide population of roughly 50 patients, a per-chromosome burden plot or a single exact
variant coordinate can identify on its own. Aggregate and categorical forms are the only publishable
derivatives; see the L5 guardrail in [../src/l5_report/run.py](../src/l5_report/run.py).

**What is not patient-derived, and why it looks like it might be.** The reference cache
(`$MVA_REF_ROOT`, default `~/.cache/mva-track2/reference`) holds whole public releases — STRING,
Reactome, and ClinVar's GRCh38 VCF. It is **not** in the register below and must not enter a purge
attestation: nothing in it derives from the subject, and every file is a complete public download
with no query attached, so not even a gene symbol was disclosed to obtain it. It contains a
`.vcf.gz`, which is why it lives outside the repository — `src/refcache.py` refuses a location
inside it, and `tests/test_smoke.py` still forbids any genomic file in the working tree.

---

## Custody register

Status values: `planned` (registered, not yet created) · `active` (holds bytes) · `purged`
(deleted and verified empty).

| ID | Location | Holds | Purge method | Status |
|---|---|---|---|---|
| `C1` | `$MVA_DATA_ROOT/raw` (WSL ext4, outside the repo) | The gated dataset as downloaded — VCF, tabix index, phenotype `.docx`, README (5 files, 303 MB). FASTQ deliberately not fetched. | `rm -rf`, then `fstrim` on the mount | **active** |
| `C2` | `$MVA_DATA_ROOT/hf-cache` (`HF_HOME`) | Hub `blobs/` (the real bytes), `snapshots/` symlinks, `trees/`, `refs/`, `xet/` chunk cache, `.incomplete` partials. Currently 468 KB / 4 files — `--local-dir` wrote real files without duplicating blobs here. | `hf cache rm` + `hf cache prune`, then `rm -rf` of the root | **active** |
| `C3` | `<repo>/results/` | Pipeline intermediates derived from C1 — L0 variant calls, per-chromosome BAF and depth vectors, burden; L2 channel D's evidence carries the subject's HPO ids | `rm -rf` contents, keep `.gitkeep` | **active** |
| `C4` | `<repo>/notebooks/` | Saved cell outputs, **if** any cell ever renders patient-derived rows | Strip outputs; verify no genomic content in tracked `.ipynb` | planned |
| `C5` | WSL distro VHDX free space | Remnant blocks from deleted files (deletion is not overwriting) | `fstrim -av` inside WSL after C1–C3 | planned |
| `C6` | Shell history, terminal scrollback, editor workspace state | Only if a record is ever printed | Prevented rather than purged — see below | n/a |
| `C7` | Assistant session transcripts — Claude Code's local `.jsonl` store, plus any claude.ai web sessions held server-side | Prompts and tool output from development. Measured 2026-09-18: variant identifiers derived from the subject are present in the local store (30.6 MB, 9 files). The web sessions were not measured and are enumerated via the account's **Export data**. | Delete the local store; delete the conversations in the account, which removes them from history immediately and from back-end storage within 30 days | **active** |

| `C8` | The organizers' Track 1 submission store — `SageBio/mva-hackathon-2026-leaderboard`, a private HF dataset | The submitted CSV: the subject's candidate variant coordinates in GRCh38, as the Track 1 answer format requires (decision D21). The local copy lives under `results/` and is covered by `C3`. | **Not purgeable by this project.** The organizers issued the data and hold the validated answer; this copy travels back to its source. Recorded for completeness, not because it can be deleted. | live once submitted |

**`C8` is the one row this project cannot act on**, and it is listed for exactly that reason. A
register that shows only the copies we can delete would overstate the reach of the deletion
commitment in [COMPLIANCE.md](../COMPLIANCE.md). The mitigation is upstream of custody: nothing is
disclosed that the recipient did not already hold, the destination is private (the dataset returns
HTTP 401 unauthenticated and the public leaderboard shows scores, not variants), and no coordinate
reaches the repository or the written report.

`C6` is mitigated at the source: L0's guardrail forbids logging patient-derived values at INFO,
and **no pipeline channel sends patient data to an external API** — which is a narrower claim
than it first reads, and deliberately so. It is about the channels. It says nothing about
development assistance, which is `C7` and is a separate surface with its own measurement and
its own terms (decision D17). `C6` is listed because acknowledging an unpurgeable surface
honestly is worth more than omitting it, and `C7` exists because the original wording invited
the reader to conclude something broader than it ever claimed.

### Two traps this register exists to catch

**The Hugging Face cache splits data from its own directory listing.** `snapshots/` contains
*symlinks*; the actual bytes live in `blobs/`, named by hash. Deleting a snapshot folder deletes
pointers and leaves the payload. The default cache root is `~/.cache/huggingface` — outside any
data directory — and `hf_xet` maintains a further chunk cache under it. Interrupted downloads leave
`.incomplete` blobs that `hf cache rm` does not remove; only `hf cache prune` does.
See the [caching guide](https://huggingface.co/docs/huggingface_hub/en/guides/manage-cache).

**Mitigation, applied before the first download:** `HF_HOME` is set to `C2`, inside the custody
root, so every cache surface lands in one registered tree instead of scattering into the home
directory. This is why the environment must be configured before the download starts, not after.

**Where the environment file lives (2026-09-07).** `~/.config/mva/env`, mode `600`, **outside the
repository**. Credentials inside the repo tree are read by editors, file-watchers, and agents that
surface changed files, which put a token into a session transcript three times before this was
fixed. Moving the file out removes that channel rather than relying on `.gitignore` to hold the
line — `.gitignore` was working correctly each time; it was never the failure. `.env.example` in the
repo remains a secret-free template.

**`results/` is gitignored, which makes it feel safe.** It is not patient data in git, but it is
patient-derived data on disk, and it is purged on the same schedule as `C1`.

---

## What the attestation must never contain

The attestation is committed to this repository and emailed to the organizers. It must therefore
carry no patient-derived content by construction:

- **No filenames.** Dataset filenames may embed a sample or subject identifier. **Confirmed on
  2026-09-02:** the VCF filename embeds a lab accession and a sequencer flowcell identifier. This
  rule is not hypothetical for this dataset.
- **No per-file digests.** A file hash is a fingerprint: it lets a holder confirm they have the same
  file, which is a re-identification aid, not a neutral integrity check.
- **No variant coordinates, depths, phenotype text, or counts of anything subject-specific.**

**What identifies the destroyed data instead: the dataset's Hub revision SHA.** That names the exact
snapshot you were granted, is verifiable by the organizers against their own records, and reveals
nothing about the child. Aggregate file counts and total bytes per location are recorded, because
they evidence completeness without describing content.

---

## Attestation format

Emitted by [../src/purge.py](../src/purge.py) as `docs/purge_attestation_<UTC date>.json`, with a
rendered Markdown companion for the email body. Schema version `1.0`:

```json
{
  "schema_version": "1.0",
  "generated_utc": "2026-__-__T__:__:__Z",
  "hackathon_close_date": "2026-10-24",
  "deletion_deadline_utc": "2026-11-23T23:59:00Z",
  "dataset": {
    "repo_id": "SageBio/mva-hackathon-2026-data",
    "repo_type": "dataset",
    "revision_sha": "<full 40-char Hub commit hash>",
    "identification_note": "Identified by Hub revision only. No filenames or per-file digests are recorded, by design."
  },
  "environment": {
    "hostname": "<host>",
    "platform": "<os and kernel>",
    "wsl_distro": "<distro or null>",
    "pipeline_commit": "<git SHA of this repo at purge time>",
    "purge_tool_version": "1.0"
  },
  "locations": [
    {
      "id": "C1",
      "path": "<absolute path>",
      "method": "rm -rf; fstrim",
      "before": { "exists": true, "file_count": 0, "total_bytes": 0 },
      "after":  { "exists": false, "file_count": 0, "total_bytes": 0 },
      "verified_utc": "2026-__-__T__:__:__Z"
    }
  ],
  "checks": [
    { "name": "custody_register_complete", "result": "pass", "detail": "every register row reached a terminal state" },
    { "name": "repo_scan_genomic_extensions", "result": "pass", "matches": 0 },
    { "name": "hf_cache_empty", "result": "pass", "repos_remaining": 0, "bytes_remaining": 0 },
    { "name": "incomplete_blobs_removed", "result": "pass", "matches": 0 }
  ],
  "attested_by": {
    "name": "<attestor>",
    "role": "<role>",
    "statement": "I verified each registered location was emptied on the date recorded above."
  }
}
```

## Verification performed at purge

1. Record `before` state per registered location — existence, file count, total bytes.
2. Delete by the method named in the register.
3. Re-stat each location and record `after`. A location is `purged` only when it is absent, or
   present and empty.
4. Run `hf cache ls` and `hf cache prune` and confirm zero repos and zero `.incomplete` blobs remain.
5. Re-run `tests/test_smoke.py::test_no_patient_data_files_present`, which scans the whole repository
   for genomic extensions.
6. `fstrim -av` so deleted blocks are discarded rather than merely unlinked.
7. Emit the attestation, commit it, and attach it to the confirmation email.

## Honest limits

This attestation proves that the locations tracked in the register were emptied and verified. It
cannot prove that no copy ever existed anywhere else — no destruction record can. Its credibility
rests entirely on the custody surface being deliberately narrow and registered from the start, which
is the argument for keeping all data in one tree under `$MVA_DATA_ROOT` and setting `HF_HOME` inside
it before the first byte is fetched.

## Open items

- **Which date "Hackathon close" means is not settled by the organizers' timeline**, and the two
  readings are a month apart. Taken here as **submission close, 2026-10-24 23:59 UTC**, giving a
  deletion deadline of **2026-11-23 23:59 UTC** — the earlier and therefore safer reading of an
  obligation this project made. The alternative reading, the announced end of the event
  (2026-11-25), would put the deadline at 2026-12-25.

  **The consequence is worth seeing before it arrives:** Track 2 expert-panel judging runs
  2026-10-24 to 2026-11-24 and winners are announced 2026-11-25, so on this reading the data is
  deleted *during* the judging window — one day before it ends. If the panel asks a question that
  needs the pipeline re-run against the real data, that is no longer possible, and the answer is
  that it is no longer possible. Do not extend the deadline to keep the option; ask the organizers
  which date they mean, and if they confirm the later one, record that as a new decision with
  their wording. See [mngmt/decisions.md](../mngmt/decisions.md) D12.
- Confirm with the organizers whether they want the attestation JSON attached to the confirmation
  email, or only the rendered summary.
