# Data

No patient data lives in this repo, ever. The dataset is gated:
https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data

Data lives outside the repository entirely, under `$MVA_DATA_ROOT` on WSL-native ext4 (see
`.env.example`), so git cannot reach it even if `.gitignore` were broken. Every location holding
patient-derived bytes is registered in [docs/data_custody.md](docs/data_custody.md).

**Pinned revision:** `59e322d27f399006b398d366d33e703e48a29914` — the identifier used in the
destruction attestation, and the one to re-download against for reproducibility.

## Inventory (verified 2026-09-02)

13 objects in the repo: **11 data files** plus `README.md` and `.gitattributes`.

| Role | Type | Count | Size | Downloaded |
|---|---|---|---|---|
| Raw reads | `.fastq.gz` | 8 | 84.67 GB | **No** — deliberately not fetched |
| Small variants, single sample | `.vcf.gz` | 1 | 0.32 GB | Yes |
| Tabix index for the VCF | `.vcf.gz.tbi` | 1 | 2.3 MB | Yes |
| Clinical phenotype narrative | `.docx` | 1 | 17 KB | Yes |
| Dataset documentation | `.md` | 1 | 5 KB | Yes |

**Filenames are not recorded here on purpose.** They embed a lab accession and a sequencer
flowcell identifier. Roles, types, and sizes evidence the inventory without carrying an identifier
into a committed file — the same rule the attestation follows.

## Findings that change the architecture

**There is no BAM.** The dataset ships raw reads and a called VCF, not alignments. Module
docstrings that say "WGS VCF (+ BAM)" describe an input that does not exist. Producing a BAM means
aligning 84.67 GB of FASTQ — reference download, index build, and a multi-hour-to-multi-day
alignment — which is out of scope for the schedule.

Consequence for L0's aneuploidy-burden feature: it is derived from the **VCF**, not from BAM depth.
The VCF carries `FORMAT/AD` (per-allele depths) and `FORMAT/DP`, so per-chromosome **B-allele
frequency at heterozygous sites** gives both detection and a mosaic-fraction estimate, and
median `DP` per chromosome serves as a secondary depth signal. Allele ratios at a single locus are
internally normalized, so the GC-bias and mappability corrections a depth-based method needs do not
arise. See [src/l0_genomics/run.py](src/l0_genomics/run.py) for the implementation and its citation.

The FASTQ files remain available at the pinned revision if alignment is ever revisited. Not
downloading them keeps 85 GB out of the custody register and out of the deletion obligation.

**There is no expression data** — resolved, and this closes the original open question. Channel C
therefore runs on a **PROXY** signature (LINCS L1000 knockdown of the causal gene, or a curated
aneuploidy-response gene set). The substitution and its limits travel with every row, figure, and
table derived from that channel; see
[src/l2_channels/channel_c_signature.py](src/l2_channels/channel_c_signature.py).

## Deletion

All data is deleted within 30 days of Hackathon close, evidenced by an attestation from
[src/purge.py](src/purge.py) and confirmed by email to the organizers. See
[COMPLIANCE.md](COMPLIANCE.md) and [docs/data_custody.md](docs/data_custody.md).
