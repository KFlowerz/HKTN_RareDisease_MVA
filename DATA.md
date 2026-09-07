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

**The VCF is GRCh38 with no `chr` prefix, and declares 2,580 contigs** (verified 2026-09-07):
`GCA_000001405.15_..._plus_hs38d1`, GATK `VariantFiltration`-hardened, single sample, ~5.0M records
carrying `FORMAT/AD` and `FORMAT/DP`. Two consequences for L0: annotation resources using `chr`
must be translated in one documented place, and per-chromosome statistics must be restricted to the
primary contigs so decoys and alts do not enter the burden vector. Heterozygous-site density at
`DP >= 10` is ample on every autosome, so the per-sample autosomal distribution supplies the diploid
baseline directly — no external calibration constant is needed.

**The phenotype document is already HPO-coded** (verified 2026-09-07), so Channel D needs no text
mining. *The terms themselves are not recorded here.* They are clinical facts about a living child
and a specific combination is identifying in a population of roughly 50 patients; they are parsed
from `data_dir` at runtime and never committed. See [COMPLIANCE.md](COMPLIANCE.md).

**No causal gene is named anywhere in the dataset.** The reconciliation step L0 was specified around
has no counterpart — there is no Track-1 answer to check against — so L0 must make the call
independently and ship the evidence that supports it, including what it rejected.

**There is no expression data** — resolved, and this closes the original open question. Channel C
therefore runs on a **PROXY** signature (LINCS L1000 knockdown of the causal gene, or a curated
aneuploidy-response gene set). The substitution and its limits travel with every row, figure, and
table derived from that channel; see
[src/l2_channels/channel_c_signature.py](src/l2_channels/channel_c_signature.py).

## Deletion

All data is deleted within 30 days of Hackathon close, evidenced by an attestation from
[src/purge.py](src/purge.py) and confirmed by email to the organizers. See
[COMPLIANCE.md](COMPLIANCE.md) and [docs/data_custody.md](docs/data_custody.md).
