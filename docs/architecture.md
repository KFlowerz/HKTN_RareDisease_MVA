# Architecture

> **Placeholder — paste the architecture analysis here.**
>
> The headings below mirror the L0–L5 package structure so the full write-up has a home. Each section
> should end up covering: purpose, inputs, outputs, method, failure modes, and the guardrail that keeps
> the layer inside scope. Keep this document and the module docstrings in sync.

## Design premise: multi-channel evidence

Parallel weak signals combined by consensus, chosen because any single repurposing method is fragile on a
hyper-rare disease. Convergence across independent channels is the signal; no single channel score is
trusted on its own.

*paste the architecture analysis here*

## L0 — `l0_genomics`

Ingest the single-sample WGS VCF; annotate variants (VEP / OpenCRAVAT-style); focus on biallelic
SAC-gene variants; reconcile with the Track-1 validated causal variant. Innovation hook:
per-chromosome **B-allele frequency** as an **aneuploidy burden** feature.

**Input correction (2026-09-02).** The design originally assumed a BAM and derived burden from
per-chromosome read depth. The gated dataset ships raw reads and a called VCF, with no alignments
(dataset revision `59e322d2…`; see [../DATA.md](../DATA.md)). Burden is therefore computed from
`FORMAT/AD` at heterozygous sites: heterozygous BAF clusters near 0.5 in disomic regions, splits
under a mosaic gain or loss, and the magnitude of the split estimates the mosaic fraction
(Conlin et al., 2010; Loh et al., 2018).

This is a stronger feature than the one it replaces, not a concession. An allele ratio is measured
between two alleles at one locus, so depth, GC content, and mappability affect numerator and
denominator identically and cancel — the corrections a depth-based caller requires never arise. And
the statistic estimates *mosaic fraction*, which is the quantity MVA is named for, rather than mere
over- or under-representation. Median normalized `FORMAT/DP` per chromosome is retained as a weaker
secondary signal, never reported alone.

*paste the architecture analysis here*

## L1 — `l1_target`

Causal gene → protein/complex → interactome disease module (Reactome / STRING / Open Targets). Upstream
target set (restore mitotic fidelity) vs. downstream target set (buffer aneuploidy stress /
chemoprevention).

*paste the architecture analysis here*

## L2 — `l2_channels`

Five parallel candidate generators: A knowledge graph, B network proximity, C signature reversion on a
proxy signature, D phenotype/HPO, E literature prior.

*paste the architecture analysis here*

## L3 — `l3_integrate`

RxCUI identity harmonization, rank aggregation (RRA/Borda) rewarding cross-channel convergence, and the
Claude-in-the-loop reasoning step (rationale synthesis, active search for contradicting evidence,
calibrated confidence).

*paste the architecture analysis here*

## L4 — `l4_validate`

In-silico validation, pediatric safety triage with hard genotoxic exclusion, and the blinded internal
benchmark. Drug-annotation sourcing and licensing live in
[../src/l4_validate/sources.md](../src/l4_validate/sources.md).

*paste the architecture analysis here*

## L5 — `l5_report`

Rubric-aligned figures, candidate tables, video assets.

*paste the architecture analysis here*

## Judging rubric this architecture serves

| Criterion | Weight | Where it is addressed |
|---|---|---|
| Scientific Rigor | 35% | Multi-channel consensus, proxy-signature honesty, blinded benchmark (L2/L3/L4) |
| Impact | 25% | Safety-triaged shortlist of approved drugs for a disease with no therapy (L4/L5) |
| Innovation | 25% | Aneuploidy burden and mosaic fraction from B-allele frequency (L0), gene-level KG anchoring (L2A), Claude-in-the-loop contradiction search (L3) |
| Scalability | 15% | Config-driven, disease-agnostic layer boundaries; nothing about MVA hardcoded |
