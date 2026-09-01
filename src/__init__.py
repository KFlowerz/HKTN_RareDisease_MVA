"""MVA Hackathon 2026 — Track 2 (Drug Repurposing).

A multi-channel, in-silico drug repurposing pipeline for mosaic variegated aneuploidy
(MVA, OMIM 257300). Layers L0-L5 run in order; see ``src.pipeline``.

Everything this package produces is a **computational hypothesis**. It is not a clinical
recommendation, not an efficacy claim, and not a substitute for clinical judgement.
"""

__all__ = [
    "l0_genomics",
    "l1_target",
    "l2_channels",
    "l3_integrate",
    "l4_validate",
    "l5_report",
    "pipeline",
]
