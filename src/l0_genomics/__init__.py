"""L0 -- Genomics: WGS VCF to causal gene, variant effect, and aneuploidy burden.

Burden comes from B-allele frequency, not read depth: the dataset ships no BAM.
See ``run.py`` and ``DATA.md``.
"""

from .run import run

__all__ = ["run"]
