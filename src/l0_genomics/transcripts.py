"""L0 -- transcript policy and the variant-call record (decision D5).

Purpose
    Decide which transcripts may make a variant "loss-of-function", and fix the shape of
    the record every L0 variant call carries, so the policy that produced a call travels
    with it. Decided 2026-09-11 (``mngmt/decisions.md`` D5) from the review and data in
    ``docs/research/transcript-policy-lof.md``.

Inputs
    ``config["annotator"]``: ``transcript_policy`` and ``mane_release``.

Outputs
    :func:`snpeff_passes` -- the annotation passes to run, each with its snpEff flags;
    :func:`lof_tier` -- the tier a variant earns from those passes;
    :class:`VariantCall` -- the per-variant record L0 emits.

Guardrail
    **MANE Select makes the call; nothing else can.** A variant that is LoF only on a
    non-MANE transcript is kept as tier ``non_mane``, naming the transcript -- never
    dropped, never promoted. LoF confined to a subset of a gene's transcripts is common and
    enriched for false positives [macarthur2012] doi:10.1126/science.1215040
    [singerberk2023] doi:10.1016/j.ajhg.2023.08.005, while all but 67 of 33,736
    PubMed-supported pathogenic variants map to MANE Select [pozo2022]
    doi:10.1038/s41525-022-00329-6 -- rare enough to demote, not rare enough to discard in
    an n=1 search.

    snpEff ``-canon`` is refused by name: it selects the longest CDS, not MANE
    [snpeffdocs], and on the SAC panel it picks a different transcript for BUB1, BUB1B and
    CEP192 (``results/transcript_policy/picks.tsv``, seed=42).

    A missing or unknown policy raises. There is no default: a silently defaulted policy is
    indistinguishable from a chosen one.

    ``hgvs_c`` / ``hgvs_p`` identify the child's variant. A :class:`VariantCall` is written
    under ``results_dir`` only, and L5 reports it in aggregate or categorical form, never
    verbatim.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Policies L0 can run, each as ``((pass_name, snpeff_flags), ...)`` with the call-making
#: pass first. Adding one needs a decision in mngmt/decisions.md.
TRANSCRIPT_POLICIES = {
    # D5: MANE Select makes the call; an all-transcript pass keeps non-MANE LoF, flagged.
    "mane_select_tiered": (("primary", ("-tag", "MANE_Select")), ("secondary", ())),
}

#: Refused by name, so the error says why rather than just "unknown".
REFUSED_POLICIES = {
    "canon": "snpEff -canon selects the longest CDS, not MANE Select (decision D5)",
    "all_transcripts": "LoF on any transcript as the call admits partial-LoF false "
                       "positives (decision D5)",
}

PASSES = ("primary", "secondary")
LOF_TIERS = ("mane_select", "non_mane", "not_lof")
IMPACTS = ("HIGH", "MODERATE", "LOW", "MODIFIER")

#: The only pass that may produce each LoF tier. ``not_lof`` may come from either.
_TIER_PASS = {"mane_select": "primary", "non_mane": "secondary"}


def transcript_policy(config: dict) -> str:
    """Return the configured policy name, refusing a missing, refused, or unknown one."""
    try:
        policy = config["annotator"]["transcript_policy"]
    except (KeyError, TypeError) as exc:
        raise ValueError(
            "config['annotator']['transcript_policy'] is required (decision D5); "
            "there is no default"
        ) from exc
    if policy in REFUSED_POLICIES:
        raise ValueError(f"transcript_policy {policy!r} is refused: {REFUSED_POLICIES[policy]}")
    if policy not in TRANSCRIPT_POLICIES:
        raise ValueError(
            f"unknown transcript_policy {policy!r}; known: {sorted(TRANSCRIPT_POLICIES)}"
        )
    return policy


def snpeff_passes(config: dict) -> tuple:
    """Return ``((pass_name, snpeff_flags), ...)`` for the configured policy, primary first."""
    return TRANSCRIPT_POLICIES[transcript_policy(config)]


def lof_tier(*, primary_lof: bool, secondary_lof: bool) -> str:
    """Tier a variant from whether each pass called it LoF.

    Raises:
        ValueError: If the primary pass calls LoF but the secondary does not. The
            all-transcript pass includes the MANE Select transcript, so that combination
            means the passes were run or read inconsistently.
    """
    if primary_lof and not secondary_lof:
        raise ValueError(
            "primary (MANE Select) LoF without secondary (all-transcript) LoF: the "
            "secondary pass includes MANE Select, so the passes are inconsistent"
        )
    if primary_lof:
        return "mane_select"
    if secondary_lof:
        return "non_mane"
    return "not_lof"


@dataclass(frozen=True)
class VariantCall:
    """One annotated variant as L0 emits it.

    ``transcript``, ``mane_release`` and ``lof_tier`` are required by decision D5, and
    ``annotation_pass`` / ``transcript_policy`` record how the call was produced, so it can
    be regenerated and disagreed with.
    """

    gene: str
    transcript: str          # the transcript the consequence was read from
    annotation_pass: str     # "primary" (MANE Select) | "secondary" (all transcripts)
    lof_tier: str            # "mane_select" | "non_mane" | "not_lof"
    mane_release: str        # MANE release the database's tags were checked against
    transcript_policy: str   # the policy that ran
    consequence: str         # snpEff effect, e.g. "stop_gained"
    impact: str              # snpEff impact
    hgvs_c: str = ""         # identifying -- results_dir only
    hgvs_p: str = ""         # identifying -- results_dir only

    def __post_init__(self) -> None:
        if not self.transcript:
            raise ValueError("VariantCall.transcript is required (decision D5)")
        if not self.mane_release:
            raise ValueError("VariantCall.mane_release is required (decision D5)")
        if self.annotation_pass not in PASSES:
            raise ValueError(f"annotation_pass {self.annotation_pass!r} not in {PASSES}")
        if self.lof_tier not in LOF_TIERS:
            raise ValueError(f"lof_tier {self.lof_tier!r} not in {LOF_TIERS}")
        if self.impact not in IMPACTS:
            raise ValueError(f"impact {self.impact!r} not in {IMPACTS}")
        if self.transcript_policy not in TRANSCRIPT_POLICIES:
            raise ValueError(f"transcript_policy {self.transcript_policy!r} is not a known policy")
        required = _TIER_PASS.get(self.lof_tier)
        if required and self.annotation_pass != required:
            raise ValueError(
                f"lof_tier {self.lof_tier!r} can only come from the {required} pass, "
                f"not {self.annotation_pass!r}"
            )
