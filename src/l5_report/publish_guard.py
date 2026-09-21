"""L5 -- the publication boundary.

Purpose
    Check everything L5 is about to write for anything that must not be published, and
    refuse rather than write it.

Inputs
    Rendered strings -- HTML, table rows, figure captions, JSON blobs.

Outputs
    Nothing on success. A :class:`ValueError` naming the offending text on failure.

Guardrail
    **This is the last gate before a file a judge will read.** Everywhere upstream the
    same patterns are checked against things that leave the *process*; here they are
    checked against things that leave the *project*. The two are not the same boundary and
    this one is stricter, because the report is the artifact that gets shared.

    It fails closed and it checks the **rendered** form, not the fields, so a value that
    reached the page through an unexpected path is still caught.

    The identifier patterns are imported from :mod:`src.l3_integrate.claude_reasoning`
    rather than restated, so the two boundaries cannot drift apart. What is added here is
    report-specific: the per-chromosome aneuploidy detail that is identifying in a
    ~50-patient population, and the endpoint wording that D4 requires.
"""

from __future__ import annotations

import re

from ..l3_integrate.claude_reasoning import FORBIDDEN_IN_PACKET

#: Per-chromosome aneuploidy results are **not published**, in any form.
#:
#: The pipeline's B-allele-frequency method is an innovation claim and L5 describes it.
#: What L5 never shows is this child's result from it. MVA case series publish karyotypes,
#: the worldwide population is around fifty, and a per-chromosome pattern -- which
#: contigs, how many, at what mosaic fraction -- is close to a fingerprint. The method is
#: publishable; the per-contig finding is not, and the report says that it is withheld
#: rather than omitting it silently.
WITHHELD_BURDEN_FIELDS = (
    "per_contig", "mosaic_fraction_if_gain", "mosaic_fraction_if_loss",
    "deconvolved_shift", "n_flagged_windows", "longest_contiguous_flagged",
    "excess_over_baseline", "median_abs_dev", "inferred_sex", "contamination_screen",
)

#: Report-surface patterns, on top of the shared identifier set.
#:
#: Deliberately about *values*, not field names. An earlier version matched
#: ``per_contig`` and ``mosaic_fraction_if_gain`` as text, which rejected
#: :func:`strip_withheld`'s own output -- the report names the fields it is withholding,
#: because withholding silently is indistinguishable from never having looked. Field
#: names are checked structurally instead, by :func:`assert_no_burden_detail`.
FORBIDDEN_IN_REPORT = FORBIDDEN_IN_PACKET + (
    (re.compile(r"\bkaryotype\s*[:=]\s*\d"), "an explicit karyotype"),
    (re.compile(r"\b4[5-9],\s*(?:XX|XY)\b"), "an aneuploid karyotype string"),
    (re.compile(r"\brs\d{4,}\b"), "a dbSNP identifier"),
)

#: D4: the endpoint is *secondary* prevention. The bare word overstates it, so the report
#: may use it only where the qualifier is already present in the same sentence.
ENDPOINT_WORD = re.compile(r"\bchemoprevention\b", re.I)
ENDPOINT_QUALIFIER = re.compile(r"\bsecondary\b", re.I)


def check_text(text: str, where: str) -> None:
    """Raise if ``text`` carries anything that must not be published.

    Args:
        text: The rendered string about to be written.
        where: A path or label naming what is being written, for the message.

    Raises:
        ValueError: On the first forbidden pattern found.
    """
    for pattern, description in FORBIDDEN_IN_REPORT:
        found = pattern.search(text)
        if found:
            raise ValueError(
                f"refusing to publish {where}: it contains {description} "
                f"({found.group(0)!r}). Nothing published may re-identify the child or "
                "family (CLAUDE.md constraint 7 and the L5 guardrail)."
            )
    _check_endpoint_wording(text, where)


def _check_endpoint_wording(text: str, where: str) -> None:
    """Every sentence naming the endpoint must qualify it as *secondary* prevention (D4).

    Checked per sentence rather than per document: a qualifier in the introduction does
    not travel with a table caption that a reader screenshots, and this report is about a
    child.
    """
    for sentence in re.split(r"(?<=[.!?])\s+|\n", text):
        if ENDPOINT_WORD.search(sentence) and not ENDPOINT_QUALIFIER.search(sentence):
            raise ValueError(
                f"refusing to publish {where}: {sentence.strip()[:160]!r} names the "
                "endpoint without qualifying it as SECONDARY prevention. Decision D4 "
                "fixed the endpoint as recurrence and second-primary risk reduction, not "
                "prevention of a first cancer, and the unqualified word overstates it."
            )


def strip_withheld(burden: dict) -> dict:
    """The publishable part of L0's aneuploidy burden artifact: the method, not the result.

    Returns the method description, its citations and the sensitivity envelope -- what the
    pipeline *can* detect -- with every per-chromosome field removed. The innovation claim
    is that allele ratios alone resolve mosaic aneuploidy without alignments; making that
    claim needs the sensitivity, not the subject's chromosomes.

    Args:
        burden: the parsed ``results/l0_genomics/aneuploidy_burden.json``.
    """
    sensitivity = dict(burden.get("sensitivity") or {})
    envelope = {
        key: sensitivity[key]
        for key in ("min_detectable_mosaic_fraction",
                    "min_detectable_mosaic_fraction_conservative", "z_threshold")
        if key in sensitivity
    }
    return {
        "method": burden.get("method", {}),
        "sensitivity_envelope": envelope,
        "withheld": sorted(WITHHELD_BURDEN_FIELDS),
        "withheld_reason": (
            "Per-chromosome aneuploidy results are withheld. MVA case series publish "
            "karyotypes and the worldwide population is around fifty patients, so which "
            "chromosomes are involved, how many, and at what mosaic fraction is close to "
            "an identifier. The method and what it can resolve are described; this "
            "child's result from it is not published."
        ),
    }


def assert_no_burden_detail(payload, where: str) -> None:
    """Raise if any **mapping key** anywhere in ``payload`` is a withheld burden field.

    Structural rather than textual, so the report can name the fields it withholds
    without tripping the check that withholds them. What is forbidden is a per-chromosome
    *result* travelling into a published file, which always arrives as a key with a value
    under it.
    """
    stack = [payload]
    while stack:
        current = stack.pop()
        if isinstance(current, dict):
            for key, value in current.items():
                if isinstance(key, str) and key in WITHHELD_BURDEN_FIELDS:
                    raise ValueError(
                        f"refusing to publish {where}: it carries the per-chromosome "
                        f"aneuploidy field {key!r}. In a ~50-patient population that "
                        "pattern is close to an identifier; publish the method and its "
                        "sensitivity, never the subject's result (L5 guardrail)."
                    )
                stack.append(value)
        elif isinstance(current, (list, tuple)):
            stack.extend(current)


def assert_publishable(payload, where: str) -> None:
    """Raise if anything in ``payload`` must not be published.

    Accepts a string or any JSON-serialisable structure. Text patterns are checked against
    the serialised form, so a value nested somewhere unexpected is still caught; burden
    fields are checked structurally.
    """
    import json

    if not isinstance(payload, str):
        assert_no_burden_detail(payload, where)
    text = payload if isinstance(payload, str) else json.dumps(payload, ensure_ascii=False,
                                                               default=str)
    check_text(text, where)
