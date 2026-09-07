"""L3 -- Claude-in-the-loop reasoning over the aggregated candidates.

Purpose
    For each top-ranked candidate, synthesize a mechanistic rationale from the
    per-channel evidence, **actively search for contradicting evidence**, and emit a
    calibrated confidence. This is the layer that turns five numeric rankings into an
    argument a clinician or judge can actually evaluate.

Inputs
    The aggregated candidate table (RxCUI-keyed, per-channel ranks, provenance) and the
    L1 target sets. Uses the ``anthropic`` package.

Outputs
    Per candidate, written under ``config["results_dir"]/l3/reasoning/``:
      - ``rationale``      -- the mechanistic argument, with citations
      - ``contradicting``  -- evidence found *against* the candidate
      - ``confidence``     -- calibrated, with the reasoning behind the number
      - ``model``, ``prompt_hash`` -- so a judge can reproduce or audit the call

Guardrail
    **No patient data leaves this machine.** Prompts are built from gene symbols, drug
    identifiers, public pathway terms, and channel scores only -- never the subject's
    variants, phenotype narrative, aneuploidy burden vector, or any identifier. This is
    the only layer that talks to an external API, so the boundary is enforced here.

    The contradiction search is not optional decoration. A rationale that only argues
    *for* a candidate is advocacy, not evidence, and would inflate confidence exactly
    where the pipeline is least reliable.

    Every claim must be traceable to a channel output or a resolved citation. Do not let
    the model introduce mechanisms that no channel surfaced -- that is fabrication with a
    persuasive voice.

    Confidence is a hypothesis-strength score. It must never be presented as, or
    convertible into, a probability of clinical benefit.

Notes on the API surface (current as of this scaffold)
    - Model: ``claude-opus-5``. Do not append a date suffix.
    - Thinking is on by default; ``output_config={"effort": ...}`` controls depth. The
      old ``thinking={"type": "enabled", "budget_tokens": N}`` form returns a 400.
    - ``temperature`` / ``top_p`` / ``top_k`` are rejected -- steer with the prompt.
    - Check ``response.stop_reason == "refusal"`` before reading ``response.content``.
"""

from __future__ import annotations

#: Fallback only. The model is config-driven -- read it with :func:`resolve_model` so a
#: judge re-running the pipeline gets the model the report names, and so the choice is
#: recorded in `results/_manifest.json` rather than buried in source. Do not read this
#: constant directly.
DEFAULT_MODEL = "claude-opus-5"


def resolve_model(config: dict) -> str:
    """Return the model id for this run.

    Args:
        config: Parsed pipeline configuration; reads ``model``.

    Returns:
        The configured model id, or :data:`DEFAULT_MODEL` when unset.

    Note:
        Model ids are complete as written -- never append a date suffix. Changing the
        model changes the rationales, so the value is persisted per candidate and in the
        run manifest; a rationale whose model is unknown is not reproducible evidence.
    """
    return str(config.get("model") or DEFAULT_MODEL)


def reason_over_candidates(config: dict) -> None:
    """Synthesize rationale, contradicting evidence, and confidence per candidate.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir`` and ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: implement as follows.
    #  1. Build the evidence packet per candidate from the L3 aggregate: per-channel
    #     ranks, the channel-specific explanations (KG meta-paths, proximity z, proxy
    #     signature id + caveat, HPO terms, graded citations). Assert the packet contains
    #     no patient-derived field before it is serialized -- fail closed.
    #  2. Call `client.messages.create(model=MODEL, ...)` with a stable system prompt
    #     carrying the MVA context and the reasoning contract, and `cache_control` on it
    #     so the shared prefix is cached across candidates; put the per-candidate packet
    #     after the breakpoint.
    #  3. Use structured outputs (`output_config={"format": {...}}`) so rationale /
    #     contradicting_evidence / confidence come back schema-valid rather than parsed
    #     out of prose. Set `effort` to "high" -- this is the rubric's Scientific Rigor.
    #  4. Run a second, adversarial pass whose only job is to argue against the
    #     candidate, and feed its output back before the confidence is finalized.
    #  5. Persist the model id, the prompt hash, and the raw response per candidate.
    #     Handle `stop_reason == "refusal"` explicitly; never silently drop a candidate.
    raise NotImplementedError("claude_reasoning.reason_over_candidates is a scaffold stub")
