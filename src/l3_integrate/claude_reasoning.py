"""L3 -- Local-model reasoning over the aggregated candidates.

Purpose
    For each top-ranked candidate, synthesize a mechanistic rationale from the
    per-channel evidence, **actively search for contradicting evidence**, and emit a
    calibrated confidence. This is the layer that turns five numeric rankings into an
    argument a clinician or judge can actually evaluate.

Inputs
    The aggregated candidate table (parent-molecule-keyed, per-channel ranks, provenance)
    and the L1 target sets. Runs a **local** open-weights model; no hosted API, no
    ``anthropic`` package, no network call of any kind.

Outputs
    Per candidate, written under ``config["results_dir"]/l3/reasoning/``:
      - ``rationale``      -- the mechanistic argument, with citations
      - ``contradicting``  -- evidence found *against* the candidate
      - ``confidence``     -- computed from a rubric the model fills, not asked for
                              directly; see the guardrail
      - ``model``, ``prompt_hash`` -- so a judge can reproduce or audit the call

Guardrail
    **Nothing leaves this machine.** Decision D17 (``mngmt/decisions.md``) puts the
    reasoning step on a local model precisely so this stops being a claim about somebody's
    terms of service and becomes a property of where the process runs. The packet is still
    asserted free of patient-derived fields before it is serialized -- fail closed -- because
    the guardrail must hold if a later session ever reconsiders the runtime.

    The contradiction search is not optional decoration. A rationale that only argues
    *for* a candidate is advocacy, not evidence, and would inflate confidence exactly
    where the pipeline is least reliable.

    Every claim must be traceable to a channel output or a resolved citation. Do not let
    the model introduce mechanisms that no channel surfaced -- that is fabrication with a
    persuasive voice. At this model size that is **enforced, not instructed**: every
    citation key in the output is checked against the supplied packet and dropped if it is
    not there.

    Confidence is a hypothesis-strength score. It must never be presented as, or
    convertible into, a probability of clinical benefit. A 7B model is poorly calibrated,
    so it is never asked for the number -- it fills a rubric (evidence grade, channel
    count, contradiction found or not) and the score is computed here, where it can be
    tested.

Notes on the local runtime (D17)
    - 6 GB of VRAM caps this at a **7-8B instruct model at Q4_K_M**, roughly 4.5-5 GB of
      weights plus KV cache. Expect one to two hours for a full pass over the survivors,
      which D1 already accounts for: this is a batch pipeline, not an application.
    - Determinism comes from temperature 0 and a fixed seed, both from ``config``. The
      resolved model id is persisted per candidate and in the run manifest -- a rationale
      whose model is unknown is not reproducible evidence.
    - Structured output uses **constrained decoding** (a grammar), since there is no
      server-side schema enforcement locally. Parse-and-retry is the fallback, never
      free-text scraping.
"""

from __future__ import annotations

#: Fallback only. The model is config-driven -- read it with :func:`resolve_model` so a
#: judge re-running the pipeline gets the model the report names, and so the choice is
#: recorded in `results/_manifest.json` rather than buried in source. Do not read this
#: constant directly.
#:
#: Provisional: this names the *class* D17 settled on -- a 7-8B Q4_K_M instruct model that
#: fits 6 GB -- not a build that has been benchmarked here. The exact file and runtime are
#: fixed when the step is implemented, and whatever runs is what the manifest records.
DEFAULT_MODEL = "qwen2.5-7b-instruct-q4_k_m"


def resolve_model(config: dict) -> str:
    """Return the model id for this run.

    Args:
        config: Parsed pipeline configuration; reads ``model``.

    Returns:
        The configured model id, or :data:`DEFAULT_MODEL` when unset.

    Note:
        Changing the model changes the rationales, so the value is persisted per candidate
        and in the run manifest.
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
    #  2. Load the local model once and reuse it across candidates; a per-candidate reload
    #     dominates the runtime at this size. Keep the MVA context and the reasoning
    #     contract in a stable system prefix so the KV cache survives between candidates.
    #  3. Constrain generation with a grammar so rationale / contradicting_evidence /
    #     the confidence rubric come back schema-valid rather than parsed out of prose.
    #  4. Run a second, adversarial pass whose only job is to argue against the
    #     candidate, and feed its output back before the rubric is finalized.
    #  5. Validate every citation key against the packet; drop unknown keys and record
    #     that they were dropped. Compute the confidence score here from the rubric.
    #  6. Persist the model id, the prompt hash, and the raw response per candidate.
    #     Never silently drop a candidate -- a refusal or a parse failure is recorded.
    raise NotImplementedError("claude_reasoning.reason_over_candidates is a scaffold stub")
