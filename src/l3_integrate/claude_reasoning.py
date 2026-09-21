"""L3 -- Local-model reasoning over the aggregated candidates.

Purpose
    For each surviving candidate, synthesize a mechanistic rationale from the per-channel
    evidence, **actively search for contradicting evidence**, and emit a calibrated
    confidence. This is the layer that turns numeric rankings into an argument a clinician
    or judge can actually evaluate.

Inputs
    The aggregated candidate table (parent-molecule-keyed, per-channel ranks, provenance)
    and the L1 target sets. Runs a **local** open-weights model behind an
    OpenAI-compatible endpoint -- by default llama.cpp's server bound to loopback.

Outputs
    Per candidate, written under ``config["results_dir"]/l3/reasoning/``:
      - ``rationale``      -- the mechanistic argument, citing only supplied keys
      - ``contradicting``  -- evidence found *against* the candidate
      - ``confidence``     -- computed here from the rubric the model filled; never
                              asked for directly, see the guardrail
      - ``model``, ``prompt_hash``, ``dropped_keys`` -- so a judge can reproduce or audit

Guardrail
    **Nothing leaves this machine.** Decision D17 puts the reasoning step on a local
    model. The endpoint must be loopback-bound: :func:`resolve_endpoint` refuses any host
    that is not ``127.0.0.1``/``localhost``, because a bind to any other interface makes
    this a hosted API and reopens every question D17 closed. The packet is still asserted
    free of patient-derived fields before it is serialized -- fail closed -- because the
    guardrail must hold if a later session reconsiders the runtime.

    The contradiction search is not optional decoration. A rationale that only argues
    *for* a candidate is advocacy, not evidence, and would inflate confidence exactly
    where the pipeline is least reliable.

    **Every claim must be traceable to a channel output.** At this model size that is
    enforced, not instructed: every citation key in the output is checked against the
    supplied packet and dropped if absent, and the drop is recorded rather than hidden.
    Keys come back bracketed (``[zerbib2024]``), so they are normalized before matching.

    **Confidence is computed, not asked for.** A 7B is poorly calibrated, and the spike
    of 2026-09-20 showed the miscalibration is not confined to the number: asked to grade
    evidence that was explicitly cell-line work, the model answered ``clinical``. So the
    evidence grade is **supplied** from channel E's own grading and echoed back, never
    generated, and the score is computed in Python where it can be tested. The model is
    asked only what it can do -- does the supplied evidence support the mechanism, is
    there a contradiction in it, which supplied keys are load-bearing.

    Confidence is a hypothesis-strength score. It must never be presented as, or
    convertible into, a probability of clinical benefit.

Notes on the local runtime (D17, measured 2026-09-20)
    - **Vulkan, not CUDA.** Current llama.cpp CUDA builds no longer ship Pascal kernels
      and abort with ``invalid device function``; CUDA 13 dropped Pascal outright. The
      Vulkan backend runs the same model and is the more portable choice.
    - On a GTX 1060 6 GB: 175-181 tok/s prompt, ~24 tok/s generation, 5,370 MiB peak at
      8192 context. A full pass over the survivors is roughly 1.6 hours -- a batch
      pipeline (D1), not an application.
    - Determinism comes from ``temperature 0`` plus a fixed seed, both from ``config``;
      three identical requests returned byte-identical output. The resolved model id is
      persisted per candidate and in the run manifest -- a rationale whose model is
      unknown is not reproducible evidence.
    - Structure comes from ``response_format: json_schema``. The in-process routes do not
      work: ``llama-completion`` exits silently in b11065, and a grammar in ``llama-cli``
      collides with the chat template's control tokens.
"""

from __future__ import annotations

from urllib.parse import urlparse

#: Fallback only. The model is config-driven -- read it with :func:`resolve_model` so a
#: judge re-running the pipeline gets the model the report names, and so the choice is
#: recorded in `results/_manifest.json` rather than buried in source. Do not read this
#: constant directly.
#:
#: This is the model the 2026-09-20 spike measured, not merely a plausible default.
DEFAULT_MODEL = "qwen2.5-7b-instruct-q4_k_m"

#: Default endpoint. Loopback by construction; see :func:`resolve_endpoint`.
DEFAULT_ENDPOINT = "http://127.0.0.1:8080"

#: The only hosts that keep D17's guarantee. A bind to anything else is a hosted API.
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1", "[::1]"})


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


def resolve_endpoint(config: dict) -> str:
    """Return the reasoning endpoint, refusing anything that is not loopback.

    Args:
        config: Parsed pipeline configuration; reads ``reasoning_endpoint``.

    Returns:
        The endpoint URL.

    Raises:
        ValueError: If the host is not a loopback address. This is D17's guarantee
            expressed as code: a local model behind a loopback socket keeps the data on
            the machine, while the same model behind a routable address does not, and
            nothing else about the call would look different.
    """
    endpoint = str(config.get("reasoning_endpoint") or DEFAULT_ENDPOINT)
    host = urlparse(endpoint).hostname
    if host not in LOOPBACK_HOSTS:
        raise ValueError(
            f"reasoning_endpoint {endpoint!r} binds host {host!r}, which is not loopback. "
            "Decision D17 requires the reasoning step to run locally; a non-loopback "
            "endpoint is a hosted API and must supersede D17 before it is used."
        )
    return endpoint


def reason_over_candidates(config: dict) -> None:
    """Synthesize rationale, contradicting evidence, and confidence per candidate.

    Args:
        config: Parsed pipeline configuration; uses ``results_dir`` and ``seed``.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: implement as follows.
    #  1. Build the evidence packet per candidate from the L3 aggregate: per-channel
    #     ranks, the channel-specific explanations (proximity z, HPO-derived terms,
    #     graded citations from channel E). Carry channel E's evidence grade INTO the
    #     packet -- it is supplied, never asked for (D17, correction 3). Assert the packet
    #     contains no patient-derived field before it is serialized -- fail closed.
    #  2. POST to `resolve_endpoint(config)` + "/v1/chat/completions" with
    #     `response_format: {"type": "json_schema", ...}`, temperature 0 and the config
    #     seed. Keep the MVA context and the reasoning contract in a stable system
    #     message so the server's prefix cache survives between candidates.
    #  3. Run a second, adversarial pass whose only job is to argue against the
    #     candidate, and feed its output back before the rubric is finalized.
    #  4. Normalize and validate every citation key against the packet; drop unknown keys,
    #     record them in `dropped_keys`, and never let a dropped key vanish silently.
    #  5. Compute the confidence score here from the rubric fields plus the supplied
    #     grade and `n_channels_supporting`. The model never sees a number.
    #  6. Persist the model id, the endpoint, the prompt hash and the raw response per
    #     candidate. A refusal or a parse failure is recorded, never a dropped candidate.
    raise NotImplementedError("claude_reasoning.reason_over_candidates is a scaffold stub")
