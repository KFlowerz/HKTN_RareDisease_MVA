"""L3 -- Local-model reasoning over the aggregated candidates.

Purpose
    For each candidate, synthesize a mechanistic rationale from the per-channel evidence,
    **actively search for contradicting evidence**, and emit a calibrated confidence.
    This is the layer that turns numeric rankings into an argument a clinician or judge
    can actually evaluate.

Inputs
    ``results/l3/candidates.tsv`` for the aggregated ranking, the per-channel tables under
    ``results/l2/`` for the channel-specific explanations, and channel E's
    ``references.json`` for verified citations. Runs a **local** open-weights model behind
    an OpenAI-compatible endpoint -- by default llama.cpp's server bound to loopback.

Outputs
    Written under ``config["results_dir"]/l3/reasoning/``:
      - ``rationales.json`` -- per candidate: rationale, contradicting evidence, the
        rubric the model filled, the computed confidence, kept and dropped citation keys
      - ``reasoning.json``  -- model, endpoint, seed, scope, counts, failures, caveats

Guardrail
    **Nothing leaves this machine.** Decision D17 puts the reasoning step on a local
    model. The endpoint must be loopback-bound: :func:`resolve_endpoint` refuses any host
    that is not loopback, because a bind to any other interface makes this a hosted API
    and reopens every question D17 closed.

    **No patient-derived field may enter a packet.** :func:`assert_publishable` checks
    every packet before it is serialized and raises rather than sending -- fail closed.
    Channel D's ``evidence.json`` names the diseases that resemble the subject and is
    never read here; only its rank reaches the packet, exactly as L3 itself does it.

    The contradiction search is not optional decoration. A rationale that only argues
    *for* a candidate is advocacy, not evidence, and would inflate confidence exactly
    where the pipeline is least reliable.

    **Every claim must be traceable to a channel output.** At this model size that is
    enforced, not instructed: every citation key in the output is checked against the
    supplied packet and dropped if absent, and the drop is recorded rather than hidden.
    Keys come back bracketed (``[DOI:...]``), so they are normalized before matching.

    **Confidence is computed, not asked for.** A 7B is poorly calibrated, and the spike of
    2026-09-20 showed the miscalibration is not confined to the number: asked to grade
    evidence that was explicitly cell-line work, the model answered ``clinical``. So the
    evidence grade is **supplied** from channel E's own counts and never generated, and
    the score is computed in :func:`confidence` where it can be tested. The model is asked
    only what it can do -- does the supplied evidence support the mechanism, is there a
    contradiction in it, which supplied keys are load-bearing.

    Confidence is a hypothesis-strength score. It must never be presented as, or
    convertible into, a probability of clinical benefit.

Notes on the local runtime (D17, measured 2026-09-20)
    - **Vulkan, not CUDA.** Current llama.cpp CUDA builds no longer ship Pascal kernels
      and abort with ``invalid device function``; CUDA 13 dropped Pascal outright.
    - On a GTX 1060 6 GB: 175-181 tok/s prompt, ~24 tok/s generation, 5,370 MiB peak at
      8192 context. A full pass over the survivors is roughly 1.6 hours -- a batch
      pipeline (D1), not an application.
    - Determinism comes from ``temperature 0`` plus a fixed seed, both from ``config``.
    - Structure comes from ``response_format: json_schema``. The in-process routes do not
      work: ``llama-completion`` exits silently in b11065, and a grammar in ``llama-cli``
      collides with the chat template's control tokens.
    - Setup is in ``docs/runtime.md``.
"""

from __future__ import annotations

import json
import logging
import re
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

LOGGER = logging.getLogger(__name__)

#: Fallback only. The model is config-driven -- read it with :func:`resolve_model` so a
#: judge re-running the pipeline gets the model the report names, and so the choice is
#: recorded in `results/_manifest.json` rather than buried in source. Do not read this
#: constant directly. This is the model the 2026-09-20 spike measured.
DEFAULT_MODEL = "qwen2.5-7b-instruct-q4_k_m"

#: Default endpoint. Loopback by construction; see :func:`resolve_endpoint`.
DEFAULT_ENDPOINT = "http://127.0.0.1:8080"

#: The only hosts that keep D17's guarantee. A bind to anything else is a hosted API.
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1", "[::1]"})

REASONING_DIR = "reasoning"
RATIONALES_FILE = "rationales.json"
REASONING_FILE = "reasoning.json"

#: Evidence grades, weakest first. Channel E counts records at each grade; the strongest
#: grade present is what travels with the claim (CLAUDE.md: "the grade travels with the
#: claim into the report"). ``ungraded`` is deliberately not a grade -- it is an absence.
GRADE_COLUMNS = (("n_clinical", "clinical"), ("n_in_vivo", "in_vivo"),
                 ("n_in_vitro", "in_vitro"))
GRADE_WEIGHT = {"clinical": 1.0, "in_vivo": 0.7, "in_vitro": 0.45, "ungraded": 0.2}

#: Patterns that must never appear in a packet. The subject's phenotype is patient data
#: (CLAUDE.md constraint 7) and a disease name from channel D's evidence is a statement
#: about which conditions resemble this child.
FORBIDDEN_IN_PACKET = (
    (re.compile(r"\bHP:\d{7}\b"), "an HPO term id"),
    (re.compile(r"\b(?:MONDO|ORPHA|OMIM):\d+\b", re.I), "a disease identifier"),
    (re.compile(r"\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y):\d{6,9}\b"), "a genomic coordinate"),
    (re.compile(r"\bc\.\d+[+\-_]?\d*[ACGT]+>[ACGT]+"), "an HGVS coding change"),
)

#: What the model is allowed to return. ``evidence_grade`` is deliberately absent -- it is
#: supplied, never generated (D17 correction 3).
RUBRIC_SCHEMA = {
    "type": "object",
    "properties": {
        "mechanism_supported": {"type": "boolean"},
        "contradiction_found": {"type": "boolean"},
        "contradiction_detail": {"type": "string"},
        "citation_keys": {"type": "array", "items": {"type": "string"}, "maxItems": 6},
        "rationale": {"type": "string"},
    },
    "required": ["mechanism_supported", "contradiction_found", "contradiction_detail",
                 "citation_keys", "rationale"],
}

SYSTEM_PROMPT = (
    "You reason over supplied evidence for a drug-repurposing shortlist in a paediatric "
    "rare disease (mosaic variegated aneuploidy). You are given a candidate, the channels "
    "that ranked it, and the literature records behind it.\n"
    "Rules:\n"
    "- Cite ONLY the citation keys supplied. Never invent one.\n"
    "- Do not introduce a mechanism no supplied record states.\n"
    "- Most aneuploidy-selectivity evidence comes from transformed cancer cells. This "
    "subject's cells are not transformed. Say so when it applies.\n"
    "- This is hypothesis generation, never a clinical recommendation or an efficacy claim."
)

ADVERSARIAL_PROMPT = (
    "Argue AGAINST this candidate using only the supplied evidence. Name the strongest "
    "reason it should not be shortlisted: weak evidence grade, a contradicting record, "
    "a mechanism that does not transfer from cancer cells to non-transformed cells, or a "
    "recorded caution. If there is no such reason in the supplied evidence, say so plainly "
    "rather than inventing one.\n"
    "Put that argument in `rationale`. Set `contradiction_found` only when a supplied "
    "record actually points the other way, and put that record's substance in "
    "`contradiction_detail` -- a limitation of the evidence is not a contradiction in it."
)

CAVEATS = (
    "Candidates are in two tiers by the KIND of evidence behind them, not by strength "
    "(decision D18). A literature-tier candidate rests on a verified, graded citation a "
    "reader can check. A network-only candidate rests on a proximity score this pipeline "
    "computed, with nothing published linking the drug to this disease. Presenting them "
    "in one ranked list would invite a reader to treat those as the same kind of claim.",
    "Network-only candidates are not sent to the model at all. Asking whether absent "
    "evidence supports a mechanism yields a paragraph that reads like analysis and states "
    "only what the pipeline already knew; the absence is recorded as the finding instead.",
    "Rationales are generated by a local open-weights model over evidence the channels "
    "already produced. The model introduces no mechanism and no citation of its own: "
    "every citation key is validated against the supplied packet and dropped otherwise.",
    "Confidence is computed in Python from the supplied evidence grade, the number of "
    "supporting channels, and the model's rubric. It is a hypothesis-strength score and "
    "is never a probability of clinical benefit.",
    "The evidence grade is supplied by channel E, not judged by the model, because a "
    "model of this size graded cell-line work as clinical during the 2026-09-20 spike.",
    "An adversarial pass runs against every candidate before its confidence is finalised. "
    "Where it found a reason to object, that reason is recorded next to the rationale.",
)


def resolve_model(config: dict) -> str:
    """Return the model id for this run.

    Args:
        config: Parsed pipeline configuration; reads ``model``.

    Returns:
        The configured model id, or :data:`DEFAULT_MODEL` when unset.
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


def read_tsv(path: Path) -> list:
    """Read a tab-separated table into a list of dicts."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def strongest_grade(row: dict) -> str:
    """Return the strongest evidence grade channel E recorded for a compound.

    Args:
        row: A channel E ``candidates.tsv`` row.

    Returns:
        ``clinical``, ``in_vivo``, ``in_vitro`` or ``ungraded``.

    Note:
        Strongest rather than modal: a single clinical record is the claim a reader will
        weigh, and burying it under a count of in-vitro records would understate it. The
        counts travel into the packet too, so the reader sees the distribution as well.
    """
    for column, grade in GRADE_COLUMNS:
        try:
            if int(row.get(column) or 0) > 0:
                return grade
        except ValueError:
            continue
    return "ungraded"


def citation_index(references: dict) -> dict:
    """Map ChEMBL id to the verified citations channel E recorded for it.

    Args:
        references: Parsed ``results/l2/channel_e_prior/references.json``.

    Returns:
        ``{chembl_id: [citation string, ...]}``, verified citations only. A citation that
        failed verification is not evidence and never reaches a packet.
    """
    index: dict = {}
    for compound in references.get("compounds", []):
        keys = [c["citation"] for c in compound.get("curated_citations", [])
                if c.get("status") == "verified" and c.get("citation")]
        if keys:
            index[compound.get("chembl_id", "")] = keys
    return index


def normalise_key(key: str) -> str:
    """Strip the brackets a model puts round a citation key."""
    return key.strip().strip("[]").strip()


def validate_citations(returned, allowed) -> tuple:
    """Split returned citation keys into those that were supplied and those invented.

    Args:
        returned: Citation keys as the model emitted them.
        allowed: The keys that were in the packet.

    Returns:
        ``(kept, dropped)``, both lists, order preserved, duplicates removed.

    Note:
        This is the enforcement behind "do not introduce a mechanism no channel cited".
        A dropped key is recorded rather than silently discarded, because a model that
        invents citations is a finding about the model, not a detail to hide.
    """
    allowed_set = set(allowed)
    kept: list = []
    dropped: list = []
    for raw in returned or []:
        key = normalise_key(str(raw))
        if not key:
            continue
        target = kept if key in allowed_set else dropped
        if key not in target:
            target.append(key)
    return kept, dropped


def confidence(rubric: dict, *, grade: str, n_channels: int, n_adverse: int,
               discriminating: bool) -> float:
    """Compute a hypothesis-strength score in [0, 1].

    Args:
        rubric: The model's filled rubric.
        grade: Evidence grade supplied by channel E, not by the model.
        n_channels: How many channels ranked this candidate.
        n_adverse: Records channel E flagged as pointing the other way.
        discriminating: Whether the supporting channels are discriminating rather than
            broad (L3's convergence classification).

    Returns:
        A score in [0, 1].

    Note:
        Deliberately arithmetic and boring. The model never sees a number and never
        produces one; every input here is either a count the pipeline computed or a
        boolean the model can answer from the supplied text. That makes the score
        testable, which a model-produced number would not be.

        It is a hypothesis-strength score: how well-supported this candidate is *within
        this pipeline's evidence*. It is not a probability of clinical benefit, and
        CLAUDE.md forbids presenting it as one.
    """
    score = GRADE_WEIGHT.get(grade, GRADE_WEIGHT["ungraded"])

    if not rubric.get("mechanism_supported", False):
        score *= 0.4
    if rubric.get("contradiction_found", False):
        score *= 0.6
    if n_adverse > 0:
        score *= 0.7

    # Cross-channel agreement is the architecture's core claim, so it earns weight -- but
    # only when the agreeing channels discriminate. Agreement among broad channels is
    # agreement that a drug exists.
    if n_channels >= 2 and discriminating:
        score = min(1.0, score * 1.25)
    elif n_channels >= 2:
        score = min(1.0, score * 1.05)

    return round(max(0.0, min(1.0, score)), 3)


def assert_publishable(packet: dict) -> None:
    """Raise if a packet contains anything patient-derived.

    Args:
        packet: The evidence packet about to be serialized.

    Raises:
        ValueError: On the first forbidden pattern found.

    Note:
        Fail closed, and check the serialized form rather than the fields, so a value
        nested somewhere unexpected is still caught. This runs before every request; the
        cost is a regex sweep over a few kilobytes and the alternative is sending the
        subject's phenotype to a process that writes logs.
    """
    blob = json.dumps(packet, ensure_ascii=False)
    for pattern, description in FORBIDDEN_IN_PACKET:
        found = pattern.search(blob)
        if found:
            raise ValueError(
                f"refusing to send an evidence packet containing {description} "
                f"({found.group(0)!r}). Patient-derived fields never leave the pipeline "
                "(CLAUDE.md constraint 7, decision D17)."
            )


#: Stated verbatim as a Tier 2 row's evidence, not as a footnote (D18).
NO_LITERATURE = (
    "No published evidence links this drug to this disease. It reaches the shortlist on "
    "network proximity alone: its targets sit near the disease module in the interactome, "
    "which is a result produced by this pipeline rather than a finding reported anywhere."
)


def evidence_tier(packet: dict) -> str:
    """Classify a candidate by the *kind* of evidence behind it (decision D18).

    Args:
        packet: The evidence packet.

    Returns:
        ``literature`` when at least one verified citation with a grade supports it,
        ``network_only`` otherwise.

    Note:
        The tiers are not a ranking. They are different kinds of claim: a graded citation
        is a finding somebody published and a reader can check, while a proximity z-score
        is a number this pipeline computed. Presenting 83 candidates in one ordered list
        invites a reader to treat those as comparable, and on the first full run 82 of 83
        were network-only -- so the distinction is most of the result, not a detail.
    """
    if packet.get("citation_keys") and packet.get("evidence_grade_supplied") != "ungraded":
        return "literature"
    return "network_only"


def build_packet(row: dict, *, channel_rows: dict, citations) -> dict:
    """Assemble one candidate's evidence packet from the channel outputs.

    Args:
        row: A row of ``results/l3/candidates.tsv``.
        channel_rows: ``{channel: row}`` for the channels that ranked this candidate.
        citations: Verified citation keys for this compound.

    Returns:
        The packet, containing only publishable fields.

    Note:
        Channel D contributes its rank and nothing else. Its ``evidence.json`` states
        which diseases resemble the subject, which is a statement about the subject;
        L3 itself reads only ``candidates.tsv`` from that channel and this does the same.
    """
    prior = channel_rows.get("prior", {})
    proximity = channel_rows.get("proximity", {})

    packet = {
        "candidate": row.get("drug_name", ""),
        "drug_type": row.get("drug_type", ""),
        "clinical_stage": row.get("clinical_stage", ""),
        "pharm_class_moa": row.get("pharm_class_moa", ""),
        "aggregate": {
            "rank": row.get("rank", ""),
            "n_channels_supporting": row.get("n_channels_supporting", ""),
            "supporting_channels": row.get("supporting_channels", ""),
            "convergence": row.get("convergence", ""),
        },
        "evidence_grade_supplied": strongest_grade(prior) if prior else "ungraded",
        "citation_keys": list(citations),
    }

    if prior:
        packet["literature_prior"] = {
            "axis": prior.get("axis", ""),
            "target": prior.get("target", ""),
            "direction": prior.get("direction", ""),
            "n_clinical": prior.get("n_clinical", "0"),
            "n_in_vivo": prior.get("n_in_vivo", "0"),
            "n_in_vitro": prior.get("n_in_vitro", "0"),
            "n_adverse_direction_records": prior.get("n_adverse_direction_records", "0"),
            "caution": prior.get("caution", ""),
        }
    if proximity:
        packet["network_proximity"] = {
            "z": proximity.get("z", ""),
            "n_targets": proximity.get("n_targets", ""),
            "n_targets_in_module": proximity.get("n_targets_in_module", ""),
            "note": proximity.get("note", ""),
        }
    if "phenotype" in (row.get("supporting_channels") or ""):
        # Rank only. See the note in this function's docstring.
        packet["phenotype_channel"] = {"rank": row.get("rank_phenotype", "")}

    return packet


def probe(endpoint: str, *, timeout: int = 5) -> None:
    """Check that a model server is actually answering, before reasoning over anything.

    Args:
        endpoint: Loopback endpoint from :func:`resolve_endpoint`.
        timeout: Seconds to wait. Deliberately short -- this is a liveness check, not a
            generation call.

    Raises:
        ValueError: If the endpoint does not answer as a model server.

    Note:
        Two measured reasons this is a ``/health`` request rather than a socket connect,
        both from WSL mirrored networking on 2026-09-21:

        - Connecting to a **dead** loopback port does not refuse, it **hangs** until the
          socket timeout. Without this probe the pipeline would stall for
          ``reasoning_timeout_seconds`` on the first candidate and again on every one
          after it.
        - Connecting to some high ports **succeeds** with nothing listening. A socket
          that opens is therefore not evidence that a server exists, so the check has to
          be a request whose answer only a model server can give.
    """
    try:
        with urllib.request.urlopen(f"{endpoint}/health", timeout=timeout) as response:
            body = json.load(response)
    except (OSError, urllib.error.URLError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(
            f"no model server answering at {endpoint} ({type(error).__name__}). "
            "Start one as described in docs/runtime.md, or set l3.reasoning: false."
        ) from error
    if body.get("status") != "ok":
        raise ValueError(
            f"{endpoint} answered /health with {body!r} rather than an ok status; "
            "it does not look like a model server. See docs/runtime.md."
        )


def _post(endpoint: str, payload: dict, *, timeout: int) -> dict:
    """POST a chat-completions request and return the parsed response."""
    request = urllib.request.Request(
        f"{endpoint}/v1/chat/completions",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def ask(packet: dict, *, endpoint: str, model: str, seed: int, adversarial: bool = False,
        timeout: int = 300) -> dict:
    """Run one reasoning pass over a packet and return the parsed rubric.

    Args:
        packet: The evidence packet. Checked by :func:`assert_publishable` first.
        endpoint: Loopback endpoint from :func:`resolve_endpoint`.
        model: Model id, recorded with the result.
        seed: Run seed, for determinism.
        adversarial: Run the pass that argues against the candidate.
        timeout: Seconds to wait.

    Returns:
        The parsed rubric.

    Raises:
        ValueError: If the packet is not publishable, or the response is not valid JSON
            despite the schema. Both are recorded by the caller, never swallowed.
    """
    assert_publishable(packet)

    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    if adversarial:
        messages.append({"role": "system", "content": ADVERSARIAL_PROMPT})
    messages.append({"role": "user", "content": json.dumps(packet, ensure_ascii=False)})

    payload = {
        "messages": messages,
        "model": model,
        "temperature": 0,
        "seed": seed,
        "response_format": {
            "type": "json_schema",
            "json_schema": {"name": "rubric", "strict": True, "schema": RUBRIC_SCHEMA},
        },
    }
    response = _post(endpoint, payload, timeout=timeout)
    content = response["choices"][0]["message"]["content"]
    try:
        return json.loads(content)
    except json.JSONDecodeError as error:
        raise ValueError(f"model returned unparseable JSON despite the schema: {error}") from error


def reason_over_candidates(config: dict) -> None:
    """Synthesize rationale, contradicting evidence, and confidence per candidate.

    Args:
        config: Parsed pipeline configuration. Uses ``results_dir``, ``seed``, ``model``
            and ``reasoning_endpoint``.

    Raises:
        FileNotFoundError: If L3 has not produced a candidate table.
        ValueError: If the endpoint is not loopback, or every candidate failed.

    Note:
        Scope: L4's survivors when that table exists, otherwise the top ``reasoning_top_n``
        of L3. The layer order runs L3 before L4, so on a first full run only the top-N is
        available; re-running this step after L4 narrows it to candidates that actually
        survived triage, which is what the dossier needs. The scope used is recorded.
    """
    results_dir = Path(config["results_dir"])
    l3_dir = results_dir / "l3"
    candidates_path = l3_dir / "candidates.tsv"
    if not candidates_path.exists():
        raise FileNotFoundError(
            f"no L3 candidate table at {candidates_path}. Run l3_integrate first."
        )

    endpoint = resolve_endpoint(config)
    model = resolve_model(config)
    seed = int(config.get("seed", 42))
    settings = dict(config.get("l3") or {})
    top_n = int(settings.get("reasoning_top_n", 100))
    timeout = int(settings.get("reasoning_timeout_seconds", 300))

    # Before reading anything else: a missing server must fail in seconds, not stall the
    # run. See probe() for why this is a request rather than a connection attempt.
    probe(endpoint, timeout=int(settings.get("reasoning_probe_seconds", 5)))

    rows = read_tsv(candidates_path)

    survivors_path = results_dir / "l4" / "survivors.tsv"
    if survivors_path.exists():
        survivor_ids = {r.get("chembl_id") for r in read_tsv(survivors_path)}
        scoped = [r for r in rows if r.get("chembl_id") in survivor_ids]
        scope = "l4_survivors"
    else:
        scoped = rows[:top_n]
        scope = f"l3_top_{top_n}"
    LOGGER.info("reasoning over %d candidate(s), scope=%s, model=%s", len(scoped), scope, model)

    l2_dir = results_dir / "l2"
    channel_tables = {}
    for channel, folder in (("prior", "channel_e_prior"), ("proximity", "channel_b_proximity")):
        path = l2_dir / folder / "candidates.tsv"
        if path.exists():
            channel_tables[channel] = {r.get("chembl_id"): r for r in read_tsv(path)}

    references_path = l2_dir / "channel_e_prior" / "references.json"
    citations_by_id = (citation_index(json.loads(references_path.read_text(encoding="utf-8")))
                       if references_path.exists() else {})

    out_dir = l3_dir / REASONING_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    rationales = []
    failures = []
    for row in scoped:
        chembl_id = row.get("chembl_id", "")
        channel_rows = {c: t[chembl_id] for c, t in channel_tables.items() if chembl_id in t}
        packet = build_packet(row, channel_rows=channel_rows,
                              citations=citations_by_id.get(chembl_id, []))

        # D18: a candidate with no literature is not sent to the model. Asking it whether
        # absent evidence supports a mechanism produces a paragraph that reads like
        # analysis and says only what the pipeline already knew. The absence is the
        # finding, and it is recorded as one.
        if evidence_tier(packet) == "network_only":
            proximity = channel_rows.get("proximity", {})
            rationales.append({
                "chembl_id": chembl_id,
                "drug_name": row.get("drug_name", ""),
                "rank": row.get("rank", ""),
                "evidence_tier": "network_only",
                "rationale": NO_LITERATURE,
                "contradicting_evidence": "",
                "contradiction_found": False,
                "mechanism_supported": None,   # not assessed; absence is not a false answer
                "evidence_grade_supplied": packet["evidence_grade_supplied"],
                "confidence": None,            # see D18: not reported where it cannot discriminate
                "network_proximity": {"z": proximity.get("z", ""),
                                      "n_targets_in_module": proximity.get(
                                          "n_targets_in_module", "")},
                "citation_keys": [],
                "dropped_citation_keys": [],
                "adversarial_citation_keys": [],
                "model_consulted": False,
            })
            continue

        try:
            rubric = ask(packet, endpoint=endpoint, model=model, seed=seed, timeout=timeout)
            against = ask(packet, endpoint=endpoint, model=model, seed=seed,
                          adversarial=True, timeout=timeout)
        except (ValueError, OSError, urllib.error.URLError, KeyError) as error:
            # Never drop a candidate silently: a candidate that could not be reasoned
            # about is a gap in the dossier, and the report must be able to say so.
            failures.append({"chembl_id": chembl_id, "drug_name": row.get("drug_name", ""),
                             "error": f"{type(error).__name__}: {error}"})
            LOGGER.warning("reasoning failed for %s: %s", row.get("drug_name", chembl_id), error)
            continue

        allowed = packet["citation_keys"]
        kept, dropped = validate_citations(rubric.get("citation_keys"), allowed)
        against_kept, against_dropped = validate_citations(against.get("citation_keys"), allowed)

        prior = channel_rows.get("prior", {})
        try:
            n_adverse = int(prior.get("n_adverse_direction_records") or 0)
        except ValueError:
            n_adverse = 0

        # The adversarial pass can raise a contradiction the first pass missed; take
        # either. Confidence must never be higher for having asked the question.
        merged = dict(rubric)
        merged["contradiction_found"] = bool(
            rubric.get("contradiction_found") or against.get("contradiction_found"))

        # The adversarial pass argues its objection in `rationale`; `contradiction_detail`
        # is reserved for a supplied record that actually points the other way. Reading
        # only the latter silently discarded every objection the pass produced, which is
        # how an adversarial pass becomes decoration.
        objection = ((against.get("contradiction_detail") or "").strip()
                     or (against.get("rationale") or "").strip())

        rationales.append({
            "chembl_id": chembl_id,
            "drug_name": row.get("drug_name", ""),
            "rank": row.get("rank", ""),
            "evidence_tier": "literature",
            "model_consulted": True,
            "rationale": rubric.get("rationale", ""),
            "contradicting_evidence": objection,
            "contradiction_found": merged["contradiction_found"],
            "mechanism_supported": bool(rubric.get("mechanism_supported")),
            "evidence_grade_supplied": packet["evidence_grade_supplied"],
            "confidence": confidence(
                merged,
                grade=packet["evidence_grade_supplied"],
                n_channels=int(row.get("n_channels_supporting") or 0),
                n_adverse=n_adverse,
                discriminating=(row.get("convergence") == "discriminating"),
            ),
            "citation_keys": kept,
            "dropped_citation_keys": sorted(set(dropped) | set(against_dropped)),
            "adversarial_citation_keys": against_kept,
        })

    # Network-only rows always succeed -- they are recorded, not reasoned -- so a bare
    # "no rationales" check would never fire even if the model failed on every candidate
    # it was actually asked about. Judge on the literature tier.
    if failures and not any(r["evidence_tier"] == "literature" for r in rationales):
        raise ValueError(
            f"every literature-tier candidate failed to reason ({len(failures)} failure(s)). "
            f"The endpoint {endpoint} answered /health but produced nothing usable; "
            "see docs/runtime.md."
        )

    assessed = [r for r in rationales if r["evidence_tier"] == "literature"]

    (out_dir / RATIONALES_FILE).write_text(
        json.dumps(rationales, indent=2, ensure_ascii=False), encoding="utf-8")
    (out_dir / REASONING_FILE).write_text(json.dumps({
        "layer": "l3_integrate.reasoning",
        "model": model,
        "endpoint": endpoint,
        "seed": seed,
        "scope": scope,
        "n_candidates": len(scoped),
        "n_reasoned": len(rationales),
        "n_failed": len(failures),
        "failures": failures,
        "n_with_dropped_citations": sum(1 for r in rationales if r["dropped_citation_keys"]),
        # A field that takes one value for every candidate carries no information, and a
        # dossier column of identical answers looks like a finding. These counts make that
        # visible in the record instead of leaving a reader to notice it.
        "by_evidence_tier": {
            "literature": sum(1 for r in rationales if r["evidence_tier"] == "literature"),
            "network_only": sum(1 for r in rationales if r["evidence_tier"] == "network_only"),
        },
        "distribution": {
            "note": ("Over the literature tier only. The network-only tier is not sent to "
                     "the model (D18), so counting it here would report an absence as an "
                     "answer."),
            "n_assessed": len(assessed),
            "mechanism_supported": sum(1 for r in assessed if r["mechanism_supported"]),
            "contradiction_found": sum(1 for r in assessed if r["contradiction_found"]),
            "with_contradicting_evidence": sum(
                1 for r in assessed if r["contradicting_evidence"]),
            "with_citations": sum(1 for r in assessed if r["citation_keys"]),
            "confidence_min": min((r["confidence"] for r in assessed), default=None),
            "confidence_max": max((r["confidence"] for r in assessed), default=None),
        },
        "caveats": list(CAVEATS),
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    LOGGER.info("L3 reasoning written: %d reasoned, %d failed", len(rationales), len(failures))
