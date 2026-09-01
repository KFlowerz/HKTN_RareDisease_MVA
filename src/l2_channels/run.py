"""L2 -- Channel fan-out.

Purpose
    Dispatch to each candidate-generation channel enabled in
    ``config["channels"]`` and collect their ranked candidate lists for L3.

Inputs
    L1 output (disease module, upstream/downstream target sets) and
    ``config["channels"]`` -- a mapping of ``{kg, proximity, signature, phenotype,
    prior} -> bool``.

Outputs
    One ranked candidate table per enabled channel, written under
    ``config["results_dir"]/l2/<channel>/``, each row carrying the channel's own score,
    its rank, and per-candidate provenance.

Guardrail
    Channels run **independently**. No channel may read another's output or scores --
    L3's rank aggregation assumes independent evidence, and leaking scores between
    channels would manufacture false consensus. A channel that fails must be recorded as
    *failed*, never silently treated as "returned no candidates": an empty list is a
    scientific claim, an exception is not.
"""

from __future__ import annotations

from . import (
    channel_a_kg,
    channel_b_proximity,
    channel_c_signature,
    channel_d_phenotype,
    channel_e_prior,
)

#: Maps the config key under ``channels:`` to the module implementing that channel.
CHANNEL_REGISTRY = {
    "kg": channel_a_kg,
    "proximity": channel_b_proximity,
    "signature": channel_c_signature,
    "phenotype": channel_d_phenotype,
    "prior": channel_e_prior,
}


def run(config: dict) -> None:
    """Execute every enabled L2 channel.

    Args:
        config: Parsed pipeline configuration. Uses ``channels``, ``results_dir``,
            ``seed``, and the L1 artifacts.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    # TODO: for each key in `config["channels"]` that is truthy, look the module up in
    # CHANNEL_REGISTRY and call its `generate(config)`. Run them independently (a thread
    # or process pool is fine -- they share no state), record per-channel wall time and
    # failure status in a manifest, and write each channel's ranked table to
    # `results_dir/l2/<channel>/`. Raise on an unknown channel key rather than ignoring
    # it, so a typo in the config cannot silently disable a line of evidence.
    raise NotImplementedError("l2_channels.run is a scaffold stub")
