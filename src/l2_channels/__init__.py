"""L2 -- Channels: five parallel, deliberately independent candidate generators.

Each channel is a weak signal on its own. Consensus across channels is the actual
evidence -- see ``src.l3_integrate``. Channels must not share intermediate scores or
read each other's output, or the independence assumption behind the aggregation is void.
"""

from . import (
    channel_a_kg,
    channel_b_proximity,
    channel_c_signature,
    channel_d_phenotype,
    channel_e_prior,
)
from .run import CHANNEL_REGISTRY, run

__all__ = [
    "CHANNEL_REGISTRY",
    "channel_a_kg",
    "channel_b_proximity",
    "channel_c_signature",
    "channel_d_phenotype",
    "channel_e_prior",
    "run",
]
