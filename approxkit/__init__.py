"""ApproxKit: bit-accurate approximate arithmetic units, error metrics and
a gate-equivalent design-space explorer.  (Unsigned units; see README.)"""
from .base import ApproxUnit
from .adders import ExactAdder, LOA, TruncatedAdder, ETA, GeAr, AMA
from .multipliers import (ExactMultiplier, BrokenArray, TruncatedMultiplier,
                          DRUM, Mitchell, BoothApprox)
from .metrics import evaluate, ErrorMetrics

__version__ = "0.1.0"

ADDERS = {"exact": ExactAdder, "loa": LOA, "trunc": TruncatedAdder,
          "eta": ETA, "gear": GeAr, "ama": AMA}
MULTIPLIERS = {"exact": ExactMultiplier, "bam": BrokenArray,
               "trunc": TruncatedMultiplier, "drum": DRUM,
               "mitchell": Mitchell, "booth": BoothApprox}


def make(kind, family, n, **params):
    """make('adder','loa',8,k=3) -> LOA(n=8,k=3)."""
    table = {"adder": ADDERS, "multiplier": MULTIPLIERS}[kind]
    return table[family](n, **params)
