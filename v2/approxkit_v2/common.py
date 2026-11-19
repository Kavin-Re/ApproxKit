"""Shared plumbing for ApproxKit v2: path hook to v1, wrap helper, Unit base.

v2 never modifies v1; it imports v1 (`approxkit`) from the parent directory.
"""
import os, sys
import numpy as np

V2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
V1_ROOT = os.path.abspath(os.path.join(V2_ROOT, ".."))
if V1_ROOT not in sys.path:
    sys.path.insert(0, V1_ROOT)
import approxkit as ak            # noqa: E402  (v1, read-only)
from approxkit.base import msb_pos  # noqa: E402,F401


def wrap_signed(v, bits):
    """Two's-complement wrap of integers to `bits` bits (returns signed int64)."""
    m = 1 << bits
    v = np.asarray(v, dtype=np.int64) & (m - 1)
    return np.where(v >= (m >> 1), v - m, v)


def mask_bits(v, bits):
    """Raw bit pattern (non-negative) of v in `bits` bits."""
    return np.asarray(v, dtype=np.int64) & ((1 << bits) - 1)


class Unit:
    """Interface shared by all v2 units (signed multipliers, dividers, ...).

    kind    : 'smul' | 'udiv' | 'urecip' | 'usq'
    inputs  : 1 or 2 operands
    lo, hi  : inclusive operand range
    out_bits: width of the raw output word (two's complement if signed)
    """
    kind = None
    family = "base"
    signed = False
    inputs = 2

    def __init__(self, n):
        self.n = int(n)

    @property
    def params(self):
        return {}

    @property
    def name(self):
        ps = ",".join(f"{k}={v}" for k, v in self.params.items())
        return f"{self.family}(n={self.n}{',' + ps if ps else ''})"

    def __repr__(self):
        return self.name

    @property
    def lo(self):
        return -(1 << (self.n - 1)) if self.signed else 0

    @property
    def hi(self):
        return (1 << (self.n - 1)) - 1 if self.signed else (1 << self.n) - 1

    @property
    def out_bits(self):
        return 2 * self.n

    def valid(self, a, b=None):
        return np.ones(np.shape(a), dtype=bool)

    def exact(self, a, b=None):  # pragma: no cover
        raise NotImplementedError

    @property
    def max_abs_exact(self):  # pragma: no cover
        raise NotImplementedError

    def _check(self, a, b):
        for t in (a, b):
            if t is not None and t.size and (t.min() < self.lo or t.max() > self.hi):
                raise ValueError(f"{self.name}: operand outside [{self.lo},{self.hi}]")

    def __call__(self, a, b=None):
        a = np.array(a, dtype=np.int64)
        if self.inputs == 2:
            b = np.array(b, dtype=np.int64)
            a, b = (np.array(t) for t in np.broadcast_arrays(a, b))
        else:
            b = None
        self._check(a, b)
        return self._approx(a, b)

    def _approx(self, a, b):  # pragma: no cover
        raise NotImplementedError
