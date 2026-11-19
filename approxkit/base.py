"""Common base class for approximate arithmetic units.

All units are *unsigned*, bit-accurate and vectorised over numpy int64 arrays.
Adders take two n-bit operands and return an (n+1)-bit sum.  Multipliers take
two n-bit operands and return a 2n-bit product.
"""
import numpy as np


def msb_pos(x, n):
    """Index of the most significant set bit of each element (-1 for 0)."""
    x = np.asarray(x, dtype=np.int64)
    pos = np.full(x.shape, -1, dtype=np.int64)
    for j in range(n):
        pos = np.where((x >> j) > 0, j, pos)
    return pos


class ApproxUnit:
    kind = None          # "adder" | "multiplier"
    family = "base"
    MAX_N = 31

    def __init__(self, n):
        n = int(n)
        if not (2 <= n <= self.MAX_N):
            raise ValueError(f"{type(self).__name__}: n must be in [2,{self.MAX_N}]")
        self.n = n

    # ---- identification -------------------------------------------------
    @property
    def params(self):
        return {}

    @property
    def name(self):
        ps = ",".join(f"{k}={v}" for k, v in self.params.items())
        return f"{self.family}(n={self.n}{',' + ps if ps else ''})"

    def __repr__(self):
        return self.name

    # ---- behaviour ------------------------------------------------------
    @property
    def out_bits(self):
        return self.n + 1 if self.kind == "adder" else 2 * self.n

    @property
    def max_exact(self):
        m = (1 << self.n) - 1
        return 2 * m if self.kind == "adder" else m * m

    def exact(self, a, b):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        return a + b if self.kind == "adder" else a * b

    def __call__(self, a, b):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        a, b = (np.array(t) for t in np.broadcast_arrays(a, b))   # writable, same shape
        lim = 1 << self.n
        if a.size and (a.min() < 0 or a.max() >= lim or b.min() < 0 or b.max() >= lim):
            raise ValueError(f"{self.name}: operands must be unsigned {self.n}-bit")
        return self._approx(a, b)

    def _approx(self, a, b):  # pragma: no cover
        raise NotImplementedError

    # ---- cost proxies (see cost.py) ------------------------------------
    def area_ge(self):
        raise NotImplementedError

    def delay_levels(self):
        """Longest carry/propagate chain in bit positions (adders only)."""
        return None
