"""Signed (two's-complement) approximate multipliers, n x n -> 2n bits.

Operands in [-2^(n-1), 2^(n-1)-1].  Results are wrapped to 2n-bit two's
complement exactly like the Verilog (no unit here ever wraps for n = 8; tests
assert this).

Units
  SignMag(unit_factory, n, sat)   sign-magnitude wrapper around any v1 UNSIGNED unit
  BaughWooley(n, k, mode, comp)   native: Baugh-Wooley partial-product matrix
        mode 'trunc' : partial-product bits in columns < k dropped (BW constant
                       bits are kept); comp=1 adds a constant equal to the rounded
                       expected value of the dropped bits (uniform inputs)
        mode 'orlow' : columns < k are reduced by OR (no carries) instead of being
                       added; columns >= k exact          (k <= n)
  SignedDRUM(n, k)               native: DRUM on two's-complement segments
"""
import numpy as np
from .common import Unit, wrap_signed, msb_pos, ak


def _bw_rows(a, b, n):
    """Baugh-Wooley partial-product rows. Yields (j, row) with row an int64 whose
    bit (i + j) is the pp bit a_i b_j (NAND'ed for exactly one of i, j = n-1),
    i.e. already placed in its column."""
    for j in range(n):
        bj = (b >> j) & 1
        row = np.zeros_like(a)
        for i in range(n):
            t = ((a >> i) & 1) & bj
            if (i == n - 1) != (j == n - 1):
                t = 1 - t
            row |= t << (i + j)
        yield j, row


def bw_const(n):
    return (1 << n) + (1 << (2 * n - 1))


def bw_expected_dropped_q4(n, k):
    """4 x E[value of the BW pp bits in columns < k] under uniform inputs (integer)."""
    e4 = 0
    for i in range(n):
        for j in range(n):
            if i + j < k:
                nand = (i == n - 1) != (j == n - 1)
                e4 += (3 if nand else 1) << (i + j)
    return e4


class SignedMultiplier(Unit):
    kind, signed = "smul", True

    @property
    def max_abs_exact(self):
        return 1 << (2 * self.n - 2)

    def exact(self, a, b=None):
        return np.asarray(a, dtype=np.int64) * np.asarray(b, dtype=np.int64)


class SignMag(SignedMultiplier):
    """Sign-magnitude wrapper around an unsigned v1 unit.

    |a|, |b| (<= 2^(n-1)) go through the unsigned unit; the sign is applied to
    the magnitude afterwards.  sat=0: the unit has width n (exact corner case
    -2^(n-1)).  sat=1: the unit has width n-1 and |x| is saturated to
    2^(n-1)-1 (so -2^(n-1) behaves as -(2^(n-1)-1); saves area, adds error at
    that single operand value).  The magnitude is clamped to 2^(2n-1)-1 before
    the sign is applied (only reachable by DRUM-2/3 style overshoot).
    """
    family = "smag"

    def __init__(self, factory, n, sat=False, label=""):
        super().__init__(n)
        self.sat = bool(sat)
        self.w = n - 1 if sat else n
        self.unit = factory(self.w)
        import re
        self.label = label or re.sub(r"n=\d+,?", "", self.unit.name).replace("()", "")

    @property
    def params(self):
        return {"unit": self.label, "sat": int(self.sat)}

    def _approx(self, a, b):
        ma, mb = np.abs(a), np.abs(b)
        if self.sat:
            lim = (1 << (self.n - 1)) - 1
            ma, mb = np.minimum(ma, lim), np.minimum(mb, lim)
        mag = np.minimum(self.unit(ma, mb), (1 << (2 * self.n - 1)) - 1)
        return np.where((a < 0) ^ (b < 0), -mag, mag)


class BaughWooley(SignedMultiplier):
    family = "bw"

    def __init__(self, n, k=0, mode="trunc", comp=0):
        super().__init__(n)
        if mode not in ("trunc", "orlow"):
            raise ValueError(mode)
        if not (0 <= k <= (n if mode == "orlow" else 2 * n)):
            raise ValueError("k out of range")
        self.k, self.mode, self.comp = int(k), mode, int(comp)

    @property
    def params(self):
        p = {"k": self.k, "mode": self.mode}
        if self.mode == "trunc":
            p["comp"] = self.comp
        return p

    def _approx(self, a, b):
        n, k = self.n, self.k
        keep = ~np.int64((1 << k) - 1)
        low = np.int64((1 << k) - 1)
        acc = np.zeros_like(a)
        orl = np.zeros_like(a)
        for _, row in _bw_rows(a, b, n):
            acc += row & keep
            orl |= row & low
        acc += bw_const(n)
        if self.mode == "trunc":
            if self.comp and k:
                acc += (bw_expected_dropped_q4(n, k) + 2) // 4
            res = acc
        else:
            res = acc | orl
        return wrap_signed(res, 2 * n)


class SignedDRUM(SignedMultiplier):
    """DRUM on two's-complement operands (no sign-magnitude conversion).

    p = index of the highest bit that differs from the sign bit (-1 for 0/-1).
    The operand needs p+2 bits (sign included).  A k-bit window (sign included)
    is kept: s = max(p+2-k, 0); seg = x >>> s (arithmetic); if s > 0 the LSB of
    seg is forced to 1 (DRUM unbiasing).  Product = seg_a * seg_b * 2^(s_a+s_b),
    wrapped to 2n bits.  2 <= k <= n; operands that fit in k bits are exact.
    """
    family = "sdrum"

    def __init__(self, n, k):
        super().__init__(n)
        if not (2 <= k <= n):
            raise ValueError("2 <= k <= n")
        self.k = int(k)

    @property
    def params(self):
        return {"k": self.k}

    def _seg(self, x):
        y = np.where(x < 0, ~x, x)
        p = msb_pos(y, self.n)
        s = np.maximum(p + 2 - self.k, 0)
        seg = x >> s
        seg = np.where(s > 0, seg | 1, seg)
        return seg, s

    def _approx(self, a, b):
        sa, ha = self._seg(a)
        sb, hb = self._seg(b)
        return wrap_signed(sa * sb * (np.int64(1) << (ha + hb)), 2 * self.n)
