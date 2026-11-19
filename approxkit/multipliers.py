"""Approximate multipliers (unsigned n x n -> 2n bits), bit-accurate."""
import numpy as np
from .base import ApproxUnit, as_int, msb_pos
from . import cost
from .cost import GE


def _m(k):
    return (1 << k) - 1


def _pp_stats(cells, n):
    """cells: iterable of column indices of every kept partial-product bit.
    Returns (bit_count, active_column_count up to the 2n-bit output)."""
    cols = list(cells)
    if not cols:
        return 0, 0
    return len(cols), 2 * n - min(cols)


class ExactMultiplier(ApproxUnit):
    """Reference AND-array multiplier."""
    kind, family = "multiplier", "exact_mul"

    def _approx(self, a, b):
        return a * b

    def area_ge(self):
        return cost.exact_array_mult_ge(self.n)


class BrokenArray(ApproxUnit):
    """Broken-array multiplier.

    Partial-product bit (row i, col j) = b_i & a_j sits at column i+j.  Rows
    i < hbl are omitted entirely (horizontal break) and all bits in columns
    < vbl are omitted (vertical break).  The rest is summed exactly.
    """
    kind, family = "multiplier", "bam"

    def __init__(self, n, vbl, hbl=0):
        super().__init__(n)
        vbl, hbl = as_int("vbl", vbl), as_int("hbl", hbl)
        if not (0 <= vbl <= 2 * self.n and 0 <= hbl <= self.n):
            raise ValueError("0<=vbl<=2n, 0<=hbl<=n")
        self.vbl, self.hbl = vbl, hbl

    @property
    def params(self):
        return {"vbl": self.vbl, "hbl": self.hbl}

    def _approx(self, a, b):
        acc = np.zeros_like(a)
        keep = ~np.int64(_m(self.vbl))
        for i in range(self.hbl, self.n):
            acc += np.where(((b >> i) & 1) > 0, (a << i) & keep, 0)
        return acc

    def _cells(self):
        n = self.n
        return [i + j for i in range(self.hbl, n) for j in range(n) if i + j >= self.vbl]

    def area_ge(self):
        m, c = _pp_stats(self._cells(), self.n)
        return cost.carry_save_ge(m, c, GE["AND2"])


class TruncatedMultiplier(BrokenArray):
    """Fixed-width style truncation: drop all partial-product bits in
    columns < k, no error compensation."""
    family = "trunc_mul"

    def __init__(self, n, k):
        super().__init__(n, vbl=k, hbl=0)
        self.k = int(k)

    @property
    def params(self):
        return {"k": self.k}


class DRUM(ApproxUnit):
    """Dynamic-range unbiased multiplier.

    Each operand is reduced to a k-bit segment starting at its leading one;
    if bits were discarded the segment LSB is forced to 1 (unbiasing).  The
    k x k exact product is shifted back by the sum of the discarded widths.
    Operands that fit in k bits are used exactly.  Requires 2 <= k <= n.
    """
    kind, family = "multiplier", "drum"

    def __init__(self, n, k):
        super().__init__(n)
        k = as_int("k", k)
        if not (2 <= k <= self.n):
            raise ValueError("2 <= k <= n")
        self.k = k

    @property
    def params(self):
        return {"k": self.k}

    def _seg(self, x):
        la = msb_pos(x, self.n)
        sh = np.maximum(la - self.k + 1, 0)
        seg = x >> sh
        seg = np.where(sh > 0, seg | 1, seg)
        return seg, sh

    def _approx(self, a, b):
        sa, ha = self._seg(a)
        sb, hb = self._seg(b)
        return (sa * sb) << (ha + hb)

    def area_ge(self):
        n, k = self.n, self.k
        seg_sh = cost.barrel_shifter_ge(k, n - k)
        out_sh = cost.barrel_shifter_ge(2 * n, 2 * (n - k))
        return (2 * cost.lod_ge(n) + 2 * seg_sh + cost.exact_array_mult_ge(k)
                + cost.ripple_adder_ge(cost.clog2(n) + 1) + out_sh)


class Mitchell(ApproxUnit):
    """Mitchell logarithmic multiplier.

    log2(x) ~ k + f with k = leading-one position and f = fractional mantissa
    (n-1 fractional bits, left-aligned).  The two logs are added in fixed
    point and converted back with the linear antilog 2^K (1+f); the result is
    truncated to an integer (floor).  Always underestimates.  Zero in -> 0.
    """
    kind, family = "multiplier", "mitchell"

    def _approx(self, a, b):
        n = self.n
        F = n - 1
        la, lb = msb_pos(a, n), msb_pos(b, n)
        nz = (a > 0) & (b > 0)
        la_, lb_ = np.maximum(la, 0), np.maximum(lb, 0)
        fa = (a - (np.int64(1) << la_)) << (F - la_)
        fb = (b - (np.int64(1) << lb_)) << (F - lb_)
        L = ((la_ + lb_) << F) + fa + fb
        K = L >> F
        val = (np.int64(1) << F) | (L & _m(F))
        res = np.where(K >= F, val << np.maximum(K - F, 0), val >> np.maximum(F - K, 0))
        return np.where(nz, res, 0)

    def area_ge(self):
        n = self.n
        F = n - 1
        return (2 * cost.lod_ge(n) + 2 * cost.barrel_shifter_ge(F, F)
                + cost.ripple_adder_ge(cost.clog2(n) + F)
                + cost.barrel_shifter_ge(2 * n, n))


class BoothApprox(ApproxUnit):
    """Radix-4 Booth multiplier with approximate partial products.

    The unsigned multiplier b is zero-extended and recoded into D = n//2+1
    digits d_i in {-2..2}.  Partial product i = d_i * a * 4^i.  Approximation:
    the *magnitude* |d_i| * a * 4^i has its bits in columns < k cleared, then
    it is negated if d_i < 0 (sign-magnitude truncation, so the error of a
    negative partial product has the opposite sign and the rounding is not
    biased through two's-complement flooring).  The sum is taken modulo 2^(2n)
    (a result that would go negative wraps, as in hardware); k = 0 -> exact.
    """
    kind, family = "multiplier", "booth_approx"

    def __init__(self, n, k):
        super().__init__(n)
        k = as_int("k", k)
        if not (0 <= k <= 2 * self.n):
            raise ValueError("0 <= k <= 2n")
        self.k = k
        self.D = self.n // 2 + 1

    @property
    def params(self):
        return {"k": self.k}

    def digits(self, b):
        out = []
        for i in range(self.D):
            b2 = (b >> (2 * i + 1)) & 1
            b1 = (b >> (2 * i)) & 1
            b0 = ((b >> (2 * i - 1)) & 1) if i > 0 else 0
            out.append(-2 * b2 + b1 + b0)
        return out

    def _approx(self, a, b):
        full = _m(2 * self.n)
        keep = ~np.int64(_m(self.k))
        acc = np.zeros_like(a)
        for i, d in enumerate(self.digits(b)):
            mag = ((np.abs(d) * a) << (2 * i)) & keep
            acc += np.where(d >= 0, mag, -mag)
        return acc & full

    def area_ge(self):
        n = self.n
        cells = []
        for i in range(self.D):
            cells += [2 * i + j for j in range(n + 2)]   # n+1 pp bits + 1 neg/sign bit
        cells = [c for c in cells if self.k <= c < 2 * n]
        m, c = _pp_stats(cells, n)
        pp_gen = GE["MUX2"] + GE["XOR2"] - 1.0           # select + conditional invert (rounded)
        enc = self.D * 6.0                               # radix-4 encoder per digit
        return enc + cost.carry_save_ge(m, c, pp_gen)
