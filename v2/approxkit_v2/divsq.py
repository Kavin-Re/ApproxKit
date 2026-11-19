"""Approximate dividers, reciprocal and squarers (unsigned, bit-accurate).

Division convention (all dividers): quotient = floor(a / b) on n bits; b = 0
returns all ones (2^n - 1) in every variant, exact included.  Metrics ignore b = 0.
Reciprocal: x in [1, 2^n - 1] -> approx of floor(2^m / x) on m+1 bits; x = 0 -> all ones.
"""
import numpy as np
from .common import Unit, msb_pos


def _ones(bits):
    return (1 << bits) - 1


# ----------------------------------------------------------------- dividers
class _Div(Unit):
    kind = "udiv"

    @property
    def out_bits(self):
        return self.n

    @property
    def max_abs_exact(self):
        return _ones(self.n)

    def valid(self, a, b=None):
        return b != 0

    def exact(self, a, b=None):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        return np.where(b == 0, _ones(self.n), a // np.maximum(b, 1))


class ExactDivider(_Div):
    family = "div_exact"

    def _approx(self, a, b):
        return self.exact(a, b)


class TruncRestoringDiv(_Div):
    """Restoring divider with the last t iterations removed.

    Only the top n-t dividend bits are shifted into the partial remainder, so
    q = floor(floor(a / 2^t) / b) * 2^t  (quotient quantised to multiples of 2^t).
    comp=1 sets quotient bit t-1 (mid-point of the quantisation step).
    t = 0 is the exact divider.  Saves t of the n subtract/compare rows.
    """
    family = "div_trunc"

    def __init__(self, n, t, comp=0):
        super().__init__(n)
        if not (0 <= t <= n - 1):
            raise ValueError("0 <= t <= n-1")
        self.t, self.comp = int(t), int(comp)

    @property
    def params(self):
        return {"t": self.t, "comp": self.comp}

    def _approx(self, a, b):
        q = ((a >> self.t) // np.maximum(b, 1)) << self.t
        if self.comp and self.t:
            q = q | (1 << (self.t - 1))
        return np.where(b == 0, _ones(self.n), q)


class MitchellDiv(_Div):
    """Mitchell logarithmic divider.

    log2 x ~ k + f (F = n-1 fraction bits).  D = (ka - kb) * 2^F + fa - fb is a
    signed fixed-point number; K = floor(D / 2^F), frac = D mod 2^F; the linear
    antilog gives (2^F | frac) * 2^(K-F), floored to an integer.  a = 0 -> 0.
    """
    family = "div_mitchell"

    def _approx(self, a, b):
        n = self.n
        F = n - 1
        la, lb = np.maximum(msb_pos(a, n), 0), np.maximum(msb_pos(b, n), 0)
        fa = (a - (np.int64(1) << la)) << (F - la)
        fb = (b - (np.int64(1) << lb)) << (F - lb)
        D = ((la - lb) << F) + fa - fb
        K = D >> F
        val = (np.int64(1) << F) | (D & _ones(F))
        sh_r = np.clip(F - K, 0, 62)
        res = np.where(K >= F, val << np.clip(K - F, 0, 62), val >> sh_r)
        res = np.where(a == 0, 0, res)
        return np.where(b == 0, _ones(n), res)


class DRUMDiv(_Div):
    """Segment (DRUM-style) divider: k-bit segments of a and b taken from their
    leading ones (LSB forced to 1 when bits were dropped), Qf = floor(sa * 2^k / sb)
    (a 2k-by-k exact division), then shifted by ha - hb - k and saturated to n bits.
    My own construction in the spirit of DRUM, not a published divider.  2 <= k <= n.
    """
    family = "div_drum"

    def __init__(self, n, k):
        super().__init__(n)
        if not (2 <= k <= n):
            raise ValueError("2 <= k <= n")
        self.k = int(k)

    @property
    def params(self):
        return {"k": self.k}

    def _seg(self, x):
        l = np.maximum(msb_pos(x, self.n), 0)
        h = np.maximum(l - self.k + 1, 0)
        s = x >> h
        return np.where(h > 0, s | 1, s), h

    def _approx(self, a, b):
        k, n = self.k, self.n
        sa, ha = self._seg(a)
        sb, hb = self._seg(b)
        qf = (sa << k) // np.maximum(sb, 1)
        sh = ha - hb - k
        q = np.where(sh >= 0, qf << np.clip(sh, 0, 62), qf >> np.clip(-sh, 0, 62))
        q = np.minimum(q, _ones(n))
        q = np.where(a == 0, 0, q)
        return np.where(b == 0, _ones(n), q)


# --------------------------------------------------------------- reciprocal
class _Recip(Unit):
    kind, inputs = "urecip", 1

    def __init__(self, n, m=None):
        super().__init__(n)
        self.m = int(m if m is not None else n)

    @property
    def params(self):
        return {"m": self.m}

    @property
    def lo(self):
        return 0

    @property
    def out_bits(self):
        return self.m + 1

    @property
    def max_abs_exact(self):
        return 1 << self.m

    def valid(self, a, b=None):
        return a != 0

    def exact(self, a, b=None):
        a = np.asarray(a, dtype=np.int64)
        return np.where(a == 0, _ones(self.m + 1), (np.int64(1) << self.m) // np.maximum(a, 1))


class ExactReciprocal(_Recip):
    family = "recip_exact"

    def _approx(self, a, b):
        return self.exact(a)


class MitchellRecip(_Recip):
    """x = 2^k (1+f):  R = floor( 2^(m-k-1) * (2 - f) ) = floor((3*2^k - x) * 2^(m-2k-1)).
    One leading-one detector, one subtract, one shifter; no adder of logs."""
    family = "recip_mitchell"

    def _approx(self, a, b):
        k = np.maximum(msb_pos(a, self.n), 0)
        v = (np.int64(3) << k) - a
        s = self.m - 2 * k - 1
        r = np.where(s >= 0, v << np.clip(s, 0, 62), v >> np.clip(-s, 0, 62))
        return np.where(a == 0, _ones(self.m + 1), r)


# ----------------------------------------------------------------- squarers
class _Sq(Unit):
    kind, inputs = "usq", 1

    @property
    def lo(self):
        return 0

    @property
    def max_abs_exact(self):
        return _ones(self.n) ** 2

    def exact(self, a, b=None):
        a = np.asarray(a, dtype=np.int64)
        return a * a


class ExactSquarer(_Sq):
    family = "sq_exact"

    def _approx(self, a, b):
        return a * a


class TruncSquarer(_Sq):
    """Folded squarer matrix: a^2 = sum_i a_i 2^(2i) + sum_{i<j} a_i a_j 2^(i+j+1).
    Bits in columns < k dropped; comp=1 adds round(E[dropped]) for uniform inputs."""
    family = "sq_trunc"

    def __init__(self, n, k, comp=0):
        super().__init__(n)
        if not (0 <= k <= 2 * n):
            raise ValueError("0 <= k <= 2n")
        self.k, self.comp = int(k), int(comp)

    @property
    def params(self):
        return {"k": self.k, "comp": self.comp}

    def comp_value(self):
        n, k = self.n, self.k
        e4 = 0
        for i in range(n):
            if 2 * i < k:
                e4 += 2 << (2 * i)
            for j in range(i + 1, n):
                if i + j + 1 < k:
                    e4 += 1 << (i + j + 1)
        return (e4 + 2) // 4

    def _approx(self, a, b):
        n, k = self.n, self.k
        acc = np.zeros_like(a)
        for i in range(n):
            ai = (a >> i) & 1
            if 2 * i >= k:
                acc += ai << (2 * i)
            for j in range(i + 1, n):
                if i + j + 1 >= k:
                    acc += (ai & ((a >> j) & 1)) << (i + j + 1)
        if self.comp and k:
            acc += self.comp_value()
        return acc & _ones(2 * n)


class DRUMSquarer(_Sq):
    """k-bit segment from the leading one (LSB forced to 1 if bits dropped),
    exact k-bit square, shifted by 2*h.  2 <= k <= n."""
    family = "sq_drum"

    def __init__(self, n, k):
        super().__init__(n)
        if not (2 <= k <= n):
            raise ValueError("2 <= k <= n")
        self.k = int(k)

    @property
    def params(self):
        return {"k": self.k}

    def _approx(self, a, b):
        l = np.maximum(msb_pos(a, self.n), 0)
        h = np.maximum(l - self.k + 1, 0)
        s = a >> h
        s = np.where(h > 0, s | 1, s)
        return (s * s) << (2 * h)


class MitchellSquarer(_Sq):
    """2*log2 x ~ 2(k+f): the antilog of the doubled log; no adder at all."""
    family = "sq_mitchell"

    def _approx(self, a, b):
        n = self.n
        F = n - 1
        l = np.maximum(msb_pos(a, n), 0)
        f = (a - (np.int64(1) << l)) << (F - l)
        L = ((l << F) + f) << 1
        K = L >> F
        val = (np.int64(1) << F) | (L & _ones(F))
        r = np.where(K >= F, val << np.clip(K - F, 0, 62), val >> np.clip(F - K, 0, 62))
        return np.where(a == 0, 0, r)
