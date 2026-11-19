"""Approximate adders (unsigned, n-bit operands, (n+1)-bit result).

Every class documents the exact bit-level behaviour; the Verilog in rtl/
implements the same definition and is checked against these models.
Approximation parameter k (or R,P) must leave at least one exact upper bit
(0 <= k <= n-1); k = 0 is the exact adder where applicable.
"""
import numpy as np
from .base import ApproxUnit, as_int
from . import cost
from .cost import GE


def _m(k):
    return (1 << k) - 1


class ExactAdder(ApproxUnit):
    """Reference ripple-carry adder."""
    kind, family = "adder", "exact_add"

    def _approx(self, a, b):
        return a + b

    def area_ge(self):
        return cost.ripple_adder_ge(self.n)

    def delay_levels(self):
        return self.n


class _LowerK(ApproxUnit):
    kind = "adder"

    def __init__(self, n, k):
        super().__init__(n)
        k = as_int("k", k)
        if not (0 <= k <= self.n - 1):
            raise ValueError("k must satisfy 0 <= k <= n-1")
        self.k = k

    @property
    def params(self):
        return {"k": self.k}


class LOA(_LowerK):
    """Lower-part-OR adder.

    Bits [0,k): sum bit = a_i OR b_i.  A single AND of bits k-1 of both
    operands is the carry-in of the exact upper (n-k)-bit adder.
    """
    family = "loa"

    def _approx(self, a, b):
        k = self.k
        if k == 0:
            return a + b
        lo = (a | b) & _m(k)
        cin = (a >> (k - 1)) & (b >> (k - 1)) & 1
        hi = (a >> k) + (b >> k) + cin
        return (hi << k) | lo

    def area_ge(self):
        k = self.k
        return k * GE["OR2"] + (GE["AND2"] if k else 0) + cost.ripple_adder_ge(self.n - k, carry_in=k > 0)

    def delay_levels(self):
        return self.n - self.k + (1 if self.k else 0)


class TruncatedAdder(_LowerK):
    """Lower k result bits forced to 0 and no carry from the lower part."""
    family = "trunc_add"

    def _approx(self, a, b):
        k = self.k
        return ((a >> k) + (b >> k)) << k

    def area_ge(self):
        return cost.ripple_adder_ge(self.n - self.k)

    def delay_levels(self):
        return self.n - self.k


class ETA(_LowerK):
    """Error-tolerant-adder-style (ETA-I flavour).

    Upper n-k bits: exact adder with no carry-in.  Lower k bits are scanned
    from bit k-1 downwards: sum_i = a_i XOR b_i until the first position
    where a_i = b_i = 1; that bit and every lower bit are set to 1.  No carry
    crosses the k boundary.
    """
    family = "eta"

    def _approx(self, a, b):
        k = self.k
        low = np.zeros_like(a)
        flag = np.zeros_like(a)
        for j in range(k - 1, -1, -1):
            aj, bj = (a >> j) & 1, (b >> j) & 1
            both = aj & bj
            bit = np.where((flag | both) > 0, 1, aj ^ bj)
            flag = flag | both
            low |= bit << j
        return (((a >> k) + (b >> k)) << k) | low

    def area_ge(self):
        # per low bit: XOR (sum), AND (both), OR (chain), OR (force-to-1)
        per = GE["XOR2"] + GE["AND2"] + 2 * GE["OR2"]
        return self.k * per + cost.ripple_adder_ge(self.n - self.k)

    def delay_levels(self):
        return max(self.n - self.k, self.k)


class GeAr(ApproxUnit):
    """Generic accuracy-configurable adder with block-wise carry speculation.

    Sub-adders have length L = R + P: R fresh result bits and P preceding
    'prediction' bits whose carry is speculated (carry-in = 0).  Block 0
    covers bits [0,L); block j covers [s_j, s_j+L) with s_j = min(j*R, n-L);
    block j contributes result bits above those already produced; the last
    block also produces the carry-out (bit n).  R=1 is the ACA-style adder
    (every sum bit from a window of P+1 bits).  L = n gives an exact adder.
    """
    kind, family = "adder", "gear"

    def __init__(self, n, R, P):
        super().__init__(n)
        R, P = as_int("R", R), as_int("P", P)
        n = self.n
        if R < 1 or P < 0 or R + P > n:
            raise ValueError("need R>=1, P>=0, R+P<=n")
        self.R, self.P, self.L = R, P, R + P
        starts = [0]
        while starts[-1] + self.L < n:
            starts.append(min(starts[-1] + R, n - self.L))
        self.starts = starts

    @property
    def params(self):
        return {"R": self.R, "P": self.P}

    def _approx(self, a, b):
        L, n = self.L, self.n
        out = np.zeros_like(a)
        top = 0
        for idx, s in enumerate(self.starts):
            t = ((a >> s) & _m(L)) + ((b >> s) & _m(L))   # L+1 bits
            last = idx == len(self.starts) - 1
            nb = (s + L - top) + (1 if last else 0)
            out |= ((t >> (top - s)) & _m(nb)) << top
            top = s + L
        return out

    def area_ge(self):
        return len(self.starts) * cost.ripple_adder_ge(self.L)

    def delay_levels(self):
        return self.L


class AMA(_LowerK):
    """Simplified-full-adder lower part ('AMA-style').

    The lower k bits are a ripple chain of an approximate FA cell, the exact
    upper part receives the carry of bit k-1.  Cells (NOT identical to any
    specific published cell; same spirit as the mirror-adder approximations):
      cell="sum_b":    Sum = b,           Cout = a        (wires only)
      cell="inv_carry": Cout = maj(a,b,cin), Sum = NOT Cout
    """
    family = "ama"

    def __init__(self, n, k, cell="sum_b"):
        super().__init__(n, k)
        if cell not in ("sum_b", "inv_carry"):
            raise ValueError("cell must be 'sum_b' or 'inv_carry'")
        self.cell = cell

    @property
    def params(self):
        return {"k": self.k, "cell": self.cell}

    @property
    def cell_id(self):
        return 0 if self.cell == "sum_b" else 1

    def _approx(self, a, b):
        k = self.k
        low = np.zeros_like(a)
        c = np.zeros_like(a)
        for j in range(k):
            aj, bj = (a >> j) & 1, (b >> j) & 1
            if self.cell == "sum_b":
                s, c = bj, aj
            else:
                c = (aj & bj) | (c & (aj | bj))
                s = 1 - c
            low |= s << j
        return (((a >> k) + (b >> k) + c) << k) | low

    def area_ge(self):
        per = 0.0 if self.cell == "sum_b" else GE["MAJ"] + GE["INV"]
        return self.k * per + cost.ripple_adder_ge(self.n - self.k, carry_in=self.k > 0)

    def delay_levels(self):
        return (self.n - self.k) + (0 if self.cell == "sum_b" else self.k)
