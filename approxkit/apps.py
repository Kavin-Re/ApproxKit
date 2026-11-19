"""Application demos: dot product, FIR filter and matrix multiply built from
approximate multipliers (and, optionally, an approximate accumulator adder).

Operands are signed integers; magnitudes go through the *unsigned* multiplier
and the sign is applied afterwards.  Accumulation uses two unsigned
accumulators (positive / negative terms), each through the supplied adder of
width `acc_bits`; their difference is the result.  With adder=None
accumulation is exact.
"""
import numpy as np


def signed_mul(mult, x, y):
    x = np.asarray(x, dtype=np.int64)
    y = np.asarray(y, dtype=np.int64)
    mag = mult(np.abs(x), np.abs(y))
    return np.where((x < 0) ^ (y < 0), -mag, mag)


class _Acc:
    def __init__(self, shape, adder):
        self.adder = adder
        self.pos = np.zeros(shape, dtype=np.int64)
        self.neg = np.zeros(shape, dtype=np.int64)
        self.mask = (1 << adder.n) - 1 if adder is not None else None

    def add(self, term):
        p = np.where(term > 0, term, 0)
        q = np.where(term < 0, -term, 0)
        if self.adder is None:
            self.pos += p
            self.neg += q
        else:
            self.pos = self.adder(self.pos, p) & self.mask
            self.neg = self.adder(self.neg, q) & self.mask

    def value(self):
        return self.pos - self.neg


def dot(x, y, mult, adder=None):
    """Dot product along the last axis of x, y (same shape)."""
    x, y = np.asarray(x), np.asarray(y)
    acc = _Acc(x.shape[:-1], adder)
    for t in range(x.shape[-1]):
        acc.add(signed_mul(mult, x[..., t], y[..., t]))
    return acc.value()


def fir(x, h, mult, adder=None):
    """Direct-form FIR, 'valid' outputs: y[i] = sum_t h[t] * x[i+len(h)-1-t]."""
    x = np.asarray(x, dtype=np.int64)
    h = np.asarray(h, dtype=np.int64)
    T = len(h)
    out_n = len(x) - T + 1
    acc = _Acc((out_n,), adder)
    for t in range(T):
        acc.add(signed_mul(mult, x[T - 1 - t:T - 1 - t + out_n], np.full(out_n, h[t])))
    return acc.value()


def matmul(A, B, mult, adder=None):
    A = np.asarray(A, dtype=np.int64)
    B = np.asarray(B, dtype=np.int64)
    acc = _Acc((A.shape[0], B.shape[1]), adder)
    for k in range(A.shape[1]):
        acc.add(signed_mul(mult, A[:, k][:, None], B[k, :][None, :]))
    return acc.value()


def accuracy(approx, exact):
    """SNR in dB (signal = exact), relative RMS error, max abs error."""
    approx = np.asarray(approx, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    err = approx - exact
    pe = float(np.mean(err ** 2))
    ps = float(np.mean(exact ** 2))
    snr = float("inf") if pe == 0 else 10 * np.log10(ps / pe)
    return {"snr_db": snr, "rel_rms": (pe / ps) ** 0.5 if ps else float("nan"),
            "max_abs_err": float(np.max(np.abs(err)))}
