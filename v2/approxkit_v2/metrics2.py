"""Error metrics for v2 units (signed multipliers, dividers, reciprocal, squarers).

e = approx - exact.  Reported over *valid* operand combinations (b != 0 for
dividers, x != 0 for reciprocal):
  ER, MED (mean |e|), MRED (mean |e|/|exact| over exact != 0), MRED16 (same over
  |exact| >= 16: removes the floor-quantisation artefact of small outputs),
  NMED = MED / max|exact|, WCE = max|e|, bias = mean e, RMSE.
Exhaustive when the operand space has <= 2^22 combinations, otherwise Monte
Carlo (uniform operands, seeded).  Monte Carlo results carry 95 % confidence
intervals from 16 independent batch means (Student t, 15 dof); exhaustive
results are exact (CI half-width 0) and WCE is exact only when exhaustive.
"""
from dataclasses import dataclass
import numpy as np

T95_15 = 2.131          # two-sided 95 % Student t, 15 degrees of freedom
EXH_LIMIT = 1 << 22


def _grid(unit):
    r = np.arange(unit.lo, unit.hi + 1, dtype=np.int64)
    if unit.inputs == 1:
        return r, None
    return np.repeat(r, r.size), np.tile(r, r.size)


def _stats(unit, a, b):
    v = unit.valid(a, b)
    a = a[v]
    b = b[v] if b is not None else None
    ex = unit.exact(a, b)
    e = unit(a, b) - ex
    ae = np.abs(e).astype(np.float64)
    aex = np.abs(ex).astype(np.float64)
    m1, m16 = aex >= 1, aex >= 16
    return dict(n=e.size, nz=int(np.count_nonzero(e)), sae=ae.sum(), se=float(e.sum()),
                se2=float((e.astype(np.float64) ** 2).sum()), wce=float(ae.max()) if e.size else 0.0,
                srel=float((ae[m1] / aex[m1]).sum()), nrel=int(m1.sum()),
                srel16=float((ae[m16] / aex[m16]).sum()), nrel16=int(m16.sum()))


def _fin(s, unit):
    return dict(ER=s["nz"] / s["n"], MED=s["sae"] / s["n"], MRED=s["srel"] / max(s["nrel"], 1),
                MRED16=s["srel16"] / max(s["nrel16"], 1) if s["nrel16"] else float("nan"),
                NMED=s["sae"] / s["n"] / unit.max_abs_exact, WCE=s["wce"],
                bias=s["se"] / s["n"], RMSE=(s["se2"] / s["n"]) ** 0.5)


def evaluate(unit, samples=1 << 20, seed=0, mode="auto"):
    """Return a flat dict: metrics, `<metric>_ci` (95 % half-width), mode, samples."""
    total = (unit.hi - unit.lo + 1) ** unit.inputs
    if mode == "auto":
        mode = "exhaustive" if total <= EXH_LIMIT else "monte_carlo"
    keys = ("ER", "MED", "MRED", "MRED16", "NMED", "bias", "RMSE")
    if mode == "exhaustive":
        a, b = _grid(unit)
        s = _stats(unit, a, b)
        out = _fin(s, unit)
        out.update({k + "_ci": 0.0 for k in keys})
        out.update(mode=mode, samples=s["n"])
    else:
        rng = np.random.default_rng(seed)
        nb, per = 16, samples // 16
        tot = None
        fins = []
        for _ in range(nb):
            a = rng.integers(unit.lo, unit.hi + 1, per, dtype=np.int64)
            b = rng.integers(unit.lo, unit.hi + 1, per, dtype=np.int64) if unit.inputs == 2 else None
            s = _stats(unit, a, b)
            fins.append(_fin(s, unit))
            tot = s if tot is None else {k: (max(tot[k], s[k]) if k == "wce" else tot[k] + s[k]) for k in tot}
        out = _fin(tot, unit)
        for k in keys:
            xs = np.array([f[k] for f in fins], dtype=np.float64)
            out[k + "_ci"] = float(T95_15 * np.nanstd(xs, ddof=1) / np.sqrt(len(xs)))
        out.update(mode=mode, samples=tot["n"])
    out.update(unit=unit.name, family=unit.family, n=unit.n)
    return out
