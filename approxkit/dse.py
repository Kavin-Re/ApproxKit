"""Design-space exploration: sweep unit parameters, compute error metrics and
the GE cost PROXY (see cost.py -- not synthesized area), extract Pareto fronts.
"""
import csv
from . import adders as A, multipliers as M
from .metrics import evaluate


def default_designs(kind, n):
    """Parameter sweep for every family at width n (always includes exact)."""
    d = []
    if kind == "adder":
        d.append(A.ExactAdder(n))
        step = max(1, n // 12)
        ks = list(range(1, n, step))
        for k in ks:
            d += [A.LOA(n, k), A.TruncatedAdder(n, k), A.ETA(n, k),
                  A.AMA(n, k, "sum_b"), A.AMA(n, k, "inv_carry")]
        for R in sorted({1, 2, 3, 4, max(1, n // 4)}):
            for P in range(1, n - R):
                d.append(A.GeAr(n, R, P))
    elif kind == "multiplier":
        d.append(M.ExactMultiplier(n))
        step = max(1, n // 8)
        for k in range(step, 2 * n - 1, step):
            d += [M.TruncatedMultiplier(n, k), M.BoothApprox(n, k)]
        for vbl in range(0, 2 * n - 2, 2 * step):
            for hbl in range(0, n // 2 + 1, max(1, n // 8)):
                if vbl or hbl:
                    d.append(M.BrokenArray(n, vbl, hbl))
        for k in range(2, n):
            d.append(M.DRUM(n, k))
        d.append(M.Mitchell(n))
    else:
        raise ValueError(kind)
    return d


def explore(units, mode="auto", samples=1 << 20, seed=0):
    """Evaluate every unit; returns list of dict rows (metrics + cost)."""
    rows = []
    for u in units:
        m = evaluate(u, mode=mode, samples=samples, seed=seed)
        r = m.row()
        r.update(family=u.family, area_ge=u.area_ge(), delay=u.delay_levels(), params=str(u.params))
        r["_metrics"] = m
        rows.append(r)
    return rows


def pareto_front(rows, cost_key="area_ge", err_key="MRED"):
    """Non-dominated rows (minimise both). Ties on both are all kept out
    except the first; result sorted by increasing cost."""
    pts = sorted(rows, key=lambda r: (r[cost_key], r[err_key]))
    front, best = [], float("inf")
    for r in pts:
        if r[err_key] < best:
            front.append(r)
            best = r[err_key]
    return front


def write_csv(rows, path):
    keys = ["unit", "family", "params", "n", "mode", "samples", "ER", "MED", "MRED", "NMED",
            "WCE", "bias", "RMSE", "area_ge", "delay"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
