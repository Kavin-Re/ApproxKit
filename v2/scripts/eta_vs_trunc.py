#!/usr/bin/env python3
"""Audit observation: in the v1 dot-product / FIR demo ETA and truncated-adder accumulators
give IDENTICAL errors (apps_accum_adders.csv: dot_rel_rms 0.001423987 for both at k=4, 0.026805747 at k=8).

Mechanism (tested in tests/test_v2.py::test_eta_equals_truncation_accumulator_mechanism):
  v1 accumulates positive and negative products in two UNSIGNED accumulators and subtracts them.
  ETA: once acc bit k-1 and the addend bit k-1 are both 1 (prob ~1/4 per add), all low-k result bits are
  forced to 1, and the all-ones state is ABSORBING (a_j = 1 -> bit = 1 whether b_j is 0 or 1).  So after a
  handful of additions both accumulators have low k bits == 2^k - 1, which cancel in pos - neg.
  Truncation: both low parts are 0, cancel as well.  The upper parts are the same ((a>>k)+(b>>k)).
  => identical outputs.  Not an ETA property of interest: an artefact of the accumulate-in-two-halves organisation.

Corrected demo (this script): (1) identical-fraction vs number of accumulated terms (saturation curve);
(2) the same three adders in a SINGLE two's-complement accumulator (what a real MAC does), where ETA, LOA
and truncation DO differ; (3) both organisations, rel-RMS with 95 % CI over independent trials.
  python scripts/eta_vs_trunc.py --out results
"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from approxkit_v2.common import ak

W = 24
MASK = (1 << W) - 1
T95 = {  # two-sided 95 % t for dof 1..30 (approx. beyond)
}


def tcrit(d):
    from math import sqrt
    return 1.96 + 2.4 / max(d, 1) + 2.0 / max(d, 1) ** 2 if d < 200 else 1.96


def products(rng, ntrial, nterm):
    x = rng.integers(-127, 128, (ntrial, nterm))
    y = rng.integers(-127, 128, (ntrial, nterm))
    return x * y


def acc_two(prod, adder):
    """v1 organisation: positive / negative unsigned accumulators, difference at the end."""
    pos = np.zeros(prod.shape[0], dtype=np.int64); neg = pos.copy()
    for t in range(prod.shape[1]):
        p = prod[:, t]
        pos = (adder(pos, np.where(p > 0, p, 0)) if adder else pos + np.where(p > 0, p, 0)) & MASK
        neg = (adder(neg, np.where(p < 0, -p, 0)) if adder else neg + np.where(p < 0, -p, 0)) & MASK
    return pos, neg


def acc_single(prod, adder):
    """single 24-bit two's-complement accumulator (what a hardware MAC does)."""
    s = np.zeros(prod.shape[0], dtype=np.int64)
    for t in range(prod.shape[1]):
        s = (adder(s, prod[:, t] & MASK) if adder else s + prod[:, t]) & MASK
    return np.where(s >> (W - 1), s - (1 << W), s)


def run_mechanism(k=4, nterm=256, ntrial=200, seed=1):
    rng = np.random.default_rng(seed)
    eta, tr = ak.ETA(W, k), ak.TruncatedAdder(W, k)
    prod = products(rng, ntrial, nterm)
    pe, ne = acc_two(prod, eta)
    pt, nt = acc_two(prod, tr)
    low = (1 << k) - 1
    res = dict(frac_both_low_saturated=float(np.mean(((pe & low) == low) & ((ne & low) == low))),
               frac_identical_diff=float(np.mean((pe - ne) == (pt - nt))))
    ps, ns = acc_two(prod[:, :2], eta); qs, ms = acc_two(prod[:, :2], tr)
    res["frac_identical_diff_short"] = float(np.mean((ps - ns) == (qs - ms)))
    se, st = acc_single(prod, eta), acc_single(prod, tr)
    res["single_acc_differs"] = float(np.mean(se != st))
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(7)
    # (1) saturation curve
    rows1 = []
    prod = products(rng, 400, 64)
    for k in (4, 8, 12):
        eta, tr = ak.ETA(W, k), ak.TruncatedAdder(W, k)
        for nt in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64):
            pe, ne = acc_two(prod[:, :nt], eta); pt, nt_ = acc_two(prod[:, :nt], tr)
            rows1.append(dict(k=k, nterm=nt, frac_identical=float(np.mean((pe - ne) == (pt - nt_)))))
    with open(os.path.join(a.out, "eta_vs_trunc_saturation.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows1[0])); w.writeheader(); w.writerows(rows1)
    # (2)+(3) error in both organisations, dot length 256, 20 independent batches of 50 vectors
    rows = []
    for k in (4, 8, 12):
        units = {"exact": None, "loa": ak.LOA(W, k), "trunc": ak.TruncatedAdder(W, k), "eta": ak.ETA(W, k)}
        for org, fn in (("two_acc(v1)", acc_two), ("single_acc", acc_single)):
            per = {u: [] for u in units}
            for b in range(20):
                prod = products(rng, 50, 256)
                ref = prod.sum(1)
                for u, ad in units.items():
                    if u == "exact":
                        continue
                    r = (lambda p, n: p - n)(*fn(prod, ad)) if org.startswith("two") else fn(prod, ad)
                    per[u].append(float(np.sqrt(np.mean((r - ref) ** 2) / np.mean(ref.astype(float) ** 2))))
            for u in ("loa", "trunc", "eta"):
                v = np.array(per[u]); ci = tcrit(19) * v.std(ddof=1) / np.sqrt(len(v))
                rows.append(dict(k=k, organisation=org, adder=u, rel_rms=v.mean(), ci95=ci, batches=len(v)))
    with open(os.path.join(a.out, "eta_vs_trunc_demo.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("mechanism (k=4, 256 terms):", run_mechanism())
    print("saturation curve (k=4):", [(r["nterm"], round(r["frac_identical"], 2)) for r in rows1 if r["k"] == 4])
    print(f"{'k':>3} {'organisation':14s} {'adder':6s} {'rel_rms':>10s} {'+-95%':>9s}")
    for r in rows:
        print(f"{r['k']:3d} {r['organisation']:14s} {r['adder']:6s} {r['rel_rms']:10.5f} {r['ci95']:9.5f}")


if __name__ == "__main__":
    main()
