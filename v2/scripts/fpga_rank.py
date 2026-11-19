#!/usr/bin/env python3
"""Rank agreement between the sky130 cost ranking (v1 flow) and GENERIC-vendor-flow FPGA cost estimates.
NOT PolarFire results: yosys synth_ice40 / synth_ecp5 / synth_xilinx(xc7) of the installed yosys 0.33 give
LUT4/LUT6/carry counts for *other* architectures; PolarFire (LUT4 + carry, different mapping) was never run.
  python scripts/fpga_rank.py --res results    -> results/rank_corr.csv (+ prints)
Reference ranking for v1 designs = sky130hd_um2 column of v1's results/rtl_synth_n8.csv (v1 flow);
the same designs re-synthesised here (costs_v1_n8.csv) are used for the FPGA columns.  Spearman rho and
Kendall tau-b with 95 % percentile-bootstrap CI (2000 resamples over designs, seed 0)."""
import argparse, csv, os, sys
import numpy as np
from scipy.stats import spearmanr, kendalltau

ap = argparse.ArgumentParser(); ap.add_argument("--res", default="results"); ap.add_argument("--v1", default="../../results/rtl_synth_n8.csv")
a = ap.parse_args()


def load(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def load_costs(res, src):
    """final CSV if the synthesis run completed, else the incrementally written .partial.csv (JSON lines)"""
    import json
    f = os.path.join(res, f"costs_{src}_n8.csv")
    if os.path.exists(f):
        return load(f)
    return [json.loads(l) for l in open(os.path.join(res, f"costs_{src}_n8.partial.csv"))]


def fl(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return np.nan


def corr(x, y, nb=2000, seed=0):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    if len(x) < 4:
        return None
    rho, tau = spearmanr(x, y)[0], kendalltau(x, y)[0]
    rng = np.random.default_rng(seed); rs, ts = [], []
    for _ in range(nb):
        i = rng.integers(0, len(x), len(x))
        if np.ptp(x[i]) == 0 or np.ptp(y[i]) == 0:
            continue
        rs.append(spearmanr(x[i], y[i])[0]); ts.append(kendalltau(x[i], y[i])[0])
    return dict(n=len(x), rho=rho, rho_lo=np.percentile(rs, 2.5), rho_hi=np.percentile(rs, 97.5),
                tau=tau, tau_lo=np.percentile(ts, 2.5), tau_hi=np.percentile(ts, 97.5))


FP = ["ice40_lut4", "ice40nc_lut4", "ecp5_lut4", "xil_lut", "gen_ge"]
out = []
v1 = [r for r in load_costs(a.res, "v1") if not r.get("synth_error")]
v1ref = {r["case"]: float(r["sky130hd_um2"]) for r in load(os.path.join(os.path.dirname(os.path.abspath(__file__)), a.v1))}
for r in v1:
    r["sky_ref"] = v1ref.get(r["key"], np.nan)
print("v1 keys matched to v1 sky130 reference:", sum(np.isfinite([fl(r, 'sky_ref') for r in v1])), "of", len(v1))
v2 = [r for r in load_costs(a.res, "v2") if not r.get("synth_error")]
for r in v2:
    r["sky_ref"] = fl(r, "sky130_um2")
sets = [("v1 all", v1), ("v1 multipliers", [r for r in v1 if "mult" in r["module"]]), ("v1 adders", [r for r in v1 if "mult" not in r["module"]]),
        ("v2 all", v2)] + [(f"v2 {g}", [r for r in v2 if r["group"] == g]) for g in ("smul", "div", "sq", "recip")]
for name, rs in sets:
    # sanity: re-synthesised sky130 vs v1 reference (same flow -> should agree)
    if name.startswith("v1"):
        c = corr([fl(r, "sky130_um2") for r in rs], [r["sky_ref"] for r in rs], nb=200)
        if c:
            out.append(dict(set=name, metric="sky130_resynth_vs_v1_ref", **c))
    for m in FP:
        c = corr([r["sky_ref"] for r in rs], [fl(r, m) for r in rs])
        if c:
            out.append(dict(set=name, metric=m, **c))
with open(os.path.join(a.res, "rank_corr.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for r in out:
    print(f"{r['set']:16s} {r['metric']:26s} n={r['n']:3d} rho={r['rho']:.3f} [{r['rho_lo']:.3f},{r['rho_hi']:.3f}] tau={r['tau']:.3f} [{r['tau_lo']:.3f},{r['tau_hi']:.3f}]")
