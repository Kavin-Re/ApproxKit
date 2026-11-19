#!/usr/bin/env python3
"""Pareto plots: application quality vs hardware cost (sky130 area, and generic-flow iCE40 LUT4 count - NOT PolarFire).
  python scripts/pareto_plots.py --res results     -> results/pareto_*.png"""
import argparse, csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(); ap.add_argument("--res", default="results"); a = ap.parse_args()
rd = lambda p: list(csv.DictReader(open(os.path.join(a.res, p))))
import json
cf = os.path.join(a.res, "costs_v2_n8.csv")
crows = rd("costs_v2_n8.csv") if os.path.exists(cf) else [json.loads(l) for l in open(os.path.join(a.res, "costs_v2_n8.partial.csv"))]
cost = {r["key"]: r for r in crows if not r.get("synth_error")}
app = rd("app_quality_v2.csv")
COL = {"bw": "#1b6ca8", "sdrum": "#d1495b", "smag": "#2e8b57", "s_exact": "#222222"}
MK = {"bw": "o", "sdrum": "s", "smag": "^", "s_exact": "*"}


def front(pts):   # pts: (cost, quality[higher better]); returns indices on Pareto front (min cost, max quality)
    idx = sorted(range(len(pts)), key=lambda i: (pts[i][0], -pts[i][1])); best = -1e18; out = []
    for i in idx:
        if pts[i][1] > best + 1e-12:
            out.append(i); best = pts[i][1]
    return out


def panel(ax, ckey, ccol, qkey, ql, cap=None):
    rows = [r for r in app if r["key"] in cost and cost[r["key"]].get(ccol) not in (None, "", "None")]
    pts = [(float(cost[r["key"]][ccol]), min(float(r[qkey]), cap) if cap else float(r[qkey])) for r in rows]
    for fam in COL:
        ii = [i for i, r in enumerate(rows) if r["family"] == fam]
        ax.scatter([pts[i][0] for i in ii], [pts[i][1] for i in ii], c=COL[fam], marker=MK[fam], s=30, label=fam, alpha=.85, linewidths=0)
    f = front(pts)
    ax.plot([pts[i][0] for i in f], [pts[i][1] for i in f], "-", c="#888", lw=1, zorder=0)
    for i in f:
        ax.annotate(rows[i]["key"], pts[i], fontsize=5.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(ckey); ax.set_ylabel(ql)
    return [rows[i]["key"] for i in f]


fronts = {}
fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
for j, (q, ql, cap) in enumerate([("blur_odd_psnr", "blur_odd PSNR (dB, cap 60)", 60), ("sharpen_ssim", "sharpen SSIM", None), ("dot_snr", "dot-product SNR (dB, cap 60)", 60)]):
    fronts[q + "|sky130"] = panel(axs[0, j], "sky130 area (um^2, v1 flow)", "sky130_um2", q, ql, cap)
    fronts[q + "|ice40"] = panel(axs[1, j], "iCE40 LUT4 (generic yosys flow; NOT PolarFire)", "ice40_lut4", q, ql, cap)
axs[0, 0].legend(fontsize=7)
fig.suptitle("ApproxKit v2: 8-bit signed multipliers, application quality vs cost (line = Pareto front; mean over 8 synthetic images / 8 batches)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(a.res, "pareto_apps_smul.png"), dpi=110); plt.close(fig)

# error vs cost for div / sq
er = {r["key"]: r for r in rd("errors_v2_n8.csv")}
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, g, t in zip(axs, ("div", "sq", "smul"), ("dividers", "squarers (+recip)", "signed multipliers")):
    for k, r in er.items():
        if (r["group"] == g or (g == "sq" and r["group"] == "recip")) and k in cost and cost[k]["sky130_um2"]:
            ax.scatter(float(cost[k]["sky130_um2"]), max(float(r["NMED"]), 1e-7), c="#1b6ca8" if r["role"] == "approx" else "#222", s=22)
            if r["role"] == "exact" or g != "smul":
                ax.annotate(k.replace("div_", "").replace("sq_", ""), (float(cost[k]["sky130_um2"]), max(float(r["NMED"]), 1e-7)), fontsize=5, xytext=(2, 2), textcoords="offset points")
    ax.set_yscale("log"); ax.set_title(t, fontsize=9); ax.set_xlabel("sky130 area (um^2)"); ax.set_ylabel("NMED (exhaustive, n=8; 1e-7 floor = exact)")
fig.tight_layout(); fig.savefig(os.path.join(a.res, "pareto_error_cost.png"), dpi=110)
import json; json.dump(fronts, open(os.path.join(a.res, "pareto_fronts.json"), "w"), indent=1)
for k, v in fronts.items():
    print(k, "->", ", ".join(v))
