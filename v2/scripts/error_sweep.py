#!/usr/bin/env python3
"""Error metrics for every v2 catalog design (exhaustive at n=8: 65 536 pairs / 256 inputs).
  python scripts/error_sweep.py --n 8 --out results   ->  results/errors_v2_n8.csv"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from approxkit_v2.catalog import catalog
from approxkit_v2.metrics2 import evaluate

ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=8); ap.add_argument("--out", default="results")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = []
for e in catalog(a.n):
    m = evaluate(e.unit())
    rows.append(dict(key=e.key, group=e.group, role=e.role, family=e.family, **{k: m[k] for k in
                ("ER", "MED", "MRED", "MRED16", "NMED", "WCE", "bias", "RMSE", "mode", "samples")}))
with open(os.path.join(a.out, f"errors_v2_n{a.n}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "designs")
for g in ("smul", "div", "recip", "sq"):
    print(f"-- {g}")
    for r in rows:
        if r["group"] == g and r["role"] == "approx" and (r["key"].count("_") < 9):
            pass
    for r in sorted((r for r in rows if r["group"] == g), key=lambda r: r["NMED"])[:: max(1, len([1 for r in rows if r['group']==g]) // 8)]:
        print(f"  {r['key']:26s} ER={r['ER']:.3f} MRED16={r['MRED16']:.4f} NMED={r['NMED']:.2e} WCE={r['WCE']:.0f} bias={r['bias']:.1f}")
