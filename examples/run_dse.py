"""Design-space exploration: error vs gate-equivalent cost PROXY.
Usage: python examples/run_dse.py [--n 8] [--out results] [--err MRED]
n <= 12 is exhaustive, wider is Monte Carlo."""
import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from approxkit import dse, evaluate
from approxkit.plotting import pareto_plot, histogram_plot

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=8)
ap.add_argument("--out", default="results")
ap.add_argument("--err", default="MRED", choices=["MRED", "NMED", "ER"])
ap.add_argument("--samples", type=int, default=1 << 20)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

for kind in ("adder", "multiplier"):
    rows = dse.explore(dse.default_designs(kind, a.n), samples=a.samples)
    front = dse.pareto_front(rows, "area_ge", a.err)
    dse.write_csv(rows, f"{a.out}/dse_{kind}_n{a.n}.csv")
    pareto_plot(rows, front, f"{a.out}/pareto_{kind}_n{a.n}_{a.err}.png",
                f"{kind}s, n={a.n}: {a.err} vs cost proxy ({rows[0]['mode']})", err_key=a.err)
    print(f"\n== {kind}s n={a.n}: {len(rows)} designs, Pareto front on (area_ge, {a.err}) ==")
    print(f"{'design':44s} {'GE':>7s} {a.err:>9s} {'NMED':>9s} {'ER':>6s} {'WCE':>8s}")
    for r in front:
        print(f"{r['unit']:44s} {r['area_ge']:7.1f} {r[a.err]:9.5f} {r['NMED']:9.5f} {r['ER']:6.3f} {r['WCE']:8d}")
    pick = [r for r in front if r["family"] != "exact_add" and r["family"] != "exact_mul"][:6]
    histogram_plot([r["_metrics"] for r in pick], f"{a.out}/errhist_{kind}_n{a.n}.png",
                   f"{kind}s n={a.n}: signed error distribution of low-cost Pareto designs")
