"""Application demos: accuracy vs approximation (and cost proxy).
Usage: python examples/demo_apps.py [--out results]"""
import argparse, os, sys, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import approxkit as ak
from approxkit import apps, dse
from approxkit.plotting import _style, FAMILY_STYLE, COST_LABEL
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results")
out = ap.parse_args().out
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(42)
mults = dse.default_designs("multiplier", 8)

# ---- 1. FIR: 31-tap windowed-sinc low-pass, Q7 coefficients, 8-bit signed samples
T, fc = 31, 0.15
t = np.arange(T) - (T - 1) / 2
h = np.sinc(2 * fc * t) * np.hamming(T)
h = np.round(h / h.sum() * 127 * 0.9).astype(np.int64)
ns = np.arange(4000)
sig = 60 * np.sin(2 * np.pi * 0.05 * ns) + 30 * np.sin(2 * np.pi * 0.31 * ns) + rng.normal(0, 8, ns.size)
x = np.clip(np.round(sig), -127, 127).astype(np.int64)
ref = apps.fir(x, h, ak.ExactMultiplier(8))
assert np.array_equal(ref, np.convolve(x, h, mode="valid"))

# ---- 2. dot product (correlator-like): 256-long, 64 random vector pairs
xd = rng.integers(-127, 128, (64, 256)); yd = rng.integers(-127, 128, (64, 256))
refd = apps.dot(xd, yd, ak.ExactMultiplier(8))

# ---- 3. matmul 16x16 uint8
A = rng.integers(0, 256, (16, 16)); B = rng.integers(0, 256, (16, 16))
refm = A @ B

rows = []
for u in mults:
    f = apps.accuracy(apps.fir(x, h, u), ref)
    d = apps.accuracy(apps.dot(xd, yd, u), refd)
    m = apps.accuracy(apps.matmul(A, B, u), refm)
    rows.append(dict(unit=u.name, family=u.family, area_ge=u.area_ge(), fir_snr_db=f["snr_db"],
                     dot_rel_rms=d["rel_rms"], matmul_rel_rms=m["rel_rms"], matmul_snr_db=m["snr_db"]))
with open(f"{out}/apps_multipliers.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ---- 4. approximate accumulator adder (exact multiplier), 24-bit accumulators
adders = [ak.ExactAdder(24)] + [c(24, k) for k in (4, 8, 12) for c in (ak.LOA, ak.TruncatedAdder, ak.ETA)] + \
         [ak.GeAr(24, 4, 8), ak.GeAr(24, 4, 12)]
arows = []
for ad in adders:
    f = apps.accuracy(apps.fir(x, h, ak.ExactMultiplier(8), ad), ref)
    d = apps.accuracy(apps.dot(xd, yd, ak.ExactMultiplier(8), ad), refd)
    arows.append(dict(unit=ad.name, area_ge=ad.area_ge(), fir_snr_db=f["snr_db"], dot_rel_rms=d["rel_rms"]))
with open(f"{out}/apps_accum_adders.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(arows[0])); w.writeheader(); w.writerows(arows)

# ---- plots
_style()
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
for ax, key, lab, yl in [(axs[0], "fir_snr_db", "FIR output SNR (dB), higher is better", False),
                         (axs[1], "dot_rel_rms", "dot-product relative RMS error", True),
                         (axs[2], "matmul_rel_rms", "16x16 matmul relative RMS error", True)]:
    for fam in dict.fromkeys(r["family"] for r in rows):
        col, mk = FAMILY_STYLE.get(fam, ("#777", "o"))
        pts = [r for r in rows if r["family"] == fam and np.isfinite(r[key]) and (r[key] > 0 or not yl)]
        ax.scatter([p["area_ge"] for p in pts], [p[key] for p in pts], c=col, marker=mk, s=22,
                   alpha=0.8, linewidths=0, label=fam)
    if yl:
        ax.set_yscale("log")
    ax.set_xlabel("cost proxy (GE) - not synthesized area"); ax.set_ylabel(lab)
    if key == "fir_snr_db":
        ax.set_ylim(0, 80)
axs[0].legend(ncol=2, fontsize=7)
fig.suptitle("8-bit approximate multipliers in applications (exact multipliers: SNR = inf, plotted off-scale/omitted)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{out}/apps_multipliers.png"); plt.close(fig)

print("FIR 31 taps, 8-bit: best SNR / cheapest design per family (area_ge, SNR dB, dot rel-RMS, matmul rel-RMS)")
for fam in dict.fromkeys(r["family"] for r in rows):
    fr = [r for r in rows if r["family"] == fam]
    picks = [min(fr, key=lambda r: r["area_ge"]), max(fr, key=lambda r: (r["fir_snr_db"], -r["area_ge"]))]
    for r in dict.fromkeys(r["unit"] for r in picks) and {r["unit"]: r for r in picks}.values():
        print(f"  {r['unit']:42s} GE={r['area_ge']:6.1f} SNR={r['fir_snr_db']:6.1f} dot={r['dot_rel_rms']:.4f} mm={r['matmul_rel_rms']:.4f}")
print("\nApproximate accumulator adders (24-bit, exact multiplier):")
for r in arows:
    print(f"  {r['unit']:38s} GE={r['area_ge']:6.1f} SNR={r['fir_snr_db']:6.1f} dB dot rel-RMS={r['dot_rel_rms']:.5f}")
