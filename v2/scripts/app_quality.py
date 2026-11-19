#!/usr/bin/env python3
"""Application-quality harness for the 8-bit SIGNED multipliers of the v2 catalog.

  * image convolution (3x3 Gaussian blur Q7 [power-of-two taps], 'blur_odd' [taps 5/21/24, none a power of two], sharpen Q4, Sobel-x) on 8 SYNTHETIC 96x96 images
    (gradients+discs, stripes, texture, text-like blocks; generated here, no external data);
    pixels are mapped to signed int8 (p-128), kernels are signed int8 -> every tap goes through the unit under test.
    Quality is measured against the output of the EXACT multiplier: PSNR (dB, capped 99) and SSIM.
  * dot product (256 terms, int8 x int8) and 31-tap FIR (v1 low-pass), output SNR vs exact.
Accumulation is exact (this harness isolates the multiplier).  CI = 95 % Student-t over images / batches.
  python scripts/app_quality.py --out results   -> results/app_quality_v2.csv
"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from approxkit_v2.catalog import catalog
from approxkit_v2.quality import psnr, ssim, snr_db, mean_ci

KERNELS = {  # (coeffs, shift, offset-mode)
    "blur": (np.array([[8, 16, 8], [16, 32, 16], [8, 16, 8]]), 7, "pix"),
    "blur_odd": (np.array([[5, 21, 5], [21, 24, 21], [5, 21, 5]]), 7, "pix"),   # no power-of-two taps except none: stresses log/DRUM units
    "sharpen": (np.array([[0, -16, 0], [-16, 80, -16], [0, -16, 0]]), 4, "pix"),
    "sobelx": (np.array([[-32, 0, 32], [-64, 0, 64], [-32, 0, 32]]), 5, "abs"),
}


def images(size=96):
    out = []
    yy, xx = np.mgrid[0:size, 0:size]
    for s in range(2):
        r = np.random.default_rng(100 + s)
        g = 40 + 170 * (xx + yy) / (2 * size)
        for _ in range(5):
            cx, cy, rad = r.integers(10, size - 10, 2).tolist() + [int(r.integers(6, 18))]
            g = np.where((xx - cx) ** 2 + (yy - cy) ** 2 < rad * rad, r.integers(20, 240), g)
        out.append(("discs%d" % s, g))
        st = 128 + 100 * np.sin(2 * np.pi * (xx * (0.05 + 0.03 * s) + yy * 0.02))
        out.append(("stripes%d" % s, st))
        tex = 128 + 25 * r.normal(size=(size, size)) + 60 * np.sin(xx / 7.0 + s) * np.cos(yy / 9.0)
        out.append(("texture%d" % s, tex))
        blk = np.full((size, size), 200.0)
        for _ in range(40):
            x0, y0 = r.integers(0, size - 8, 2); w, h = r.integers(2, 8, 2)
            blk[y0:y0 + h, x0:x0 + w] = r.integers(10, 100)
        out.append(("blocks%d" % s, blk))
    return [(n, np.clip(np.round(a), 0, 255).astype(np.int64)) for n, a in out]


def conv(unit, img, kern, shift, mode):
    s = img - 128                               # signed int8 pixels
    H, W = s.shape
    acc = np.zeros((H - 2, W - 2), dtype=np.int64)
    for dy in range(3):
        for dx in range(3):
            c = int(kern[dy, dx])
            if c == 0:
                continue
            acc += unit(s[dy:dy + H - 2, dx:dx + W - 2], np.full((H - 2, W - 2), c, dtype=np.int64))
    v = np.abs(acc) >> shift if mode == "abs" else (acc >> shift) + 128
    return np.clip(v, 0, 255)


def fir_taps():
    T, fc = 31, 0.15
    t = np.arange(T) - (T - 1) / 2
    h = np.sinc(2 * fc * t) * np.hamming(T)
    return np.round(h / h.sum() * 127 * 0.9).astype(np.int64)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ents = [e for e in catalog(8) if e.group == "smul"]
    imgs = images()
    exact = ents[0].unit()
    assert ents[0].key == "s_exact"
    ref = {(n, k): conv(exact, im, *KERNELS[k]) for n, im in imgs for k in KERNELS}
    rng = np.random.default_rng(5)
    h = fir_taps()
    firx = [np.clip(np.round(60 * np.sin(2 * np.pi * 0.05 * np.arange(1500) + p) + 30 * np.sin(2 * np.pi * 0.31 * np.arange(1500))
                             + rng.normal(0, 8, 1500)), -127, 127).astype(np.int64) for p in range(8)]
    dx = [rng.integers(-127, 128, (32, 256)) for _ in range(8)]
    dy = [rng.integers(-127, 128, (32, 256)) for _ in range(8)]

    def firp(u, x):
        n = len(x) - len(h) + 1
        acc = np.zeros(n, dtype=np.int64)
        for t in range(len(h)):
            acc += u(x[len(h) - 1 - t:len(h) - 1 - t + n], np.full(n, h[t]))
        return acc

    fref = [firp(exact, x) for x in firx]
    dref = [(x * y).sum(1) for x, y in zip(dx, dy)]
    rows = []
    for e in ents:
        u = e.unit()
        r = dict(key=e.key, role=e.role, family=e.family)
        for k in KERNELS:
            ps = [psnr(ref[(n, k)], conv(u, im, *KERNELS[k])) for n, im in imgs]
            ss = [ssim(ref[(n, k)], conv(u, im, *KERNELS[k])) for n, im in imgs]
            r[f"{k}_psnr"], r[f"{k}_psnr_ci"] = mean_ci(ps)
            r[f"{k}_ssim"], r[f"{k}_ssim_ci"] = mean_ci(ss)
        fs = [snr_db(fr, firp(u, x)) for fr, x in zip(fref, firx)]
        r["fir_snr"], r["fir_snr_ci"] = mean_ci(fs)
        ds = [snr_db(dr, np.array([u(x[i], y[i]).sum() for i in range(32)])) for dr, x, y in zip(dref, dx, dy)]
        r["dot_snr"], r["dot_snr_ci"] = mean_ci(ds)
        r["img_psnr_mean"] = float(np.mean([r[f"{k}_psnr"] for k in KERNELS]))
        r["img_ssim_mean"] = float(np.mean([r[f"{k}_ssim"] for k in KERNELS]))
        rows.append(r)
        print(f"{e.key:20s} blurPSNR={r['blur_psnr']:6.2f} sharpSSIM={r['sharpen_ssim']:.4f} sobelPSNR={r['sobelx_psnr']:6.2f} fir={r['fir_snr']:6.1f} dot={r['dot_snr']:6.1f}", flush=True)
    with open(os.path.join(a.out, "app_quality_v2.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
