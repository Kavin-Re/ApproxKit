"""Application-quality metrics implemented from scratch (numpy only): PSNR, SSIM
(Wang et al. 2004: 11x11 Gaussian window, sigma 1.5, K1=0.01, K2=0.03, 'valid' region,
mean SSIM), output SNR.  Images are 2-D arrays with values 0..255."""
import numpy as np


def psnr(ref, test, peak=255.0, cap=99.0):
    mse = np.mean((np.asarray(ref, float) - np.asarray(test, float)) ** 2)
    return cap if mse == 0 else min(cap, 10 * np.log10(peak * peak / mse))


def _gauss(size=11, sigma=1.5):
    x = np.arange(size) - (size - 1) / 2
    g = np.exp(-x * x / (2 * sigma * sigma))
    return g / g.sum()


def _filt(img, g):
    """separable 'valid' correlation with 1-D kernel g"""
    k = g.size
    h = sum(g[i] * img[:, i:img.shape[1] - k + 1 + i] for i in range(k))
    return sum(g[i] * h[i:h.shape[0] - k + 1 + i, :] for i in range(k))


def ssim(ref, test, L=255.0, size=11, sigma=1.5, K1=0.01, K2=0.03):
    x, y = np.asarray(ref, float), np.asarray(test, float)
    g = _gauss(size, sigma)
    mx, my = _filt(x, g), _filt(y, g)
    sxx = _filt(x * x, g) - mx * mx
    syy = _filt(y * y, g) - my * my
    sxy = _filt(x * y, g) - mx * my
    C1, C2 = (K1 * L) ** 2, (K2 * L) ** 2
    s = ((2 * mx * my + C1) * (2 * sxy + C2)) / ((mx * mx + my * my + C1) * (sxx + syy + C2))
    return float(s.mean())


def snr_db(exact, approx, cap=99.0):
    e = np.asarray(exact, float); d = np.asarray(approx, float) - e
    pe = np.mean(d * d)
    return cap if pe == 0 else min(cap, 10 * np.log10(np.mean(e * e) / pe))


def t95(dof):
    tab = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
           11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131, 19: 2.093, 29: 2.045}
    if dof in tab:
        return tab[dof]
    ks = sorted(tab)
    for lo, hi in zip(ks, ks[1:]):
        if lo < dof < hi:
            return tab[lo] + (tab[hi] - tab[lo]) * (dof - lo) / (hi - lo)
    return 1.96 if dof > 29 else tab[1]


def mean_ci(xs):
    xs = np.asarray(xs, float)
    if xs.size < 2:
        return float(xs.mean()), 0.0
    return float(xs.mean()), float(t95(xs.size - 1) * xs.std(ddof=1) / np.sqrt(xs.size))
