"""ApproxKit v2 tests (plain asserts; pytest-compatible; also `python tests/test_v2.py`).

Independent SCALAR (python-int, bit-by-bit) references are compared with the
vectorised models over ALL operand combinations at n = 5 and 6, plus structural
properties, metric hand-checks and the ETA-vs-truncation accumulator test.
RTL <-> model equivalence is a separate flow: scripts/rtl_regress_v2.py (iverilog).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import numpy as np
from approxkit_v2 import signed as S, divsq as D
from approxkit_v2.common import ak
from approxkit_v2.metrics2 import evaluate

NS = (5, 6)


def allpairs(u):
    r = np.arange(u.lo, u.hi + 1)
    return np.repeat(r, r.size), np.tile(r, r.size)


def w2(v, bits):
    v &= (1 << bits) - 1
    return v - (1 << bits) if v >> (bits - 1) else v


# ---------------- scalar references ----------------
def ref_bw(a, b, n, k, comp_val=0):
    A, B = a & ((1 << n) - 1), b & ((1 << n) - 1)
    tot = (1 << n) + (1 << (2 * n - 1))
    for i in range(n):
        for j in range(n):
            if i + j < k:
                continue
            t = ((A >> i) & 1) & ((B >> j) & 1)
            if (i == n - 1) != (j == n - 1):
                t ^= 1
            tot += t << (i + j)
    return w2(tot + comp_val, 2 * n)


def ref_bw_orlow(a, b, n, k):
    A, B = a & ((1 << n) - 1), b & ((1 << n) - 1)
    hi, lo = (1 << n) + (1 << (2 * n - 1)), 0
    for i in range(n):
        for j in range(n):
            t = ((A >> i) & 1) & ((B >> j) & 1)
            if (i == n - 1) != (j == n - 1):
                t ^= 1
            if i + j < k:
                lo |= t << (i + j)
            else:
                hi += t << (i + j)
    return w2(hi | lo, 2 * n)


def ref_sdrum(a, b, n, k):
    def seg(x):
        y = ~x if x < 0 else x
        p = y.bit_length() - 1          # -1 for 0/-1
        s = max(p + 2 - k, 0)
        g = x >> s
        return (g | 1 if s > 0 else g), s
    (sa, ha), (sb, hb) = seg(a), seg(b)
    return w2(sa * sb * (1 << (ha + hb)), 2 * n)


def ref_restoring(a, b, n, t, comp):
    """literal restoring division, last t iterations removed"""
    if b == 0:
        return (1 << n) - 1
    q, r = 0, 0
    for i in range(n - 1, t - 1, -1):
        r = (r << 1) | ((a >> i) & 1)
        if r >= b:
            r -= b
            q |= 1 << i
    if comp and t:
        q |= 1 << (t - 1)
    return q


def ref_sq_trunc(a, n, k, comp):
    tot = 0
    for i in range(n):
        for j in range(n):
            if ((a >> i) & 1) and ((a >> j) & 1):
                col = i + j if i == j else None
        # explicit full matrix a_i a_j 2^(i+j), folded: i==j once, i<j doubled
    for i in range(n):
        if (a >> i) & 1 and 2 * i >= k:
            tot += 1 << (2 * i)
        for j in range(i + 1, n):
            if (a >> i) & 1 and (a >> j) & 1 and i + j + 1 >= k:
                tot += 1 << (i + j + 1)
    return tot


# ---------------- tests ----------------
def test_bw_exact_is_exact_and_never_wraps():
    for n in NS + (8,):
        u = S.BaughWooley(n, 0)
        a, b = allpairs(u)
        assert np.array_equal(u(a, b), a * b)


def test_bw_trunc_and_orlow_vs_scalar():
    for n in NS:
        a, b = allpairs(S.BaughWooley(n, 0))
        for k in range(0, 2 * n + 1, 1 if n == 5 else 2):
            u = S.BaughWooley(n, k, "trunc", 0)
            exp = np.array([ref_bw(int(x), int(y), n, k) for x, y in zip(a, b)])
            assert np.array_equal(u(a, b), exp), (n, k)
            u = S.BaughWooley(n, k, "trunc", 1)
            cv = (S.bw_expected_dropped_q4(n, k) + 2) // 4
            exp = np.array([ref_bw(int(x), int(y), n, k, cv) for x, y in zip(a, b)])
            assert np.array_equal(u(a, b), exp), (n, k, "comp")
        for k in range(0, n + 1):
            u = S.BaughWooley(n, k, "orlow")
            exp = np.array([ref_bw_orlow(int(x), int(y), n, k) for x, y in zip(a, b)])
            assert np.array_equal(u(a, b), exp), (n, k, "orlow")


def test_bw_trunc_never_overestimates_without_comp():
    n = 6
    u = S.BaughWooley(n, 5)
    a, b = allpairs(u)
    # dropped bits are non-negative contributions of NAND'ed/AND'ed pp bits -> error <= 0
    assert (u(a, b) - a * b).max() <= 0


def test_sdrum_vs_scalar_and_exact_when_window_covers():
    for n in NS:
        for k in range(2, n + 1):
            u = S.SignedDRUM(n, k)
            a, b = allpairs(u)
            exp = np.array([ref_sdrum(int(x), int(y), n, k) for x, y in zip(a, b)])
            assert np.array_equal(u(a, b), exp), (n, k)
        u = S.SignedDRUM(n, n)
        a, b = allpairs(u)
        assert np.array_equal(u(a, b), a * b)      # k = n -> exact


def test_signmag_exact_and_sat():
    n = 6
    u = S.SignMag(lambda w: ak.ExactMultiplier(w), n, 0)
    a, b = allpairs(u)
    assert np.array_equal(u(a, b), a * b)
    u = S.SignMag(lambda w: ak.ExactMultiplier(w), n, 1)
    a, b = allpairs(u)
    lo = -(1 << (n - 1))
    ok = (a != lo) & (b != lo)
    assert np.array_equal(u(a, b)[ok], (a * b)[ok])
    assert (u(a, b)[~ok] != (a * b)[~ok]).any()     # the single saturated value is wrong, by design


def test_restoring_divider_vs_literal_loop():
    for n in NS:
        for t in range(0, n):
            for c in (0, 1):
                u = D.TruncRestoringDiv(n, t, c)
                a, b = allpairs(u)
                exp = np.array([ref_restoring(int(x), int(y), n, t, c) for x, y in zip(a, b)])
                assert np.array_equal(u(a, b), exp), (n, t, c)
    u = D.TruncRestoringDiv(6, 0)
    a, b = allpairs(u)
    assert np.array_equal(u(a, b), u.exact(a, b))


def test_divider_conventions_and_bounds():
    for u in (D.MitchellDiv(6), D.DRUMDiv(6, 3), D.TruncRestoringDiv(6, 2, 1), D.ExactDivider(6)):
        a, b = allpairs(u)
        q = u(a, b)
        assert (q >= 0).all() and (q <= 63).all()
        assert (q[b == 0] == 63).all()                       # b = 0 -> all ones, every variant
    # DRUM divider with a full-width window is exact on the whole space
    u = D.DRUMDiv(6, 6)
    a, b = allpairs(u)
    assert np.array_equal(u(a, b)[b > 0], u.exact(a, b)[b > 0])


def test_mitchell_exact_on_powers_of_two():
    n = 8
    p2 = 1 << np.arange(n)
    a, b = np.repeat(p2, n), np.tile(p2, n)
    assert np.array_equal(D.MitchellDiv(n)(a, b), a // b)
    assert np.array_equal(D.MitchellRecip(n, 2 * n)(p2), (1 << (2 * n)) // p2)
    assert np.array_equal(D.MitchellSquarer(n)(p2), p2 * p2)


def test_reciprocal_exact_and_error_sign():
    n = 6
    r = D.ExactReciprocal(n, 2 * n)
    x = np.arange(1, 1 << n)
    assert np.array_equal(r(x), (1 << (2 * n)) // x)
    assert r(np.array([0]))[0] == (1 << (2 * n + 1)) - 1
    m = D.MitchellRecip(n, 2 * n)
    assert np.array_equal(m(np.array([0])), r(np.array([0])))


def test_squarer_trunc_vs_scalar_and_drum_full_window():
    for n in NS:
        a = np.arange(1 << n)
        for k in range(0, 2 * n + 1):
            u = D.TruncSquarer(n, k, 0)
            exp = np.array([ref_sq_trunc(int(x), n, k, 0) for x in a])
            assert np.array_equal(u(a), exp), (n, k)
            assert (u(a) <= a * a).all()
        assert np.array_equal(D.TruncSquarer(n, 0)(a), a * a)
        assert np.array_equal(D.DRUMSquarer(n, n)(a), a * a)


def test_metrics_hand_check():
    u = S.BaughWooley(5, 4)
    a, b = allpairs(u)
    e = u(a, b) - a * b
    m = evaluate(u)
    assert m["mode"] == "exhaustive" and m["MED_ci"] == 0
    assert abs(m["MED"] - np.abs(e).mean()) < 1e-12 and m["WCE"] == np.abs(e).max()
    assert abs(m["bias"] - e.mean()) < 1e-12 and abs(m["ER"] - (e != 0).mean()) < 1e-12
    z = evaluate(S.BaughWooley(5, 0))
    assert z["MED"] == 0 and z["WCE"] == 0 and z["ER"] == 0
    # divider: b = 0 must be excluded from the metric population
    d = D.TruncRestoringDiv(5, 2)
    aa, bb = allpairs(d)
    v = bb != 0
    ee = np.abs(d(aa, bb) - d.exact(aa, bb))[v]
    assert evaluate(d)["samples"] == v.sum() and abs(evaluate(d)["MED"] - ee.mean()) < 1e-12


def test_monte_carlo_matches_exhaustive_within_ci():
    u = S.BaughWooley(8, 8, "trunc", 0)       # 65536 pairs: force MC for a check against exact value
    ex = evaluate(u, mode="exhaustive")
    mc = evaluate(u, samples=1 << 18, seed=3, mode="monte_carlo")
    assert abs(mc["MED"] - ex["MED"]) <= 3 * mc["MED_ci"] + 1e-9


def test_eta_equals_truncation_accumulator_mechanism():
    """audit observation: with positive/negative accumulators of the v1 apps harness,
    low-k bits of ETA saturate to 2^k-1 in BOTH accumulators and cancel in pos-neg."""
    from eta_vs_trunc import run_mechanism
    r = run_mechanism(k=4, nterm=256, ntrial=200, seed=1)
    assert r["frac_both_low_saturated"] > 0.99          # low 4 bits == 15 in both accumulators
    assert r["frac_identical_diff"] > 0.99              # ETA and truncation give identical dot results
    assert r["frac_identical_diff_short"] < 0.9         # but not for very short sums (T=2): mechanism is saturation
    assert r["single_acc_differs"] > 0.3                # a single two's-complement accumulator breaks the equality


def test_quality_metrics():
    from approxkit_v2.quality import psnr, ssim, snr_db
    r = np.random.default_rng(0)
    x = r.integers(0, 256, (40, 40)).astype(float)
    y = np.clip(x + r.normal(0, 10, x.shape), 0, 255)
    assert ssim(x, x) > 1 - 1e-12 and psnr(x, x) == 99.0
    mse = np.mean((x - y) ** 2)
    assert abs(psnr(x, y) - 10 * np.log10(255 ** 2 / mse)) < 1e-9
    assert 0 < ssim(x, y) < 1 and ssim(x, y) > ssim(x, np.clip(x + r.normal(0, 40, x.shape), 0, 255))
    assert abs(snr_db(x, y) - 10 * np.log10(np.mean(x ** 2) / mse)) < 1e-9
    try:                       # cross-check against scikit-image when installed
        from skimage.metrics import structural_similarity as sk
        assert abs(ssim(x, y) - sk(x, y, data_range=255, gaussian_weights=True, sigma=1.5, use_sample_covariance=False)) < 2e-3
        print("   (ssim cross-checked with scikit-image)")
    except ImportError:
        print("   (scikit-image not installed: ssim cross-check skipped)")


def test_conv_uses_unit_and_exact_matches_reference():
    from app_quality import conv, KERNELS, images
    ex = S.BaughWooley(8, 0)
    n, im = images()[0]
    k, sh, mode = KERNELS["blur_odd"]
    o = conv(ex, im, k, sh, mode)
    s = im - 128
    acc = sum(int(k[i, j]) * s[i:i + im.shape[0] - 2, j:j + im.shape[1] - 2] for i in range(3) for j in range(3))
    assert np.array_equal(o, np.clip((acc >> sh) + 128, 0, 255))
    assert not np.array_equal(conv(S.BaughWooley(8, 10), im, k, sh, mode), o)     # a bad unit changes the picture


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import time
    fails = 0
    for name, f in sorted(globals().items()):
        if name.startswith("test_"):
            t0 = time.time()
            try:
                f(); print(f"PASS {name} ({time.time()-t0:.1f}s)")
            except Exception as e:
                import traceback; traceback.print_exc(); print(f"FAIL {name}: {e!r}"); fails += 1
    print("ALL PASS" if not fails else f"{fails} FAILED"); sys.exit(1 if fails else 0)
