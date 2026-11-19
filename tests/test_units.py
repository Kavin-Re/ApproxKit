"""Plain-assert tests (pytest-compatible; also runnable: python tests/test_units.py).

Independent *scalar, bit-by-bit* reference models are compared with the
vectorised numpy units over all operand pairs at n=6 (4096 pairs).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import approxkit as ak
from approxkit import apps, dse
from approxkit.metrics import evaluate

N = 6


def bits(x, n):
    return [(x >> i) & 1 for i in range(n)]


def val(bl):
    return sum(b << i for i, b in enumerate(bl))


# ---------------- scalar references ----------------
def ref_loa(n, k):
    def f(a, b):
        A, B = bits(a, n), bits(b, n)
        low = [A[i] | B[i] for i in range(k)]
        cin = A[k - 1] & B[k - 1] if k else 0
        hi = val(A[k:]) + val(B[k:]) + cin
        return val(low) + (hi << k)
    return f


def ref_trunc(n, k):
    return lambda a, b: ((a >> k) + (b >> k)) << k


def ref_eta(n, k):
    def f(a, b):
        A, B = bits(a, n), bits(b, n)
        low, flag = [0] * k, 0
        for i in range(k - 1, -1, -1):
            if flag:
                low[i] = 1
            elif A[i] and B[i]:
                low[i] = 1
                flag = 1
            else:
                low[i] = A[i] ^ B[i]
        return val(low) + (((a >> k) + (b >> k)) << k)
    return f


def ref_ama(n, k, cell):
    def f(a, b):
        A, B = bits(a, n), bits(b, n)
        c, low = 0, []
        for i in range(k):
            if cell == "sum_b":
                s, c = B[i], A[i]
            else:
                c = 1 if A[i] + B[i] + c >= 2 else 0
                s = 1 - c
            low.append(s)
        return val(low) + (((a >> k) + (b >> k) + c) << k)
    return f


def ref_gear(n, R, P):
    L = R + P

    def f(a, b):
        # one explicit sub-adder per output block, window = bits [s, s+L)
        out = [0] * (n + 1)
        s_prev_top, j = 0, 0
        while True:
            s = min(j * R, n - L)
            t = ((a >> s) & ((1 << L) - 1)) + ((b >> s) & ((1 << L) - 1))
            tb = bits(t, L + 1)
            last = s + L >= n
            for pos in range(s_prev_top, s + L):
                out[pos] = tb[pos - s]
            if last:
                out[n] = tb[L]
                break
            s_prev_top = s + L
            j += 1
        return val(out)
    return f


def ref_bam(n, vbl, hbl):
    def f(a, b):
        tot = 0
        for i in range(hbl, n):
            for j in range(n):
                if i + j >= vbl:
                    tot += ((a >> j) & 1) * ((b >> i) & 1) << (i + j)
        return tot
    return f


def ref_drum(n, k):
    def seg(x):
        if x == 0:
            return 0, 0
        l = x.bit_length() - 1
        if l < k:
            return x, 0
        sh = l - k + 1
        return (x >> sh) | 1, sh
    def f(a, b):
        sa, ha = seg(a)
        sb, hb = seg(b)
        return (sa * sb) << (ha + hb)
    return f


def ref_mitchell(n):
    from fractions import Fraction
    F = n - 1
    def f(a, b):
        if a == 0 or b == 0:
            return 0
        la, lb = a.bit_length() - 1, b.bit_length() - 1
        xa = Fraction(a, 1 << la) - 1
        xb = Fraction(b, 1 << lb) - 1
        s = xa + xb
        if s < 1:
            r = Fraction(1 << (la + lb)) * (1 + s)
        else:
            r = Fraction(1 << (la + lb + 1)) * s
        return int(r)   # floor
    return f


def ref_booth(n, k):
    D = n // 2 + 1
    def f(a, b):
        tot = 0
        for i in range(D):
            bm1 = (b >> (2 * i - 1)) & 1 if i else 0
            d = -2 * ((b >> (2 * i + 1)) & 1) + ((b >> (2 * i)) & 1) + bm1
            mag = (abs(d) * a) << (2 * i)
            mag = (mag >> k) << k
            tot += mag if d >= 0 else -mag
        return tot % (1 << (2 * n))
    return f


def check(unit, ref, n=N):
    a = np.repeat(np.arange(1 << n), 1 << n)
    b = np.tile(np.arange(1 << n), 1 << n)
    got = unit(a, b)
    exp = np.array([ref(int(x), int(y)) for x, y in zip(a, b)])
    bad = np.flatnonzero(got != exp)
    assert bad.size == 0, f"{unit.name}: {bad.size} mismatches, e.g. a={a[bad[0]]} b={b[bad[0]]} got={got[bad[0]]} ref={exp[bad[0]]}"


def test_adders_vs_scalar_reference():
    for k in range(0, N):
        check(ak.LOA(N, k), ref_loa(N, k))
        check(ak.TruncatedAdder(N, k), ref_trunc(N, k))
        check(ak.ETA(N, k), ref_eta(N, k))
        for cell in ("sum_b", "inv_carry"):
            check(ak.AMA(N, k, cell), ref_ama(N, k, cell))
    for R in range(1, N + 1):
        for P in range(0, N - R + 1):
            check(ak.GeAr(N, R, P), ref_gear(N, R, P))


def test_multipliers_vs_scalar_reference():
    for vbl in range(0, 2 * N + 1, 2):
        for hbl in range(0, N + 1, 2):
            check(ak.BrokenArray(N, vbl, hbl), ref_bam(N, vbl, hbl))
    for k in range(2, N + 1):
        check(ak.DRUM(N, k), ref_drum(N, k))
    check(ak.Mitchell(N), ref_mitchell(N))
    check(ak.Mitchell(5), ref_mitchell(5), n=5)
    for k in range(0, 2 * N + 1):
        check(ak.BoothApprox(N, k), ref_booth(N, k))
    for n in (5, 7):   # odd widths
        check(ak.BoothApprox(n, 3), ref_booth(n, 3), n=n)


def test_exact_limits():
    for n in (5, 6):
        a = np.repeat(np.arange(1 << n), 1 << n)
        b = np.tile(np.arange(1 << n), 1 << n)
        for u in [ak.LOA(n, 0), ak.TruncatedAdder(n, 0), ak.ETA(n, 0), ak.AMA(n, 0), ak.GeAr(n, 1, n - 1),
                  ak.GeAr(n, 2, n - 2), ak.GeAr(n, n, 0)]:
            assert np.array_equal(u(a, b), a + b), u.name
        for u in [ak.BrokenArray(n, 0, 0), ak.TruncatedMultiplier(n, 0), ak.BoothApprox(n, 0),
                  ak.DRUM(n, n), ak.ExactMultiplier(n)]:
            assert np.array_equal(u(a, b), a * b), u.name


def test_known_values_and_properties():
    # LOA n=4,k=2: 7+3 -> low=0b11, cin=1, hi=1+0+1=2 -> 0b1011 = 11 (exact 10)
    assert int(ak.LOA(4, 2)(7, 3)) == 11
    # Mitchell: 3*3 -> 8 (classic 11.1% underestimate), 2^i*2^j exact
    assert int(ak.Mitchell(8)(3, 3)) == 8
    assert int(ak.Mitchell(8)(64, 32)) == 2048
    # DRUM keeps small operands exact
    assert int(ak.DRUM(8, 4)(13, 9)) == 117
    n = 8
    a = np.repeat(np.arange(1 << n), 1 << n)
    b = np.tile(np.arange(1 << n), 1 << n)
    ex = a * b
    mi = ak.Mitchell(n)(a, b)
    assert np.all(mi <= ex)                               # never overestimates
    nz = ex > 0
    assert np.max((ex[nz] - mi[nz]) / ex[nz]) < 0.126     # ~11.1% analytic + floor effect
    for k in (2, 5, 8):
        assert np.all(ak.TruncatedMultiplier(n, k)(a, b) <= ex)
        assert np.all(np.abs(ak.TruncatedAdder(n, k % n)(a, b) - (a + b)) < 2 ** (k % n + 1))


def test_range_check():
    try:
        ak.LOA(4, 2)(16, 0)
    except ValueError:
        return
    raise AssertionError("range check missing")


def test_metrics_vs_bruteforce():
    u = ak.LOA(5, 2)
    m = evaluate(u, mode="exhaustive")
    errs, rel = [], []
    for a in range(32):
        for b in range(32):
            e = int(u(a, b)) - (a + b)
            errs.append(e)
            if a + b:
                rel.append(abs(e) / (a + b))
    errs = np.array(errs)
    assert m.samples == 1024
    assert abs(m.ER - np.mean(errs != 0)) < 1e-12
    assert abs(m.MED - np.mean(np.abs(errs))) < 1e-12
    assert abs(m.MRED - np.mean(rel)) < 1e-12
    assert abs(m.NMED - np.mean(np.abs(errs)) / 62) < 1e-12
    assert m.WCE == np.abs(errs).max()
    assert abs(m.bias - errs.mean()) < 1e-12
    c, _ = m.histogram(bins=11)
    assert c.sum() == 1024
    # exact unit has zero error everywhere
    z = evaluate(ak.ExactAdder(6))
    assert z.ER == 0 and z.WCE == 0 and z.MRED == 0


def test_monte_carlo_close_to_exhaustive():
    u = ak.DRUM(8, 4)
    ex = evaluate(u, mode="exhaustive")
    mc = evaluate(u, mode="monte_carlo", samples=1 << 19, seed=1)
    assert mc.mode == "monte_carlo" and mc.samples == 1 << 19
    assert abs(mc.MRED - ex.MRED) < 0.002 and abs(mc.NMED - ex.NMED) < 0.0005
    assert mc.WCE <= ex.WCE
    # reproducible
    assert evaluate(u, mode="monte_carlo", samples=1 << 12, seed=3).MED == \
        evaluate(u, mode="monte_carlo", samples=1 << 12, seed=3).MED


def test_wide_units_run():
    for u in [ak.LOA(24, 8), ak.GeAr(24, 4, 8), ak.ETA(24, 10), ak.AMA(24, 8), ak.DRUM(16, 6), ak.Mitchell(16),
              ak.BoothApprox(16, 10), ak.BrokenArray(16, 12, 3)]:
        m = evaluate(u, samples=1 << 14)
        assert m.mode == "monte_carlo" and m.WCE >= 0


def test_cost_monotonic_and_pareto():
    for n in (8,):
        assert ak.LOA(n, 6).area_ge() < ak.LOA(n, 2).area_ge() < ak.ExactAdder(n).area_ge()
        assert ak.TruncatedMultiplier(n, 10).area_ge() < ak.TruncatedMultiplier(n, 4).area_ge() < ak.ExactMultiplier(n).area_ge()
    rows = [dict(family="x", area_ge=c, MRED=e) for c, e in [(1, 5), (2, 3), (2, 4), (3, 3), (4, 1), (5, 2)]]
    f = dse.pareto_front(rows)
    assert [(r["area_ge"], r["MRED"]) for r in f] == [(1, 5), (2, 3), (4, 1)]
    ds = dse.default_designs("adder", 8) + dse.default_designs("multiplier", 8)
    assert len(ds) > 60 and len({d.name for d in ds}) == len(ds)


def test_apps():
    rng = np.random.default_rng(0)
    x = rng.integers(-127, 128, (4, 32))
    y = rng.integers(-127, 128, (4, 32))
    ex = ak.ExactMultiplier(8)
    assert np.array_equal(apps.dot(x, y, ex), (x * y).sum(-1))
    assert np.array_equal(apps.dot(x, y, ex, adder=ak.ExactAdder(20)), (x * y).sum(-1))
    h = rng.integers(-100, 101, 7)
    sig = rng.integers(-100, 101, 50)
    assert np.array_equal(apps.fir(sig, h, ex), np.convolve(sig, h, mode="valid"))
    A = rng.integers(0, 256, (5, 6)); B = rng.integers(0, 256, (6, 3))
    assert np.array_equal(apps.matmul(A, B, ex), A @ B)
    r = apps.accuracy(apps.matmul(A, B, ak.DRUM(8, 4)), A @ B)
    assert 0.0 < r["rel_rms"] < 0.1


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for f in fns:
        f()
        print("PASS", f.__name__)
    print(f"{len(fns)} tests passed")
