"""Error metrics: exhaustive (n <= 12) or Monte Carlo (wider).

Error e = approx - exact (signed).  With N evaluated input pairs:
  ER   = P(e != 0)
  MED  = mean |e|
  MRED = mean |e| / exact  over pairs with exact != 0
  NMED = MED / (maximum exact output)
  WCE  = max |e|          (observed; a lower bound on the true worst case
                           when Monte Carlo is used)
  bias = mean e           (signed)
The full signed error PMF is kept so histograms of any binning can be drawn.
"""
from dataclasses import dataclass, field
import numpy as np

EXHAUSTIVE_MAX_BITS = 12


@dataclass
class ErrorMetrics:
    unit: str
    n: int
    mode: str
    samples: int
    ER: float
    MED: float
    MRED: float
    NMED: float
    WCE: int
    bias: float
    rmse: float
    err_values: np.ndarray = field(repr=False, default=None)
    err_counts: np.ndarray = field(repr=False, default=None)

    def row(self):
        return {"unit": self.unit, "n": self.n, "mode": self.mode, "samples": self.samples,
                "ER": self.ER, "MED": self.MED, "MRED": self.MRED, "NMED": self.NMED,
                "WCE": self.WCE, "bias": self.bias, "RMSE": self.rmse}

    def histogram(self, bins=61, normalize=False):
        """Histogram of the signed error. Returns (counts, edges).
        Range covers [-WCE, WCE]; bins have an integer width (>=1) so integer
        errors never alias against the bin grid (at most `bins` bins)."""
        w = max(int(self.WCE), 1)
        width = -(-(2 * w + 1) // bins)                      # ceil
        nb = -(-(2 * w + 1) // width)
        edges = -w - 0.5 + width * np.arange(nb + 1)
        counts, _ = np.histogram(self.err_values, bins=edges, weights=self.err_counts)
        if normalize:
            counts = counts / counts.sum()
        return counts, edges


def _accumulate(unit, a, b, state):
    ex = unit.exact(a, b)
    e = unit(a, b) - ex
    state["n"] += e.size
    state["nz"] += int(np.count_nonzero(e))
    ae = np.abs(e)
    state["sae"] += float(ae.sum())
    state["se"] += float(e.sum())
    state["se2"] += float((e.astype(np.float64) ** 2).sum())
    state["wce"] = max(state["wce"], int(ae.max()))
    m = ex != 0
    state["srel"] += float((ae[m] / ex[m]).sum())
    state["nrel"] += int(m.sum())
    v, c = np.unique(e, return_counts=True)
    for vi, ci in zip(v.tolist(), c.tolist()):
        state["pmf"][vi] = state["pmf"].get(vi, 0) + ci


def evaluate(unit, mode="auto", samples=1 << 20, seed=0, chunk=1 << 22):
    """Evaluate error metrics of `unit`.

    mode: 'exhaustive' | 'monte_carlo' | 'auto' (exhaustive iff n <= 12).
    Monte Carlo draws `samples` uniform operand pairs from a seeded RNG.
    """
    n = unit.n
    if mode == "auto":
        mode = "exhaustive" if n <= EXHAUSTIVE_MAX_BITS else "monte_carlo"
    st = dict(n=0, nz=0, sae=0.0, se=0.0, se2=0.0, wce=0, srel=0.0, nrel=0, pmf={})
    if mode == "exhaustive":
        if n > 14:
            raise ValueError("exhaustive evaluation limited to n <= 14")
        full = np.arange(1 << n, dtype=np.int64)
        rows = max(1, chunk >> n)
        for a0 in range(0, 1 << n, rows):
            av = full[a0:a0 + rows]
            A = np.repeat(av, 1 << n)
            B = np.tile(full, av.size)
            _accumulate(unit, A, B, st)
    elif mode == "monte_carlo":
        rng = np.random.default_rng(seed)
        left = samples
        while left > 0:
            m = min(left, chunk)
            A = rng.integers(0, 1 << n, m, dtype=np.int64)
            B = rng.integers(0, 1 << n, m, dtype=np.int64)
            _accumulate(unit, A, B, st)
            left -= m
    else:
        raise ValueError(mode)
    N = st["n"]
    vals = np.array(sorted(st["pmf"]), dtype=np.int64)
    cnts = np.array([st["pmf"][v] for v in vals.tolist()], dtype=np.int64)
    return ErrorMetrics(
        unit=unit.name, n=n, mode=mode, samples=N,
        ER=st["nz"] / N, MED=st["sae"] / N,
        MRED=st["srel"] / max(st["nrel"], 1),
        NMED=st["sae"] / N / unit.max_exact,
        WCE=st["wce"], bias=st["se"] / N, rmse=(st["se2"] / N) ** 0.5,
        err_values=vals, err_counts=cnts)
