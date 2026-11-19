"""Single source of truth: every v2 design = Python model + Verilog module + parameters.

Used by the RTL regression, the synthesis cost scripts, the error sweep and the
application harness, so that a design name means the same thing everywhere.
"""
import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List
from .common import V1_ROOT, V2_ROOT, ak
from . import signed as S
from . import divsq as D

V1_RTL = os.path.join(V1_ROOT, "rtl")
V2_RTL = os.path.join(V2_ROOT, "rtl")


@dataclass
class Entry:
    key: str
    group: str                     # smul | div | recip | sq
    role: str                      # exact | approx
    make: Callable                 # n -> python Unit
    module: str
    params: Dict[str, int]         # Verilog parameters (N included)
    srcs: List[str] = field(default_factory=list)   # absolute paths of .v files
    family: str = ""

    def unit(self):
        return self.make(self.params["N"])


def _v2(*names):
    return [os.path.join(V2_RTL, x + ".v") for x in names]


def _v1(*names):
    return [os.path.join(V1_RTL, x + ".v") for x in names]


def _ks(n, fr):
    return sorted({max(1, int(round(n * f))) for f in fr})


def catalog(n=8):
    E = []

    def add(*a, **k):
        E.append(Entry(*a, **k))

    # ---------------------------------------------------------- signed multipliers
    add("s_exact", "smul", "exact", lambda n: S.BaughWooley(n, 0), "s_exact_mult", {"N": n}, _v2("s_exact_mult"), "s_exact")
    add("bw_k0", "smul", "exact", lambda n: S.BaughWooley(n, 0), "bw_mult", {"N": n, "K": 0, "MODE": 0, "COMP": 0}, _v2("bw_mult"), "bw")
    for k in _ks(n, (.25, .5, .75, 1.0, 1.5)):
        for c in (0, 1):
            add(f"bw_trunc_k{k}_c{c}", "smul", "approx", (lambda n, k=k, c=c: S.BaughWooley(n, k, "trunc", c)),
                "bw_mult", {"N": n, "K": k, "MODE": 0, "COMP": c}, _v2("bw_mult"), "bw")
    for k in _ks(n, (.25, .5, .75, 1.0)):
        add(f"bw_orlow_k{k}", "smul", "approx", (lambda n, k=k: S.BaughWooley(n, k, "orlow")),
            "bw_mult", {"N": n, "K": k, "MODE": 1, "COMP": 0}, _v2("bw_mult"), "bw")
    for k in sorted({kk for kk in (2, 3, 4, 5, 6, 7) if 2 <= kk < n}):
        add(f"sdrum_k{k}", "smul", "approx", (lambda n, k=k: S.SignedDRUM(n, k)),
            "sdrum_mult", {"N": n, "K": k}, _v2("sdrum_mult"), "sdrum")
    # sign-magnitude wrappers around v1 unsigned units
    sm = [("exact", 0, 0, 0, lambda w: ak.ExactMultiplier(w), []),
          ("trunc", 1, max(1, n // 2), 0, None, ["trunc_mult", "bam_mult"]),
          ("trunc", 1, n, 0, None, ["trunc_mult", "bam_mult"]),
          ("drum", 2, 3, 0, None, ["drum_mult"]), ("drum", 2, 4, 0, None, ["drum_mult"]),
          ("drum", 2, 5, 0, None, ["drum_mult"]), ("drum", 2, 6, 0, None, ["drum_mult"]),
          ("mitchell", 3, 0, 0, lambda w: ak.Mitchell(w), ["mitchell_mult"]),
          ("booth", 4, max(1, n // 2), 0, None, ["booth_mult"]), ("booth", 4, n, 0, None, ["booth_mult"]),
          ("bam", 5, n - 2, 2, None, ["bam_mult"])]
    for fam, uid, p1, p2, mk, v1f in sm:
        if fam == "drum" and not (2 <= p1 <= n):
            continue
        if mk is None:
            mk = {"trunc": lambda w, p1=p1: ak.TruncatedMultiplier(w, p1),
                  "drum": lambda w, p1=p1: ak.DRUM(w, p1),
                  "booth": lambda w, p1=p1: ak.BoothApprox(w, p1),
                  "bam": lambda w, p1=p1, p2=p2: ak.BrokenArray(w, p1, p2)}[fam]
        for sat in ((0, 1) if (fam in ("exact", "mitchell") or (fam, p1) in (("drum", 4), ("trunc", n))) else (0,)):
            lab = fam + (f"{p1}" if fam in ("trunc", "drum", "booth") else "") + (f"_{p1}_{p2}" if fam == "bam" else "")
            add(f"smag_{lab}" + ("_sat" if sat else ""), "smul", "exact" if fam == "exact" else "approx",
                (lambda n, mk=mk, sat=sat: S.SignMag(mk, n, sat)), "smag_mult",
                {"N": n, "UNIT": uid, "P1": p1, "P2": p2, "SAT": sat}, _v2("smag_mult") + _v1(*v1f), "smag")
    # ---------------------------------------------------------- dividers
    add("div_exact", "div", "exact", lambda n: D.ExactDivider(n), "div_exact", {"N": n}, _v2("div_exact"), "div_exact")
    add("div_trunc_t0", "div", "exact", lambda n: D.TruncRestoringDiv(n, 0), "div_trunc", {"N": n, "T": 0, "COMP": 0}, _v2("div_trunc"), "div_trunc")
    for t in sorted({t for t in (1, 2, 3, 4, 5, 6) if t <= n - 1}):
        for c in (0, 1):
            add(f"div_trunc_t{t}_c{c}", "div", "approx", (lambda n, t=t, c=c: D.TruncRestoringDiv(n, t, c)),
                "div_trunc", {"N": n, "T": t, "COMP": c}, _v2("div_trunc"), "div_trunc")
    add("div_mitchell", "div", "approx", lambda n: D.MitchellDiv(n), "div_mitchell", {"N": n}, _v2("div_mitchell"), "div_mitchell")
    for k in sorted({k for k in (2, 3, 4, 5, 6) if k <= n}):
        add(f"div_drum_k{k}", "div", "approx", (lambda n, k=k: D.DRUMDiv(n, k)), "div_drum", {"N": n, "K": k}, _v2("div_drum"), "div_drum")
    # ---------------------------------------------------------- reciprocal (m = 2n)
    add("recip_exact", "recip", "exact", lambda n: D.ExactReciprocal(n, 2 * n), "recip_exact", {"N": n, "M": 2 * n}, _v2("recip_exact"), "recip_exact")
    add("recip_mitchell", "recip", "approx", lambda n: D.MitchellRecip(n, 2 * n), "recip_mitchell", {"N": n, "M": 2 * n}, _v2("recip_mitchell"), "recip_mitchell")
    # ---------------------------------------------------------- squarers
    add("sq_exact", "sq", "exact", lambda n: D.ExactSquarer(n), "sq_trunc", {"N": n, "K": 0, "COMP": 0}, _v2("sq_trunc"), "sq_exact")
    for k in _ks(n, (.25, .5, .75, 1.0, 1.25)):
        for c in (0, 1):
            add(f"sq_trunc_k{k}_c{c}", "sq", "approx", (lambda n, k=k, c=c: D.TruncSquarer(n, k, c)),
                "sq_trunc", {"N": n, "K": k, "COMP": c}, _v2("sq_trunc"), "sq_trunc")
    for k in sorted({k for k in (3, 4, 5, 6) if k <= n}):
        add(f"sq_drum_k{k}", "sq", "approx", (lambda n, k=k: D.DRUMSquarer(n, k)), "sq_drum", {"N": n, "K": k}, _v2("sq_drum"), "sq_drum")
    add("sq_mitchell", "sq", "approx", lambda n: D.MitchellSquarer(n), "sq_mitchell", {"N": n}, _v2("sq_mitchell"), "sq_mitchell")
    # unique keys
    seen = set()
    out = []
    for e in E:
        if e.key not in seen:
            seen.add(e.key)
            out.append(e)
    return out
