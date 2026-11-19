#!/usr/bin/env python3
"""v2 RTL regression: Python model -> exhaustive vectors -> self-checking iverilog bench.

Every operand combination (2^(2n) for two-input units, 2^n for one-input units,
INCLUDING b = 0 for dividers and x = 0 for the reciprocal) is applied and the
raw output bits (two's complement for signed units) are compared with the model.

  python scripts/rtl_regress_v2.py --out build --widths 8,6,5
  python scripts/rtl_regress_v2.py --out build --mutations      # bench must FAIL on broken RTL
"""
import argparse, os, re, subprocess, sys, shutil, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from approxkit_v2.catalog import catalog, V2_RTL
from approxkit_v2.common import mask_bits

TB = """`timescale 1ns/1ps
module tb;
  localparam N = {n}; localparam OW = {ow}; localparam NV = {nv};
  reg [{wordw}-1:0] mem [0:NV-1];
  reg  [N-1:0] a, b;
  wire [OW-1:0] y;
  reg  [OW-1:0] expv;
  integer i, errs;
  {mod} #({plist}) dut ({ports});
  initial begin
    $readmemh("{vec}", mem);
    errs = 0;
    for (i = 0; i < NV; i = i + 1) begin
      a = mem[i][OW{bsh} +: N];
      b = {bexpr};
      expv = mem[i][OW-1:0];
      #1;
      if (y !== expv) begin
        errs = errs + 1;
        if (errs <= 5) $display("MISMATCH a=%0d b=%0d got=%0d expected=%0d", a, b, y, expv);
      end
    end
    if (errs == 0) $display("RESULT PASS %0d vectors", NV);
    else           $display("RESULT FAIL %0d of %0d vectors", errs, NV);
    $finish;
  end
endmodule
"""


def vectors(unit):
    n = unit.n
    r = np.arange(unit.lo, unit.hi + 1, dtype=np.int64)
    if unit.inputs == 2:
        a, b = np.repeat(r, r.size), np.tile(r, r.size)
    else:
        a, b = r, None
    y = mask_bits(unit(a, b), unit.out_bits)
    return a, b, y


def write_case(out, e, unit):
    n, ow = unit.n, unit.out_bits
    a, b, y = vectors(unit)
    two = unit.inputs == 2
    word = (mask_bits(a, n) << (n + ow if two else ow)) | ((mask_bits(b, n) << ow) if two else 0) | y
    wordw = (2 * n + ow) if two else (n + ow)
    hexw = (wordw + 3) // 4
    vec = os.path.join(out, e.key + f"_n{n}.hex")
    with open(vec, "w") as f:
        f.write("\n".join(format(int(v), f"0{hexw}x") for v in word) + "\n")
    plist = ", ".join(f".{k}({v})" for k, v in e.params.items())
    tbp = os.path.join(out, f"tb_{e.key}_n{n}.v")
    with open(tbp, "w") as f:
        f.write(TB.format(n=n, ow=ow, nv=len(a), mod=e.module, plist=plist, vec=os.path.abspath(vec),
                          wordw=hexw * 4, bsh=" + N" if two else "", ports=".a(a), .b(b), .y(y)" if two else ".a(a), .y(y)",
                          bexpr="mem[i][OW +: N]" if two else "0"))
    return vec, tbp


def run_case(out, e, unit, srcs=None):
    vec, tbp = write_case(out, e, unit)
    exe = os.path.join(out, f"{e.key}_n{unit.n}.vvp")
    c = subprocess.run(["iverilog", "-g2012", "-Wall", "-o", exe, tbp] + (srcs or e.srcs), capture_output=True, text=True)
    if c.returncode:
        return "COMPILE_ERROR " + (c.stderr.strip().splitlines() or [""])[0][:160], vec
    r = subprocess.run(["vvp", exe], capture_output=True, text=True, timeout=900)
    m = re.search(r"RESULT (PASS|FAIL).*", r.stdout)
    st = m.group(0) if m else "NO_RESULT"
    os.remove(exe)
    if "PASS" in st and os.path.exists(vec):
        os.remove(vec)
    return st, vec


MUTATIONS = [  # (entry key, rtl file, old, new)
    ("bw_trunc_k4_c0", "bw_mult.v", "if ((i == N-1) != (j == N-1)) t = ~t;", "if (0) t = ~t;"),
    ("bw_orlow_k4", "bw_mult.v", "y = (MODE == 1) ? (acc | orl) : acc;", "y = acc;"),
    ("sdrum_k4", "sdrum_mult.v", "if (sa > 0) fa[0] = 1'b1;", ""),
    ("smag_drum4_sat", "smag_mult.v", "assign ua = ma[N-1] ? {W{1'b1}} : ma[W-1:0];", "assign ua = ma[W-1:0];"),
    ("div_trunc_t2_c1", "div_trunc.v", "if (r >= {1'b0, b})", "if (r > {1'b0, b})"),
    ("div_mitchell", "div_mitchell.v", "d  = ((la - lb) << F) + fa - fb;", "d  = ((la - lb) << F) + fa + fb;"),
    ("div_drum_k3", "div_drum.v", "if (hb > 0) sb[0] = 1'b1;", ""),
    ("recip_mitchell", "recip_mitchell.v", "v = (3 << k) - a;", "v = (2 << k) - a;"),
    ("sq_trunc_k6_c1", "sq_trunc.v", "compval = (e4 + 2) / 4;", "compval = e4 / 4;"),
    ("sq_mitchell", "sq_mitchell.v", "L2 = ((l << F) + f) << 1;", "L2 = ((l << F) + f);"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="build")
    ap.add_argument("--widths", default="8,6,5")
    ap.add_argument("--only", default="")
    ap.add_argument("--mutations", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    results = []
    if a.mutations:
        n = 5
        cat = {e.key: e for e in catalog(n)}
        nbad = 0
        for key, fn, old, new in MUTATIONS:
            e = cat[key]
            src = open(os.path.join(V2_RTL, fn)).read()
            assert old in src, (key, old)
            mdir = os.path.join(a.out, "mut"); os.makedirs(mdir, exist_ok=True)
            mp = os.path.join(mdir, fn); open(mp, "w").write(src.replace(old, new))
            srcs = [mp if os.path.basename(s) == fn else s for s in e.srcs]
            st, _ = run_case(mdir, e, e.unit(), srcs)
            ok = "FAIL" in st
            nbad += (not ok)
            print(f"{'DETECTED' if ok else 'MISSED  '} mutation of {fn:16s} on {key:18s}: {st[:60]}", flush=True)
        print(f"\nMUTATION SUMMARY: {len(MUTATIONS) - nbad}/{len(MUTATIONS)} detected")
        sys.exit(1 if nbad else 0)
    t0 = time.time()
    for n in [int(w) for w in a.widths.split(",")]:
        for e in catalog(n):
            if a.only and a.only not in e.key:
                continue
            u = e.unit()
            st, _ = run_case(a.out, e, u)
            results.append((f"{e.key}@n{n}", st))
            print(f"{st[:48]:48s} {e.key} n={n} [{e.module}]", flush=True)
    npass = sum("PASS" in s for _, s in results)
    nfail = len(results) - npass
    print(f"\nSUMMARY: {len(results)} cases, {npass} PASS, {nfail} FAIL/ERROR  ({time.time() - t0:.0f} s)")
    json.dump(results, open(os.path.join(a.out, "rtl_results.json"), "w"), indent=1)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
