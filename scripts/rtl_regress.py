#!/usr/bin/env python3
"""RTL regression: Python models -> exhaustive vectors -> self-checking Verilog
testbenches (iverilog) -> optional yosys gate counts.

  python scripts/rtl_regress.py --out build            # generate + simulate + (if yosys) synth
  python scripts/rtl_regress.py --out build --no-sim   # only write vectors/testbenches
  python scripts/rtl_regress.py --widths 8 --synth

Every testbench applies ALL 2^(2N) operand pairs and compares the DUT output
bit-for-bit with the Python model (exhaustive, not sampled).
Synthesis numbers are yosys generic-gate counts after `synth` + `abc`; they are
NOT technology-mapped area and use a plain-adder `+` structure.
"""
import argparse, os, re, subprocess, sys, shutil, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import approxkit as ak

RTL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rtl")


def configs(n):
    """(verilog_module, {param: value}, python_unit) for width n."""
    c = []
    ks = sorted({1, 2, n // 2, n - 2, n - 1} & set(range(1, n)))
    c.append(("exact_adder", {"N": n}, ak.ExactAdder(n)))
    for k in ks:
        c += [("loa_adder", {"N": n, "K": k}, ak.LOA(n, k)),
              ("trunc_adder", {"N": n, "K": k}, ak.TruncatedAdder(n, k)),
              ("eta_adder", {"N": n, "K": k}, ak.ETA(n, k))]
        for cell, name in ((0, "sum_b"), (1, "inv_carry")):
            c.append(("ama_adder", {"N": n, "K": k, "CELL": cell}, ak.AMA(n, k, name)))
    for R, P in [(1, 1), (1, n // 2), (2, 2), (2, 1), (3, 1), (3, 2), (n // 2, n // 2), (n - 1, 1), (n, 0), (1, n - 1)]:
        if R >= 1 and R + P <= n:
            c.append(("gear_adder", {"N": n, "R": R, "P": P}, ak.GeAr(n, R, P)))
    c.append(("exact_mult", {"N": n}, ak.ExactMultiplier(n)))
    for vbl, hbl in [(0, 0), (n // 2, 0), (n - 2, 2), (n, n // 2 - 1), (2 * n - 2, n // 2)]:
        c.append(("bam_mult", {"N": n, "VBL": vbl, "HBL": hbl}, ak.BrokenArray(n, vbl, hbl)))
    for k in (0, 1, n - 2, n, 2 * n - 4):
        c.append(("trunc_mult", {"N": n, "K": k}, ak.TruncatedMultiplier(n, k)))
    for k in range(2, n):
        c.append(("drum_mult", {"N": n, "K": k}, ak.DRUM(n, k)))
    c.append(("mitchell_mult", {"N": n}, ak.Mitchell(n)))
    for k in sorted({0, 1, 3, n - 3, n, n + 2, 2 * n - 4}):
        c.append(("booth_mult", {"N": n, "K": k}, ak.BoothApprox(n, k)))
    return c


TB = """`timescale 1ns/1ps
module tb;
  localparam N = {n};
  localparam OW = {ow};
  localparam NV = {nv};
  reg [N+N+OW-1:0] mem [0:NV-1];
  reg  [N-1:0] a, b;
  wire [OW-1:0] y;
  reg  [OW-1:0] expv;
  integer i, errs;
  {mod} #({plist}) dut (.a(a), .b(b), .y(y));
  initial begin
    $readmemh("{vec}", mem);
    errs = 0;
    for (i = 0; i < NV; i = i + 1) begin
      a = mem[i][OW+N +: N];
      b = mem[i][OW +: N];
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


def tag(mod, params):
    return mod + "_" + "_".join(f"{k}{v}" for k, v in params.items())


def write_case(out, mod, params, unit):
    n = unit.n
    a = np.repeat(np.arange(1 << n, dtype=np.int64), 1 << n)
    b = np.tile(np.arange(1 << n, dtype=np.int64), 1 << n)
    y = unit(a, b)
    ow = unit.out_bits
    word = (a << (n + ow)) | (b << ow) | y
    hexw = (n + n + ow + 3) // 4
    t = tag(mod, params)
    vec = os.path.join(out, t + ".hex")
    with open(vec, "w") as f:
        f.write("\n".join(format(int(v), f"0{hexw}x") for v in word) + "\n")
    plist = ", ".join(f".{k}({v})" for k, v in params.items())
    tbp = os.path.join(out, "tb_" + t + ".v")
    with open(tbp, "w") as f:
        f.write(TB.format(n=n, ow=ow, nv=len(a), mod=mod, plist=plist, vec=os.path.abspath(vec)))
    return t, vec, tbp


def _yosys(mod, params, tail):
    src = os.path.join(RTL, mod + ".v")
    srcs = [src] + ([os.path.join(RTL, "bam_mult.v")] if mod == "trunc_mult" else [])
    cp = " ".join(f"-set {k} {v}" for k, v in params.items())
    script = f"read_verilog {' '.join(srcs)}; chparam {cp} {mod}; synth -top {mod} -flatten; {tail}; opt_clean; stat"
    if "-liberty" in tail:
        script += " -liberty " + tail.split("-liberty ")[1].split()[0]
    r = subprocess.run(["yosys", "-p", script], capture_output=True, text=True)
    return r.stdout + r.stderr


def synth_gates(mod, params, liberty=None):
    """Generic-gate cell histogram (+ optional liberty-mapped area in um^2)."""
    txt = _yosys(mod, params, "abc -g AND,NAND,OR,NOR,XOR,XNOR,ANDNOT,ORNOT,MUX")
    cells = {m.group(1): int(m.group(2)) for m in re.finditer(r"^\s+(\$_[A-Z0-9]+_)\s+(\d+)\s*$", txt, re.M)}
    area = None
    if liberty:
        t2 = _yosys(mod, params, f"abc -liberty {liberty}")
        m = re.search(r"Chip area for module.*?:\s*([0-9.]+)", t2)
        area = float(m.group(1)) if m else None
    return cells, area


GE = {"$_AND_": 1.5, "$_NAND_": 1.0, "$_OR_": 1.5, "$_NOR_": 1.0, "$_XOR_": 2.5, "$_XNOR_": 2.5,
      "$_ANDNOT_": 1.5, "$_ORNOT_": 1.5, "$_NOT_": 0.5, "$_MUX_": 2.5}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="build")
    ap.add_argument("--widths", default="8,7,5")
    ap.add_argument("--no-sim", action="store_true")
    ap.add_argument("--synth", action="store_true", help="yosys gate counts (needs yosys)")
    ap.add_argument("--liberty", default=None, help="optional liberty file: also report mapped chip area (um^2)")
    ap.add_argument("--only", default="", help="substring filter on case tag")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    have_iv = shutil.which("iverilog") and shutil.which("vvp") and not a.no_sim
    have_y = shutil.which("yosys") and a.synth
    results, synth, seen = [], [], set()
    for n in [int(w) for w in a.widths.split(",")]:
        for mod, params, unit in configs(n):
            t = tag(mod, params)
            if (a.only and a.only not in t) or t + str(n) in seen:
                continue
            seen.add(t + str(n))
            t, vec, tbp = write_case(a.out, mod, params, unit)
            status = "NOT_SIMULATED"
            if have_iv:
                srcs = [os.path.join(RTL, mod + ".v")] + ([os.path.join(RTL, "bam_mult.v")] if mod == "trunc_mult" else [])
                exe = os.path.join(a.out, t + ".vvp")
                c = subprocess.run(["iverilog", "-g2012", "-Wall", "-o", exe, tbp] + srcs, capture_output=True, text=True)
                if c.returncode:
                    status = "COMPILE_ERROR " + c.stderr.strip().splitlines()[0][:120] if c.stderr.strip() else "COMPILE_ERROR"
                else:
                    r = subprocess.run(["vvp", exe], capture_output=True, text=True, timeout=600)
                    m = re.search(r"RESULT (PASS|FAIL).*", r.stdout)
                    status = m.group(0) if m else "NO_RESULT"
                    if "FAIL" in status:
                        status += " | " + " ; ".join(l for l in r.stdout.splitlines() if "MISMATCH" in l)[:200]
                    os.remove(exe)
                os.remove(vec) if "PASS" in status else None
            results.append((t, status))
            print(f"{status.split(' |')[0][:60]:34s} {t}", flush=True)
            if have_y and n == 8:
                cells, area = synth_gates(mod, params, a.liberty)
                ge = sum(GE.get(k, 0) * v for k, v in cells.items())
                synth.append({"case": t, "cells": sum(cells.values()), "ge_generic": ge,
                              "proxy_ge": unit.area_ge(), "liberty_area_um2": area, "cells_by_type": cells})
                print(f"    yosys: {sum(cells.values())} generic cells, {ge:.1f} GE-weighted | proxy {unit.area_ge():.1f}" + (f" | mapped {area} um^2" if area else ""), flush=True)
    npass = sum("PASS" in s for _, s in results)
    nfail = sum(("FAIL" in s or "ERROR" in s or "NO_RESULT" in s) for _, s in results)
    print(f"\nSUMMARY: {len(results)} cases, {npass} PASS, {nfail} FAIL/ERROR, "
          f"{sum(s == 'NOT_SIMULATED' for _, s in results)} not simulated")
    with open(os.path.join(a.out, "rtl_results.json"), "w") as f:
        json.dump({"results": results, "synth": synth}, f, indent=1)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
