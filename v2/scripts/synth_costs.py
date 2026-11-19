#!/usr/bin/env python3
"""Synthesis cost table: sky130 (yosys+abc, liberty) and three GENERIC vendor flows
of the installed yosys (iCE40, ECP5, Xilinx-7 style).  NONE of these is PolarFire.

  python scripts/synth_costs.py --n 8 --liberty /path/sky130hd_tt_du.lib --out results

Outputs  results/costs_v2_n8.csv   (v2 designs, from approxkit_v2.catalog)
         results/costs_v1_n8.csv   (the 61 v1 designs, same RTL, re-synthesised here)
Columns: sky130_um2 | gen_ge (generic gates, v1 GE weights) | ice40_lut4, ice40_carry |
         ice40nc_lut4 (iCE40 with -nocarry: carry chains disabled) | ecp5_lut4, ecp5_ccu2c |
         xil_lut (LUT1-6), xil_carry4, xil_muxf.
Flows: ice40 'synth_ice40 -flatten', ecp5 'synth_ecp5 -nodsp -flatten', xilinx
'synth_xilinx -family xc7 -nodsp -flatten' (DSP inference disabled so '*' is built from fabric).
"""
import argparse, csv, json, threading, importlib.util, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from approxkit_v2.catalog import catalog
from approxkit_v2.common import V1_ROOT

GE = {"$_AND_": 1.5, "$_NAND_": 1.0, "$_OR_": 1.5, "$_NOR_": 1.0, "$_XOR_": 2.5, "$_XNOR_": 2.5,
      "$_ANDNOT_": 1.5, "$_ORNOT_": 1.5, "$_NOT_": 0.5, "$_MUX_": 2.5}


def run_yosys(srcs, module, params, tail, timeout=140):
    cp = " ".join(f"-set {k} {v}" for k, v in params.items())
    script = f"read_verilog {' '.join(srcs)}; chparam {cp} {module}; {tail}"
    r = subprocess.run(["yosys", "-p", script], capture_output=True, text=True, timeout=timeout)
    txt = r.stdout + r.stderr
    if r.returncode:
        raise RuntimeError("yosys failed: " + "\n".join(txt.splitlines()[-6:]))
    return txt


def last_stat(txt):
    blk = txt.rsplit("Printing statistics", 1)[-1]
    return {m.group(1): int(m.group(2)) for m in re.finditer(r"^\s+(\S+)\s+(\d+)\s*$", blk, re.M)
            if m.group(1) != "cells"}, blk


def costs_for(module, params, srcs, liberty, xil=True, nocarry=True):
    out = {}
    # sky130
    if liberty:
        t = run_yosys(srcs, module, params, f"synth -top {module} -flatten; abc -liberty {liberty}; opt_clean; stat -liberty {liberty}")
        m = re.search(r"Chip area for module.*?:\s*([0-9.]+)", t)
        out["sky130_um2"] = float(m.group(1)) if m else None
    # generic gates
    t = run_yosys(srcs, module, params, f"synth -top {module} -flatten; abc -g AND,NAND,OR,NOR,XOR,XNOR,ANDNOT,ORNOT,MUX; opt_clean; stat")
    c, _ = last_stat(t)
    out["gen_ge"] = sum(GE.get(k, 0) * v for k, v in c.items())
    out["gen_cells"] = sum(v for k, v in c.items() if k.startswith("$_"))
    # iCE40
    t = run_yosys(srcs, module, params, f"synth_ice40 -flatten -top {module}; stat")
    c, _ = last_stat(t)
    out["ice40_lut4"], out["ice40_carry"] = c.get("SB_LUT4", 0), c.get("SB_CARRY", 0)
    if nocarry:
        t = run_yosys(srcs, module, params, f"synth_ice40 -flatten -nocarry -top {module}; stat")
        c, _ = last_stat(t)
        out["ice40nc_lut4"] = c.get("SB_LUT4", 0)
    # ECP5
    t = run_yosys(srcs, module, params, f"synth_ecp5 -nodsp -flatten -top {module}; stat")
    c, _ = last_stat(t)
    out["ecp5_lut4"], out["ecp5_ccu2c"] = c.get("LUT4", 0), c.get("CCU2C", 0)
    if not xil:
        return out
    # Xilinx-7 style
    t = run_yosys(srcs, module, params, f"synth_xilinx -family xc7 -nodsp -flatten -top {module}; stat")
    c, _ = last_stat(t)
    out["xil_lut"] = sum(v for k, v in c.items() if re.fullmatch(r"LUT[1-6]", k))
    out["xil_carry4"] = c.get("CARRY4", 0)
    out["xil_muxf"] = sum(v for k, v in c.items() if k.startswith("MUXF"))
    return out


def v1_cases(n):
    spec = importlib.util.spec_from_file_location("v1_rtl_regress", os.path.join(V1_ROOT, "scripts", "rtl_regress.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    res = []
    for m, params, unit in mod.configs(n):
        srcs = [os.path.join(V1_ROOT, "rtl", m + ".v")] + ([os.path.join(V1_ROOT, "rtl", "bam_mult.v")] if m == "trunc_mult" else [])
        res.append((mod.tag(m, params), m, params, srcs, unit))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--liberty", default=os.environ.get("APPROXKIT_LIBERTY", ""),
                    help="sky130 (or other) liberty file for the mapped-area column; "
                         "defaults to $APPROXKIT_LIBERTY, column left empty if unset")
    ap.add_argument("--out", default="results")
    ap.add_argument("--which", default="v2,v1")
    ap.add_argument("--groups", default="", help="v2 only: comma list of groups (smul,div,recip,sq); default all")
    ap.add_argument("--v1-mult-only", action="store_true", help="time-box: only v1 multiplier designs")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--xilinx", action="store_true", help="also run synth_xilinx (slow: ~9 s/design); v2 designs only")
    a = ap.parse_args()
    lib = a.liberty if a.liberty and os.path.exists(a.liberty) else None
    if lib is None:
        print("WARNING: liberty not found; sky130 column left empty")
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    if "v2" in a.which:
        for e in catalog(a.n):
            if a.groups and e.group not in a.groups.split(','):
                continue
            jobs.append(("v2", e.key, e.module, e.params, e.srcs, e.group, e.role, e.unit().name))
    if "v1" in a.which:
        for t, m, p, s, u in v1_cases(a.n):
            if a.v1_mult_only and 'mult' not in m:
                continue
            jobs.append(("v1", t, m, p, s, u.kind, "exact" if "exact" in u.family or t.endswith("N8_K0") else "approx", u.name))

    lock = threading.Lock()

    def work(j):
        src, key, m, p, s, grp, role, uname = j
        t0 = time.time()
        try:
            c = costs_for(m, p, s, lib, xil=(a.xilinx and src == 'v2'), nocarry=(src == 'v2'))
        except Exception as ex:      # e.g. yosys 0.33 aborts silently on div_trunc T=0 (q[T-1] elaboration)
            print(f"{src} {key:28s} SYNTH FAILED: {str(ex)[-80:]!r}", flush=True)
            return dict(source=src, key=key, module=m, group=grp, role=role, unit=uname, n=a.n, synth_error="failed")
        row = dict(source=src, key=key, module=m, group=grp, role=role, unit=uname, n=a.n, **c)
        print(f"{src} {key:28s} sky={c.get('sky130_um2')} ice40={c['ice40_lut4']}+{c['ice40_carry']}c ecp5={c['ecp5_lut4']}+{c['ecp5_ccu2c']} xil={c.get('xil_lut')}  ({time.time()-t0:.1f}s)", flush=True)
        with lock, open(os.path.join(a.out, f"costs_{src}_n{a.n}.partial.csv"), "a") as pf:
            pf.write(json.dumps(row) + "\n")
        return row

    with ThreadPoolExecutor(a.jobs) as ex:
        rows = list(ex.map(work, jobs))
    rows = [r for r in rows if r is not None]
    allk = []
    allk.append('synth_error')
    [allk.append(k) for r in rows for k in r if k not in allk]
    os.makedirs(a.out, exist_ok=True)
    for src in ("v2", "v1"):
        rs = [r for r in rows if r["source"] == src]
        if rs:
            with open(os.path.join(a.out, f"costs_{src}_n{a.n}.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=allk, restval='', extrasaction='ignore'); w.writeheader(); w.writerows(rs)


if __name__ == "__main__":
    main()
