#!/usr/bin/env python3
"""Compare two CSV files cell by cell; numeric cells may differ by a relative
tolerance (floating-point summation order can change the last digits between
numpy builds), everything else must match exactly.

  python scripts/compare_csv.py new.csv results/dse_adder_n8.csv [--rtol 1e-9]

Exit status 0 if equivalent, 1 otherwise.  Standard library only.
"""
import argparse
import csv
import math
import sys


def same(x, y, rtol):
    if x == y:
        return True
    try:
        fx, fy = float(x), float(y)
    except ValueError:
        return False
    if math.isnan(fx) and math.isnan(fy):
        return True
    return math.isclose(fx, fy, rel_tol=rtol, abs_tol=1e-12)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("new")
    ap.add_argument("reference")
    ap.add_argument("--rtol", type=float, default=1e-9)
    a = ap.parse_args()
    with open(a.new, newline="") as f:
        new = list(csv.reader(f))
    with open(a.reference, newline="") as f:
        ref = list(csv.reader(f))
    if len(new) != len(ref):
        print(f"row count differs: {len(new)} vs {len(ref)}")
        return 1
    bad = 0
    for i, (rn, rr) in enumerate(zip(new, ref)):
        if len(rn) != len(rr) or not all(same(x, y, a.rtol) for x, y in zip(rn, rr)):
            bad += 1
            if bad <= 5:
                print(f"row {i} differs:\n  new: {rn}\n  ref: {rr}")
    print(f"{a.new} vs {a.reference}: {len(ref)} rows, {bad} differing")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
