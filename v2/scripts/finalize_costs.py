#!/usr/bin/env python3
"""Turn the incremental *.partial.csv (JSON lines written by synth_costs.py) into the final CSV summaries
(dedupe by key, last row wins).  Used because the synthesis run was time-boxed and stopped before its own end-of-run write.
  python scripts/finalize_costs.py results/costs_v2_n8.partial.csv results/costs_v2_n8.csv"""
import csv, json, sys
rows = {}
for l in open(sys.argv[1]):
    r = json.loads(l); rows[r["key"]] = r
keys = ["source", "key", "module", "group", "role", "unit", "n"]
rest = []
for r in rows.values():
    for k in r:
        if k not in keys and k not in rest:
            rest.append(k)
with open(sys.argv[2], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys + rest, restval=""); w.writeheader(); w.writerows(rows.values())
print(len(rows), "rows ->", sys.argv[2])
