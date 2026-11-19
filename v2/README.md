# ApproxKit v2: signed multipliers, dividers, squarers and application quality

A follow-on study built on the core library in [`../approxkit`](../approxkit). It imports
the core unchanged (the core is never modified from here) and adds units and experiments
that the unsigned core does not cover.

| | Core (`approxkit/`) | v2 (`v2/approxkit_v2/`) |
|---|---|---|
| Units | Unsigned adders and multipliers | Signed multipliers (Baugh-Wooley truncated / OR-low, signed DRUM, sign-magnitude wrappers around the core units), unsigned dividers, reciprocal, squarers |
| Metrics | ER, MED, MRED, NMED, WCE, bias, RMSE | Same, plus MRED16 and 95 % confidence intervals for Monte Carlo |
| Cost | Gate-equivalent proxy, yosys generic gates / sky130 | yosys iCE40 / ECP5 LUT and carry counts, rank correlation against sky130 |
| Applications | FIR, dot product, matmul | Image convolution (PSNR / SSIM), FIR, dot product, with confidence intervals |
| Installed by `pip install` | Yes | No; run from this directory |

## Layout

| Path | Contents |
|---|---|
| `approxkit_v2/` | `signed.py`, `divsq.py`, `metrics2.py`, `quality.py`, `catalog.py` (the 75-design n = 8 catalogue), `common.py` |
| `rtl/` | Verilog for every v2 unit |
| `scripts/` | `error_sweep.py`, `rtl_regress_v2.py`, `synth_costs.py`, `fpga_rank.py`, `app_quality.py`, `pareto_plots.py`, `eta_vs_trunc.py` |
| `tests/test_v2.py` | Scalar bit-level references, structural properties, metric checks |
| `results/` | CSVs, logs and plots behind every claim |
| [`RESULTS_V2.md`](RESULTS_V2.md) | Findings, a claim-to-script-to-evidence table, and limits |

## Run

```bash
python -m pytest v2/tests                 # from the repository root
cd v2
python scripts/error_sweep.py --n 8 --out results
python scripts/rtl_regress_v2.py --out build              # needs iverilog
python scripts/rtl_regress_v2.py --out build --mutations  # testbench must fail on injected faults
APPROXKIT_LIBERTY=/path/to/sky130_fd_sc_hd__tt_025C_1v80.lib \
  python scripts/synth_costs.py --n 8 --out results       # needs yosys; liberty optional
```

All hardware cost numbers are estimates from generic yosys flows; none is a vendor-tool or
board result. Read `RESULTS_V2.md` section 6 ("Limits and not done") before reusing them.
