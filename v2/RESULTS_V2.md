# ApproxKit v2 - results, claim-to-script map, caveats

Everything below was produced by the scripts in `scripts/` (run from `v2/`; yosys 0.33, iverilog 12.0, Python/numpy).
Raw numbers are in `results/*.csv`; this file contains no hand-typed measurement that is not in a CSV or log.
All hardware cost numbers are **estimates from generic yosys flows**; **nothing here is a PolarFire result**.

## 1. Claim -> script -> evidence

| # | Claim | Script | Evidence file | CI / exactness | Caveat |
|---|-------|--------|---------------|----------------|--------|
| C1 | Python models of all 75 n=8 designs (signed BW trunc/OR-low, signed DRUM, sign-magnitude wrappers, dividers, reciprocal, squarers) are bit-identical to the Verilog | `scripts/rtl_regress_v2.py` (iverilog, self-checking bench, ALL operand combinations incl. b=0 / x=0) | `results/rtl_regress_v2.log`: 285 cases (n=8:75, n=7:74, n=6:71, n=5:65), 285 PASS, 0 FAIL | exhaustive, no sampling | The RTL regression log is from the first part of the build and was not re-run afterwards. Mutation mode (`--mutations`) exists but its result is not part of this evidence. |
| C2 | Models agree with independent scalar bit-by-bit references (BW matrix, OR-low, signed DRUM, literal restoring-division loop, squarer matrix) over all operand pairs at n=5,6 | `tests/test_v2.py` | `results/test_v2.log`: 15 tests, ALL PASS | exhaustive | the Mitchell/DRUM dividers and Mitchell squarer/reciprocal have no second independent reference, only structural checks (exact on powers of two, bounds, b=0 convention) plus the RTL equivalence of C1 |
| C3 | Error metrics (ER, MED, MRED, MRED16, NMED, WCE, bias, RMSE) for every design, exhaustive at n=8 | `scripts/error_sweep.py`, `approxkit_v2/metrics2.py` | `results/errors_v2_n8.csv` | exact (CI half-width 0); Monte-Carlo path (16 batch means, t-95 %) is only used above 2^22 combinations and is checked against the exhaustive value in `test_monte_carlo_matches_exhaustive_within_ci` | MRED is over exact != 0; MRED16 over |exact| >= 16 to remove floor-quantisation artefacts of small quotients; divider metrics exclude b=0 |
| C4 | Generic-flow FPGA cost estimates and their rank agreement with sky130 | `scripts/synth_costs.py`, `scripts/fpga_rank.py` | `results/costs_*.csv`, `results/rank_corr.csv` | Spearman/Kendall with 95 % bootstrap CI over designs | see section 3 |
| C5 | Application quality (image convolution PSNR/SSIM, FIR SNR, dot-product SNR) and Pareto fronts for the signed multipliers | `scripts/app_quality.py`, `scripts/pareto_plots.py`, `approxkit_v2/quality.py` | `results/app_quality_v2.csv`, `results/pareto_apps_smul.png`, `results/pareto_fronts.json` | 95 % Student-t over 8 synthetic images (image metrics) / 8 batches (FIR, dot) | see section 4 |
| C6 | Audit observation: ETA and truncated-adder accumulators give identical dot-product errors because the low bits of both accumulators saturate and cancel | `scripts/eta_vs_trunc.py`, `tests/test_v2.py::test_eta_equals_truncation_accumulator_mechanism` | `results/eta_vs_trunc_saturation.csv`, `results/eta_vs_trunc_demo.csv` | 400 trials (saturation curve), 20 batches x 50 vectors (error table) | see section 5 |

## 2. What was built and what passes

* Signed multipliers (n=8): native Baugh-Wooley with truncated or OR-reduced low columns (`bw_trunc_k*_c{0,1}`, `bw_orlow_k*`), native signed DRUM (`sdrum_k*`, two's-complement segments), sign-magnitude wrappers around the 61 v1 unsigned units (`smag_*`, with optional saturation of -2^(n-1)). 37 designs.
* Unsigned dividers (truncated restoring `div_trunc_t*_c*`, Mitchell, segment/"DRUM-style" `div_drum_k*` - my own construction, not a published divider), Mitchell reciprocal, squarers (truncated folded matrix, DRUM, Mitchell). All with Python model + Verilog.
* Evidence: RTL<->model exhaustive equivalence 285/285 (C1), scalar-reference tests 15/15 (C2), `results/test_v2.log`.
* Error metrics, exhaustive at n=8: `results/errors_v2_n8.csv` (75 designs). Examples (NMED / WCE): `bw_trunc_k4_c1` 4.8e-4 / 37, `bw_orlow_k4` 2.0e-4 / 34, `sdrum_k5` 1.4e-2 / 1984, `div_mitchell` 2.6e-4 / 11, `recip_mitchell` 1.5e-3 / 2731 (bias +97), `sq_trunc_k6_c1` 4.0e-4 / 91, `sq_mitchell` 1.2e-2 / 4096 (bias -780).
* The compensated (c1) variants roughly halve MRED16 and remove most of the bias of plain truncation at equal area (e.g. k=4: bias -12.25 -> -0.25, MRED16 1.8 % -> 1.1 %) - from `errors_v2_n8.csv`.

## 3. Generic-vendor-flow FPGA cost vs sky130 ranking (C4)

**Labelled estimates, NOT PolarFire.** Flows: `synth_ice40`, `synth_ecp5 -nodsp`, (`synth_xilinx` is implemented in the script behind `--xilinx` but was NOT run for the reported tables - too slow on the build machine), all yosys 0.33, `-flatten`, DSP inference disabled. "LUT4" columns count LUTs only; carry cells (SB_CARRY / CCU2C) are in separate columns and are not in the LUT count. `ice40nc_lut4` = iCE40 with `-nocarry` (adders mapped to LUTs).
Reference ranking: v1 flow sky130 (yosys + abc + sky130_fd_sc_hd tt, `results/rtl_synth_n8.csv` of v1; re-synthesis here reproduces it exactly, rho = 1.000 on the 25 v1 multipliers).

| set (n designs) | iCE40 LUT4 | iCE40 -nocarry | ECP5 LUT4 | generic gate-equivalents |
|---|---|---|---|---|
| v1 multipliers (25) | rho 0.50 [0.02, 0.82] | not run | rho 0.55 [0.13, 0.80] | 0.99 [0.95, 1.00] |
| v2 all (71) | 0.73 [0.62, 0.81] | 0.99 [0.97, 0.99] | 0.73 [0.59, 0.82] | 1.00 [0.99, 1.00] |
| v2 signed multipliers (37) | 0.53 [0.21, 0.73] | 0.97 [0.92, 0.98] | 0.65 [0.36, 0.82] | 0.99 [0.96, 1.00] |
| v2 dividers (16) | 0.93 [0.74, 1.00] | 0.99 [0.91, 1.00] | 0.98 [0.85, 1.00] | 0.99 [0.92, 1.00] |
| v2 squarers + recip (16) | 0.55 [0.01, 0.87] | 0.97 [0.85, 1.00] | 0.39 [-0.16, 0.83] | 0.98 [0.87, 1.00] |

(Spearman rho, 95 % bootstrap CI over designs; Kendall tau-b is in `results/rank_corr.csv`.)

Findings: (1) A LUT count that ignores the carry chain tracks sky130 area only moderately (rho 0.5-0.7 for multipliers); once adders are forced into LUTs (`-nocarry`) or gates are counted, agreement is ~0.97-0.99. So "which design is cheaper" on a carry-chain FPGA depends heavily on how much of the saving sits in carry-chain adders, and sky130 area is an unreliable proxy for it. (2) Shifter-heavy designs (DRUM, Mitchell, segment divider) look much worse in LUT count than in sky130 area (e.g. `sdrum_k4`: 1260 um2 vs 498 LUT4 + 216 carry, against `bw_trunc_k4_c1` 1703 um2 vs 137 LUT4 + 8 carry; exact `s_exact` 2381 um2 / 181 LUT4). Barrel shifters cost LUTs, not carry cells. (3) None of this says anything about PolarFire (LUT4 + dedicated carry, different mapper, different DSP/LSRAM math blocks).

Limits: (a) v1 side is the 25 v1 multipliers only (v1 adders were not synthesised in this time-box); (b) 4 v2 designs have no cost: `div_trunc_t0` (yosys 0.33 aborts silently on `q[T-1]` with T=0 - it is the exact divider; `div_exact` is the exact reference), `div_drum_k2/k5/k6` (the 2k-by-k `/` unrolls to thousands of cells; per-tool timeout 140 s or run stopped for time; `div_drum_k3` took 313 s and gives 1671 LUT4 + 1168 carry, which is largely an artefact of how yosys maps `/`, so the DRUM-divider FPGA cost is not trustworthy); `div_drum_k4` finished but see the same caveat. (c) Mixed-effects: rank correlations across heterogeneous families are driven by between-family gaps; within-family rankings were not analysed separately.

## 4. Application quality (C5)

Setup in `scripts/app_quality.py`: 8 synthetic 96x96 images (generated in the script), 3x3 kernels with signed int8 taps so every tap goes through the unit under test: `blur` (Q7 Gaussian, power-of-two taps), `blur_odd` (taps 5/21/24 - no power-of-two taps, deliberately), `sharpen`, `sobelx`; FIR (31 taps, 8 signals x 1470 outputs) and dot product (8 batches x 32 vectors x 256 terms). Reference = exact signed multiplier output; PSNR capped at 99 dB when outputs are identical; SSIM implemented in `approxkit_v2/quality.py` (Gaussian 11x11, sigma 1.5, cross-checked against scikit-image to 2e-3 in `test_quality_metrics`).
Pareto fronts: `results/pareto_apps_smul.png`, `results/pareto_fronts.json` (cost = sky130 area, or iCE40 LUT4 - labelled generic).

* On the Pareto front in at least one of the 6 panels: only the Baugh-Wooley truncation / OR-low family (`bw_trunc_k{2,4,6,8,12}_c1`, `bw_orlow_k6`, `bw_trunc_k4_c0`) plus, in some panels, `smag_trunc8(_sat)`, `smag_bam_6_2`, `smag_booth8` and (sky130-cost panels only) `sdrum_k2/k4`. No signed-DRUM point is on any iCE40-cost front. The leftmost front point (`bw_trunc_k12_c1`, near-zero quality) is just the cheapest design.
* Quality at moderate cost (`blur_odd`, mean over 8 images, 95 % CI): `bw_trunc_k4_c1` 56.7 +/- 3.2 dB, `bw_trunc_k6_c1` 43.8 +/- 3.0 dB, `sdrum_k5` 38.7 +/- 2.2 dB, `smag_drum5` 47.5 +/- 5.8 dB; dot-product SNR: `bw_trunc_k4_c1` 54.6 dB, `sdrum_k5` 23.8 dB, `smag_mitchell` 26.8 dB (see CSV for CIs; FIR/dot CIs are over 8 batches).
* The compensation constant (c1) is worth about 3-25 dB depending on metric and k (`bw_trunc_k4`: `blur_odd` 48.2 -> 56.7 dB, FIR 24.9 -> 33.3 dB, dot 29.5 -> 54.6 dB; `bw_trunc_k6`: 35.3 -> 43.8, 8.9 -> 16.0, 13.1 -> 41.1).
* Caveat that changes conclusions: with power-of-two kernel taps (`blur`, `sobelx`, `sharpen` mostly), Mitchell and OR-low units are exact (PSNR = 99 cap) because one operand is a power of two; `blur_odd` was added to avoid that flattering result. The "PSNR 99" rows are identical outputs, not measured PSNR.
* Caveat: 8 synthetic images only (no natural photographs, no external data); CIs reflect image-to-image spread of this synthetic set, not a population of natural images. Accumulation is exact; approximate adders are not combined with approximate multipliers here.

## 5. ETA vs truncation identical dot-product error (C6)

Reproduced and explained (v1 `results/apps_accum_adders.csv`: dot rel-RMS 0.001424 for both at k=4, 0.026805 at k=8, 0.38954 at k=12).
Mechanism: v1's dot product keeps two UNSIGNED accumulators (positive and negative products) and subtracts them. In ETA, once accumulator bit k-1 and addend bit k-1 are both 1 (probability about 1/4 per addition), all lower result bits become 1, and the all-ones state is absorbing. After a few additions both accumulators have low k bits = 2^k - 1, which cancel in the difference; truncation has 0 in both. Upper parts are identical ((a>>k)+(b>>k)). `results/eta_vs_trunc_saturation.csv` (k=4, 400 trials): fraction of identical results = 0.20 (1 term), 0.05 (2), 0.57 (8), 0.96 (16), 1.00 (32+). Tests: `test_eta_equals_truncation_accumulator_mechanism` (256 terms: both low nibbles == 15 in 100 % of 200 trials, identical outputs 100 %; 2 terms: identical only 4.5 %; single accumulator: differs 100 %).
Corrected demo (`scripts/eta_vs_trunc.py`, `results/eta_vs_trunc_demo.csv`, rel-RMS +/- 95 % CI, 20 batches x 50 vectors, 256 terms):

| k | organisation | LOA | truncation | ETA |
|---|---|---|---|---|
| 4 | two accumulators (v1) | 0.00088 +/- 0.00004 | 0.00145 +/- 0.00009 | 0.00145 +/- 0.00009 (bit-identical) |
| 4 | single 24-bit two's-complement | 0.0031 +/- 0.0002 | 0.0196 +/- 0.0010 | 0.0194 +/- 0.0010 |
| 8 | two accumulators (v1) | 0.0137 +/- 0.0007 | 0.0270 +/- 0.0009 | 0.0270 +/- 0.0009 (bit-identical) |
| 8 | single accumulator | 0.0148 +/- 0.0011 | 0.384 +/- 0.019 | 0.381 +/- 0.019 |

Interpretation: the identity is an artefact of the two-unsigned-accumulator organisation, not a property of the adders. In a single two's-complement accumulator ETA and truncation are no longer bit-identical, but their RMS errors are statistically indistinguishable (CIs overlap) and both are far worse than LOA, because negative addends are floored in truncation (bias) and ETA saturates to ones. The v1 conclusion "ETA ~ truncation, LOA better" holds, but the v1 evidence for ETA was uninformative, and the cost-benefit of ETA (which is larger than truncation: 212 vs 184 GE at k=4 in v1's proxy) is therefore negative in this application.

## 6. Limits and not done

* No PolarFire/Libero results; vendor-flow numbers are generic yosys estimates (section 3). `synth_xilinx` not run for the tables.
* Cost analysis covers 71 of 75 v2 designs and 25 of 61 v1 designs (v1 adders omitted for time).
* Mutation testing of the RTL bench (`rtl_regress_v2.py --mutations`) has not been run; the claim that the bench detects broken RTL is untested.
* The RTL regression log is from the first part of the build (not re-run); the Python-side tests were re-run (`results/test_v2.log`).
* Only n=8 cost/quality; n=5-7 used for equivalence testing only. No formal (z3/yosys-smtbmc) equivalence.
* 8 synthetic images only; no natural-image evaluation; no PSNR on colour or video.
* Unit-model RTL for dividers uses behavioural `/` (restoring-loop for `div_trunc`), so FPGA/ASIC divider costs depend on yosys's mapping of `/` and `%`.
