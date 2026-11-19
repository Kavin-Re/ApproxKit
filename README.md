# ApproxKit

Bit-accurate Python models, error metrics, a Pareto design-space explorer and
Verilog implementations of approximate adders and multipliers. Written to be
reused across projects (correlators / MAC arrays, small neural-net datapaths,
filters). Dependencies: `numpy` (matplotlib only for plots). Plain Python
asserts for tests (pytest also works if installed). MIT licence (change if you prefer).

```
approxkit/      package: adders.py multipliers.py metrics.py cost.py dse.py apps.py plotting.py
examples/       run_dse.py (Pareto + plots), demo_apps.py (FIR / dot / matmul)
rtl/            Verilog-2001 units (one module per file)
scripts/        rtl_regress.py: Python models -> exhaustive vectors -> iverilog self-check -> yosys
tests/          test_units.py (scalar reference models, metrics, apps)
results/        CSVs, PNGs and logs produced by the runs described below
```

## Quick start

```python
import numpy as np, approxkit as ak
add = ak.LOA(8, k=3)                  # 8-bit lower-part-OR adder, 3 approximate bits
mul = ak.DRUM(8, k=4)                 # 8-bit DRUM with 4-bit segments
add(200, 100), mul(np.array([13, 200]), np.array([9, 77]))   # numpy-vectorised, unsigned

m = ak.evaluate(mul)                  # exhaustive for n <= 12, Monte Carlo above
m.ER, m.MED, m.MRED, m.NMED, m.WCE    # see definitions below
counts, edges = m.histogram(bins=61)  # signed-error distribution
mul.area_ge()                         # gate-equivalent PROXY (not synthesized area)

ak.make("multiplier", "booth", 16, k=10)   # factory by family name
```
Run `python examples/run_dse.py --n 8` (or `--n 16`), `python examples/demo_apps.py`,
`python tests/test_units.py`, `python scripts/rtl_regress.py --out build --synth` (needs iverilog, yosys).

## Units (all unsigned; adders n -> n+1 bits, multipliers n x n -> 2n bits)

| Class | Parameters | Definition (bit-exact, same in Python and Verilog) |
|---|---|---|
| `ExactAdder`, `ExactMultiplier` | n | references |
| `LOA` | k | low k sum bits = a OR b; carry into the exact upper part = a[k-1] AND b[k-1] |
| `TruncatedAdder` | k | upper n-k bits added exactly, low k result bits = 0, no carry from the low part |
| `ETA` (ETA-I flavour) | k | upper part exact with no carry-in; low part scanned from bit k-1 down: sum = a XOR b until the first bit where a=b=1, which and all lower bits become 1 |
| `GeAr` (ACA-style for R=1) | R, P | sub-adders of length L=R+P, each producing R new bits from P predictor bits (carry-in speculated 0); last block anchored at the MSB; L=n is exact |
| `AMA` (simplified FA cells) | k, cell | low k bits use a simplified full-adder cell (`sum_b`: Sum=b, Cout=a; `inv_carry`: Cout=maj, Sum=NOT Cout). These two cells are *my choice in the spirit of* mirror-adder approximations, not copies of a specific published cell |
| `TruncatedMultiplier` | k | partial-product bits in columns < k dropped, no compensation |
| `BrokenArray` | vbl, hbl | rows i < hbl dropped and columns < vbl dropped, rest summed exactly |
| `DRUM` | k | k-bit segment from the leading one with LSB forced to 1 if bits were discarded; exact k x k product, shifted back; operands < 2^k exact |
| `Mitchell` | - | log2 x ~ K + f, add in fixed point (n-1 fraction bits), linear antilog, floor. Always underestimates |
| `BoothApprox` | k | radix-4 Booth (unsigned input handled with n//2+1 digits); each partial-product *magnitude* has columns < k cleared, then it is negated if the digit is negative; sum mod 2^2n (k=0 exact) |

Constraints: `0 <= k <= n-1` for adders (Verilog: `1 <= k <= n-1`); n <= 31
(int64 arithmetic). Signed use: apply sign-magnitude around the unsigned unit
(`approxkit.apps.signed_mul`).

## Error metrics

With e = approx - exact over the evaluated pairs (uniform inputs):
**ER** = P(e != 0); **MED** = mean |e|; **MRED** = mean |e|/exact over pairs with exact != 0;
**NMED** = MED / (max exact output: 2(2^n-1) for adders, (2^n-1)^2 for multipliers);
**WCE** = max |e|; also reported: signed mean error (`bias`) and RMSE.
Exhaustive for n <= 12 (up to n = 14 can be forced); otherwise Monte Carlo
(default 2^20 seeded uniform pairs). Monte-Carlo WCE is only a *lower bound* of the true
worst case. Uniform inputs are rarely the real distribution of your application:
re-run `evaluate` on your own operand statistics (pass arrays through the unit directly)
before trusting MRED/NMED for a design.

## Cost model (a proxy)

`approxkit/cost.py`: NAND2 = 1 GE, AND/OR 1.5, XOR 2.5, MUX 2.5, HA 4, FA 9.5.
Adders are ripple structures; multipliers use carry-save bit counting
(#FA = kept partial-product bits - active columns); shifters/leading-one detectors are
heuristics. **This is a ranking proxy, not synthesized area**: no wiring, no timing, no sizing.
`delay_levels()` gives the longest carry chain for adders. Sanity check against
synthesis is in the next section; the proxy is clearly wrong for some structures
(it understates Mitchell's shifters and ETA's control logic, and overstates GeAr's cost
for the small overlaps).

## Results from this build

**Tests.** `tests/test_units.py`: 10/10 pass (every unit against an independent scalar bit-level
reference over all operand pairs at n=6, plus exactness limits, known values, Mitchell
never overestimates, metrics vs brute force, Monte-Carlo vs exhaustive, apps).

**Verilog vs Python (iverilog 12.0, exhaustive).** `scripts/rtl_regress.py`: 173/173 cases pass
(61 configurations at n = 8, the same families at n = 7 and n = 5, every one of the 2^(2n) operand pairs compared bit-for-bit, including odd widths and
non-divisible GeAr windows). The testbench was verified to fail when the RTL is deliberately broken
(two mutations tried). Log: `results/rtl_regress.log`.

**Synthesis (yosys 0.33, n = 8, `synth -flatten` then `abc`).** Two numbers per design in
`results/rtl_synth_n8.csv`: generic-gate count / GE-weighted count (same GE weights as the proxy),
and, mapped to the sky130 high-density liberty used on the build machine, the cell area in um^2
(yosys+abc only: no place&route, no timing). These are *synthesis estimates, not silicon area*.
Rank agreement with the proxy (Spearman): 0.94 (adders) / 0.92 (multipliers) against the sky130
area, 0.70 / 0.88 against generic gate counts. Important caveat: `bam_mult`, `trunc_mult` and
`booth_mult` accumulate partial products with a sequential chain of `+`, whereas `exact_mult` uses `*`
(yosys builds a macc/carry-save structure). The fair baseline for those three is `bam_mult` with
VBL=0/HBL=0 (2054 um^2), not `exact_mult` (1927 um^2). A hand-built compressor tree would shift every
multiplier number.

**Pareto highlights, n = 8, exhaustive, error = MRED** (sky130 area from synthesis above):

| design | MRED | NMED | sky130 area (um^2) | exact reference |
|---|---|---|---|---|
| LOA k=2 | 0.34 % | 0.12 % | 184 | exact adder 208 |
| ETA k=2 | 0.37 % | 0.14 % | 168 | |
| LOA k=4 | 1.5 % | 0.56 % | 141 | |
| truncated mult k=6 | 2.6 % | 0.12 % | 1292 | array-with-adder-chain 2054 |
| truncated mult k=8 | 9.8 % | 0.69 % | 905 | |
| Booth-approx k=8 | 4.8 % | 0.24 % | 1132 | |
| DRUM-4 | 5.9 % | 1.4 % | 1295 | |
| DRUM-6 | 1.3 % | 0.36 % | 1647 | |
| Mitchell | 3.8 % | 0.93 % | 1848 | no saving at n = 8 |

Findings (n = 8): the best low-error multiplier designs on both the proxy and the synthesized area are
the column-truncation/broken-array/Booth-approx family (DRUM reaches the synthesized Pareto front
only at isolated points: k = 2, 6, 7), and Mitchell gives no area benefit at 8 bits in this implementation (its barrel shifters and
detectors cost about as much as the exact multiplier), though its relative error is data-independent
and one-sided. Adders: LOA and ETA win at small k; the `sum_b` AMA cell is the cheapest at larger k
but errs on almost every input (ER ~ 94 % at k=4), ETA and LOA have similar MRED at equal k
(ETA's WCE is about 2x larger at k=4) and ETA's control chain costs more proxy area. GeAr buys delay (carry chain = L), not area: its cost is
higher than the exact ripple adder, and a missed carry costs a large error (2^L-scale), which
is visible in the accumulator demo. Full sweeps: `results/dse_*_n8.csv`, `results/dse_*_n16.csv`,
Pareto plots `results/pareto_*.png`, histograms `results/errhist_*.png`.
(The n = 16 sweep is Monte Carlo, 2^20 samples, seed 0.)

**Application demos** (`results/apps_*.csv/png`, 8-bit multipliers; 31-tap FIR on a noisy two-tone
signal, 256-long dot products, 16x16 uint8 matmul; accuracy vs exact): e.g. Mitchell gives 31 dB FIR SNR
(dot-product relative RMS 5 %, matmul 3.4 %); truncated k=1 gives 66 dB; DRUM-2 only 9.8 dB. With an exact
multiplier and 24-bit accumulators, LOA k=4 gives 47 dB FIR SNR vs 43.5 dB for truncated/ETA k=4, and
GeAr(R=4,P=8) is very poor on dot products (relative RMS 0.49) because missed carries are large. FIR operands here
span only 7 bits, so DRUM-7 is trivially exact there.
Sign handling: magnitudes go through the unsigned unit, accumulation uses separate positive and
negative accumulators.

## Verilog (`rtl/`)

`exact_adder`, `loa_adder`, `trunc_adder`, `eta_adder`, `gear_adder`, `ama_adder`,
`exact_mult`, `bam_mult`, `trunc_mult`, `drum_mult`, `mitchell_mult`, `booth_mult`.
All have ports `a`, `b`, `y` and parameters `N` plus the unit parameter names above (upper case).
Plain Verilog-2001 (`$clog2`), purely combinational, written for clarity rather than timing; the
multipliers use behavioural `for` loops and `+` chains that synthesis tools restructure.
Regenerate vectors/testbenches: `python scripts/rtl_regress.py --out build [--widths 8,7,5] [--no-sim] [--synth] [--liberty file.lib]`.

## Honest limits

* Unsigned only; signed operation via sign-magnitude wrappers (extra hardware not modelled).
* Cost = proxy; do not quote GE or the sky130 numbers as area for a PolarFire/other-FPGA flow. On an
  FPGA, carry chains and LUT mapping change the picture entirely (an exact adder is already fast and cheap
  there, and many approximations save nothing); synthesize on your actual target.
* No delay/power model beyond the adder carry-chain length.
* Metrics assume uniform independent operands; real signals (e.g. correlator inputs, NN activations) differ.
* Booth unit uses a sign-magnitude truncation definition of mine; the two `AMA` cells are mine; the
  other units follow the usual descriptions of each technique but exact variants in papers differ in details
  (e.g. compensation constants), so numbers will not match published tables.
* No literature numbers are quoted anywhere in this repository.

## References (only those I am sure of)

* J. N. Mitchell Jr., "Computer Multiplication and Division Using Binary Logarithms", IRE Trans. Electronic Computers, 1962.
* S. Hashemi, R. I. Bahar, S. Reda, "DRUM: A Dynamic Range Unbiased Multiplier for Approximate Applications", ICCAD 2015.

The remaining techniques (lower-part-OR adders, error-tolerant adders, generic accuracy-configurable
adders / almost-correct adders, mirror-adder approximations, broken-array multipliers, approximate Booth
multipliers) are described from general knowledge of the techniques without citation.
