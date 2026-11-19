"""Transparent gate-equivalent (GE) cost proxy.

THIS IS A PROXY, NOT SYNTHESIZED AREA.  Every number below is a round,
textbook-style relative weight (NAND2 = 1 GE).  It is not taken from any
specific standard-cell library, and the structural formulas (ripple adders,
carry-save bit counting, barrel shifters) are deliberate simplifications.
Use it to *rank* designs inside one family or across families at the same
bit width, not to quote absolute area.  Cross-check with a real synthesis
flow (see scripts/rtl_regress.py --synth) before drawing hardware conclusions.
"""
import math

GE = {
    "INV": 0.5, "NAND2": 1.0, "NOR2": 1.0, "AND2": 1.5, "OR2": 1.5,
    "XOR2": 2.5, "MUX2": 2.5,
    "HA": 4.0,    # XOR2 + AND2
    "FA": 9.5,    # 2*XOR2 + 2*AND2 + OR2
    "MAJ": 6.0,   # AND2 + OR2 + AND2 + OR2
}


def clog2(x: int) -> int:
    return 0 if x <= 1 else int(math.ceil(math.log2(x)))


def ripple_adder_ge(bits: int, carry_in: bool = False) -> float:
    """Ripple-carry adder: HA at the bottom (or FA if there is a carry-in)."""
    if bits <= 0:
        return 0.0
    if carry_in:
        return bits * GE["FA"]
    return GE["HA"] + (bits - 1) * GE["FA"]


def barrel_shifter_ge(width: int, max_shift: int) -> float:
    """Log shifter: ceil(log2(max_shift+1)) stages of `width` MUX2."""
    return width * clog2(max_shift + 1) * GE["MUX2"]


def lod_ge(n: int) -> float:
    """Leading-one detector + priority encoder, heuristic 3 GE per input bit."""
    return 3.0 * n


def carry_save_ge(m_bits: int, n_cols: int, pp_gen_ge: float) -> float:
    """Partial-product array cost from bit counting.

    `m_bits` partial-product bits are generated (pp_gen_ge each) and reduced
    to one bit per column over `n_cols` columns; every full adder removes one
    bit, so FA count = m_bits - n_cols.  Wiring and sign-extension are ignored.
    """
    if m_bits <= 0:
        return 0.0
    return m_bits * pp_gen_ge + max(m_bits - n_cols, 0) * GE["FA"]


def exact_array_mult_ge(n: int) -> float:
    return carry_save_ge(n * n, 2 * n, GE["AND2"])
