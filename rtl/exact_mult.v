// ApproxKit: exact array reference multiplier (baseline)
// Bit-accurate to approxkit.ExactMultiplier (see README). Verilog-2001, no vendor primitives.
module exact_mult #(parameter N = 8) (input [N-1:0] a, b, output [2*N-1:0] y);
  assign y = a * b;
endmodule
