// ApproxKit: truncated multiplier: partial-product columns < K dropped (bam with HBL=0)
// Bit-accurate to approxkit.TruncatedMultiplier (see README). Verilog-2001, no vendor primitives.
module trunc_mult #(parameter N = 8, parameter K = 6) (input [N-1:0] a, b, output [2*N-1:0] y);
  bam_mult #(.N(N), .VBL(K), .HBL(0)) u (.a(a), .b(b), .y(y));
endmodule
