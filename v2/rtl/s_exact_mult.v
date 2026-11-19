// ApproxKit v2: exact signed multiplier reference (behavioural '*'), N x N -> 2N two's complement
module s_exact_mult #(parameter N = 8) (input [N-1:0] a, b, output [2*N-1:0] y);
  wire signed [N-1:0] sa = a, sb = b;
  wire signed [2*N-1:0] p = sa * sb;
  assign y = p;
endmodule
