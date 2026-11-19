// ApproxKit v2: exact reciprocal floor(2^M / x) on M+1 bits, x = 0 -> all ones
module recip_exact #(parameter N = 8, parameter M = 8) (input [N-1:0] a, output [M:0] y);
  wire [M:0] num = {1'b1, {M{1'b0}}};
  assign y = (a == 0) ? {(M+1){1'b1}} : num / a;
endmodule
