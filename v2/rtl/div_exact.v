// ApproxKit v2: exact divider reference, q = floor(a/b), b = 0 -> all ones
module div_exact #(parameter N = 8) (input [N-1:0] a, b, output [N-1:0] y);
  assign y = (b == 0) ? {N{1'b1}} : a / b;
endmodule
