// ApproxKit: simplified-FA-cell lower part; CELL=0: Sum=b,Cout=a ; CELL=1: Cout=maj, Sum=~Cout; 1 <= K <= N-1
// Bit-accurate to approxkit.AMA (see README). Verilog-2001, no vendor primitives.
module ama_adder #(parameter N = 8, parameter K = 3, parameter CELL = 0) (input [N-1:0] a, b, output [N:0] y);
  wire [K:0]   c;
  wire [K-1:0] lo;
  assign c[0] = 1'b0;
  genvar j;
  generate for (j = 0; j < K; j = j + 1) begin : g
    if (CELL == 0) begin : c0
      assign lo[j]  = b[j];
      assign c[j+1] = a[j];
    end else begin : c1
      assign c[j+1] = (a[j] & b[j]) | (c[j] & (a[j] | b[j]));
      assign lo[j]  = ~c[j+1];
    end
  end endgenerate
  wire [N-K:0] hi = {1'b0, a[N-1:K]} + {1'b0, b[N-1:K]} + c[K];
  assign y = {hi, lo};
endmodule
