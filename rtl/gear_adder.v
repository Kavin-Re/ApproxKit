// ApproxKit: GeAr / ACA-style block-carry adder; R>=1, P>=0, R+P<=N
// Bit-accurate to approxkit.GeAr (see README). Verilog-2001, no vendor primitives.
module gear_adder #(parameter N = 8, parameter R = 2, parameter P = 2) (input [N-1:0] a, b, output [N:0] y);
  localparam L  = R + P;
  localparam NB = (N - L + R - 1) / R + 1;       // number of sub-adders
  genvar j;
  generate for (j = 0; j < NB; j = j + 1) begin : blk
    localparam S  = (j * R < N - L) ? j * R : N - L;                    // window start
    localparam SP = (j == 0) ? 0 : (((j-1) * R < N - L) ? (j-1) * R : N - L);
    localparam PT = (j == 0) ? 0 : SP + L;                              // first new bit
    localparam W  = S + L - PT;                                         // new bits
    wire [L:0] t = {1'b0, a[S +: L]} + {1'b0, b[S +: L]};
    assign y[PT +: W] = t[(PT - S) +: W];
    if (j == NB - 1) begin : last
      assign y[N] = t[L];
    end
  end endgenerate
endmodule
