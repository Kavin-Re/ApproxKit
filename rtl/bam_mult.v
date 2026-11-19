// ApproxKit: broken-array multiplier: rows i<HBL dropped, columns <VBL dropped
// Bit-accurate to approxkit.BrokenArray (see README). Verilog-2001, no vendor primitives.
module bam_mult #(parameter N = 8, parameter VBL = 6, parameter HBL = 0) (input [N-1:0] a, b, output reg [2*N-1:0] y);
  localparam [2*N-1:0] MASK = {(2*N){1'b1}} << VBL;
  reg [2*N-1:0] pp;
  integer i;
  always @* begin
    y = 0;
    for (i = HBL; i < N; i = i + 1) begin
      pp = b[i] ? ({{N{1'b0}}, a} << i) : {(2*N){1'b0}};
      y  = y + (pp & MASK);
    end
  end
endmodule
