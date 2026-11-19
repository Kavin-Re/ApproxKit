// ApproxKit: radix-4 Booth multiplier, partial-product magnitudes truncated below column K (K=0 exact)
// Bit-accurate to approxkit.BoothApprox (see README). Verilog-2001, no vendor primitives.
module booth_mult #(parameter N = 8, parameter K = 4) (input [N-1:0] a, b, output reg [2*N-1:0] y);
  localparam D = N / 2 + 1;
  localparam [2*N-1:0] MASK = {(2*N){1'b1}} << K;
  reg [2*D:0]   bb;
  reg [2:0]     w;
  reg [2*N-1:0] A, mag;
  integer i;
  always @* begin
    A  = {{N{1'b0}}, a};
    bb = {{(2*D-N){1'b0}}, b, 1'b0};
    y  = 0;
    for (i = 0; i < D; i = i + 1) begin
      w = bb[2*i+2 -: 3];                 // {b[2i+1], b[2i], b[2i-1]}
      case (w)
        3'b001, 3'b010: mag = A;
        3'b011:         mag = A << 1;
        3'b100:         mag = A << 1;
        3'b101, 3'b110: mag = A;
        default:        mag = 0;
      endcase
      mag = (mag << (2*i)) & MASK;
      if (w == 3'b100 || w == 3'b101 || w == 3'b110) y = y - mag;
      else                                            y = y + mag;
    end
  end
endmodule
