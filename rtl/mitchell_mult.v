// ApproxKit: Mitchell logarithmic multiplier (floor); N >= 2
// Bit-accurate to approxkit.Mitchell (see README). Verilog-2001, no vendor primitives.
module mitchell_mult #(parameter N = 8) (input [N-1:0] a, b, output reg [2*N-1:0] y);
  localparam F  = N - 1;
  localparam LW = $clog2(2*N) + F + 1;
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin
      msb = 0;
      for (m = 0; m < N; m = m + 1) if (x[m]) msb = m;
    end
  endfunction
  integer la, lb, kk;
  reg [LW-1:0]  fa, fb, lg;
  reg [2*N-1:0] val;
  always @* begin
    la = msb(a);  lb = msb(b);
    if (a == 0 || b == 0) y = 0;
    else begin
      fa = a - (1 << la);   fa = fa << (F - la);        // fraction, F bits, left aligned
      fb = b - (1 << lb);   fb = fb << (F - lb);
      lg = ((la + lb) << F) + fa + fb;                  // fixed-point log sum
      kk = lg >> F;
      val = (1 << F) | (lg & ((1 << F) - 1));
      y = (kk >= F) ? (val << (kk - F)) : (val >> (F - kk));
    end
  end
endmodule
