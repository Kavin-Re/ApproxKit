// ApproxKit v2: Mitchell squarer (doubled log, linear antilog, floor).  Matches MitchellSquarer.
module sq_mitchell #(parameter N = 8) (input [N-1:0] a, output reg [2*N-1:0] y);
  localparam F = N - 1;
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin msb = 0; for (m = 0; m < N; m = m + 1) if (x[m]) msb = m; end
  endfunction
  integer l, f, L2, kk;
  reg [2*N:0] val;
  always @* begin
    l  = msb(a);
    f  = (a - (1 << l)) << (F - l);
    L2 = ((l << F) + f) << 1;
    kk = L2 >> F;
    val = (1 << F) | (L2 & ((1 << F) - 1));
    y = (a == 0) ? 0 : ((kk >= F) ? (val << (kk - F)) : (val >> (F - kk)));
  end
endmodule
