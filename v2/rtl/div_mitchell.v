// ApproxKit v2: Mitchell logarithmic divider (floor), N >= 2.  Matches approxkit_v2.divsq.MitchellDiv.
module div_mitchell #(parameter N = 8) (input [N-1:0] a, b, output reg [N-1:0] y);
  localparam F = N - 1;
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin msb = 0; for (m = 0; m < N; m = m + 1) if (x[m]) msb = m; end
  endfunction
  integer la, lb, fa, fb, d, kk;
  reg [2*N:0] val;
  always @* begin
    la = msb(a); lb = msb(b);
    fa = (a - (1 << la)) << (F - la);
    fb = (b - (1 << lb)) << (F - lb);
    d  = ((la - lb) << F) + fa - fb;
    kk = d >>> F;                                   // arithmetic shift of a signed integer: floor
    val = (1 << F) | (d & ((1 << F) - 1));
    if (b == 0)      y = {N{1'b1}};
    else if (a == 0) y = 0;
    else             y = (kk >= F) ? (val << (kk - F)) : (val >> (F - kk));
  end
endmodule
