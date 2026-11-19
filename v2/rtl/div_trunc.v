// ApproxKit v2: restoring divider with the last T iterations removed (T = 0 exact).
// COMP=1 sets quotient bit T-1.  b = 0 -> all ones.  Matches approxkit_v2.divsq.TruncRestoringDiv.
module div_trunc #(parameter N = 8, parameter T = 0, parameter COMP = 0) (input [N-1:0] a, b, output reg [N-1:0] y);
  reg [N:0] r;
  reg [N-1:0] q;
  integer i;
  always @* begin
    r = 0; q = 0;
    for (i = N - 1; i >= T; i = i - 1) begin
      r = {r[N-1:0], a[i]};
      if (r >= {1'b0, b}) begin r = r - {1'b0, b}; q[i] = 1'b1; end
    end
    if (COMP != 0 && T > 0) q[T-1] = 1'b1;
    y = (b == 0) ? {N{1'b1}} : q;
  end
endmodule
