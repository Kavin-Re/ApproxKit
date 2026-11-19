// ApproxKit v2: DRUM squarer, 2 <= K <= N.  Matches approxkit_v2.divsq.DRUMSquarer.
module sq_drum #(parameter N = 8, parameter K = 4) (input [N-1:0] a, output reg [2*N-1:0] y);
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin msb = 0; for (m = 0; m < N; m = m + 1) if (x[m]) msb = m; end
  endfunction
  integer h;
  reg [K-1:0] s;
  reg [2*K-1:0] pk;
  always @* begin
    h = msb(a) - K + 1; if (h < 0) h = 0;
    s = a >> h; if (h > 0) s[0] = 1'b1;
    pk = s * s;
    y = pk;
    y = y << (2*h);
  end
endmodule
