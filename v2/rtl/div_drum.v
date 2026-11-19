// ApproxKit v2: segment (DRUM-style) divider, 2 <= K <= N.  Matches approxkit_v2.divsq.DRUMDiv.
module div_drum #(parameter N = 8, parameter K = 4) (input [N-1:0] a, b, output reg [N-1:0] y);
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin msb = 0; for (m = 0; m < N; m = m + 1) if (x[m]) msb = m; end
  endfunction
  integer ha, hb, sh;
  reg [K-1:0] sa, sb;
  reg [2*N+K:0] qf, q;
  always @* begin
    ha = msb(a) - K + 1; if (ha < 0) ha = 0;
    hb = msb(b) - K + 1; if (hb < 0) hb = 0;
    sa = a >> ha; if (ha > 0) sa[0] = 1'b1;
    sb = b >> hb; if (hb > 0) sb[0] = 1'b1;
    qf = ({{(2*N+1){1'b0}}, sa} << K) / ((sb == 0) ? 1 : sb);
    sh = ha - hb - K;
    q  = (sh >= 0) ? (qf << sh) : (qf >> (-sh));
    if (b == 0)      y = {N{1'b1}};
    else if (a == 0) y = 0;
    else             y = (q > {N{1'b1}}) ? {N{1'b1}} : q[N-1:0];
  end
endmodule
