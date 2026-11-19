// ApproxKit v2: signed DRUM on two's-complement operands, 2 <= K <= N (K counts the sign bit).
// Matches approxkit_v2.signed.SignedDRUM.
module sdrum_mult #(parameter N = 8, parameter K = 4) (input [N-1:0] a, b, output reg [2*N-1:0] y);
  function integer lead;                  // highest bit differing from the sign bit, -1 if none
    input [N-1:0] x;
    reg [N-1:0] t;
    integer m;
    begin
      t = x[N-1] ? ~x : x;
      lead = -1;
      for (m = 0; m < N; m = m + 1) if (t[m]) lead = m;
    end
  endfunction
  integer sa, sb;
  reg signed [N-1:0] fa, fb;
  reg signed [K-1:0] ga, gb;
  reg signed [2*K-1:0] pk;
  reg signed [2*N-1:0] pe;
  always @* begin
    sa = lead(a) + 2 - K; if (sa < 0) sa = 0;
    sb = lead(b) + 2 - K; if (sb < 0) sb = 0;
    fa = $signed(a) >>> sa;  if (sa > 0) fa[0] = 1'b1;
    fb = $signed(b) >>> sb;  if (sb > 0) fb[0] = 1'b1;
    ga = fa[K-1:0];  gb = fb[K-1:0];
    pk = ga * gb;
    pe = pk;
    y  = pe <<< (sa + sb);
  end
endmodule
