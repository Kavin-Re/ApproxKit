// ApproxKit v2: Mitchell reciprocal  R = floor((3*2^k - x) * 2^(M-2k-1)),  x=0 -> all ones.
module recip_mitchell #(parameter N = 8, parameter M = 8) (input [N-1:0] a, output reg [M:0] y);
  function integer msb;
    input [N-1:0] x;
    integer m;
    begin msb = 0; for (m = 0; m < N; m = m + 1) if (x[m]) msb = m; end
  endfunction
  integer k, s;
  reg [N+M+3:0] v, r;
  always @* begin
    k = msb(a);
    v = (3 << k) - a;
    s = M - 2*k - 1;
    r = (s >= 0) ? (v << s) : (v >> (-s));
    y = (a == 0) ? {(M+1){1'b1}} : r[M:0];
  end
endmodule
