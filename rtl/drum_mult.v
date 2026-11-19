// ApproxKit: DRUM-K dynamic-range unbiased multiplier; 2 <= K < N
// Bit-accurate to approxkit.DRUM (see README). Verilog-2001, no vendor primitives.
module drum_mult #(parameter N = 8, parameter K = 4) (input [N-1:0] a, b, output reg [2*N-1:0] y);
  function integer msb;               // leading-one index, -1 if x == 0
    input [N-1:0] x;
    integer m;
    begin
      msb = -1;
      for (m = 0; m < N; m = m + 1) if (x[m]) msb = m;
    end
  endfunction
  integer la, lb, sa, sb;
  reg [K-1:0]   ga, gb;
  reg [2*K-1:0] pk;
  always @* begin
    la = msb(a);  lb = msb(b);
    sa = (la >= K) ? la - K + 1 : 0;
    sb = (lb >= K) ? lb - K + 1 : 0;
    ga = a >> sa;  if (sa > 0) ga[0] = 1'b1;
    gb = b >> sb;  if (sb > 0) gb[0] = 1'b1;
    pk = ga * gb;
    y  = pk;
    y  = y << (sa + sb);
  end
endmodule
