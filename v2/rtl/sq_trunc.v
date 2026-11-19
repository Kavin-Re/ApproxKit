// ApproxKit v2: folded squarer matrix, pp bits in columns < K dropped; COMP adds round(E[dropped]).
// K = 0 is the exact squarer (explicit folded matrix).  Matches approxkit_v2.divsq.TruncSquarer.
module sq_trunc #(parameter N = 8, parameter K = 0, parameter COMP = 0) (input [N-1:0] a, output reg [2*N-1:0] y);
  function integer compval;
    input integer dummy;
    integer i, j, e4;
    begin
      e4 = 0;
      for (i = 0; i < N; i = i + 1) begin
        if (2*i < K) e4 = e4 + (2 << (2*i));
        for (j = i + 1; j < N; j = j + 1) if (i + j + 1 < K) e4 = e4 + (1 << (i + j + 1));
      end
      compval = (e4 + 2) / 4;
    end
  endfunction
  localparam integer CV = (COMP != 0 && K != 0) ? compval(0) : 0;
  integer i, j;
  always @* begin
    y = CV;
    for (i = 0; i < N; i = i + 1) begin
      if (2*i >= K) y = y + (a[i] << (2*i));
      for (j = i + 1; j < N; j = j + 1)
        if (i + j + 1 >= K) y = y + ((a[i] & a[j]) << (i + j + 1));
    end
  end
endmodule
