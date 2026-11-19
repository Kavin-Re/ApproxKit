// ApproxKit v2: Baugh-Wooley signed multiplier with approximate partial products.
// Matches approxkit_v2.signed.BaughWooley.  K = approximation column bound.
//   MODE 0 'trunc': pp bits in columns < K dropped (BW constants kept); COMP=1 adds the
//                   rounded expected value of the dropped bits.  0 <= K <= 2N
//   MODE 1 'orlow': columns < K reduced with OR (no carries), columns >= K exact.  K <= N
// K = 0 is the exact signed multiplier (explicit BW matrix + chain of adds).
module bw_mult #(parameter N = 8, parameter K = 0, parameter MODE = 0, parameter COMP = 0)
               (input [N-1:0] a, b, output reg [2*N-1:0] y);
  localparam [2*N-1:0] KEEP = {(2*N){1'b1}} << K;
  localparam [2*N-1:0] LOWM = ~KEEP;
  function integer compval;               // (4*E[dropped] + 2) / 4, integer arithmetic
    input integer dummy;
    integer ii, jj, e4;
    begin
      e4 = 0;
      for (ii = 0; ii < N; ii = ii + 1)
        for (jj = 0; jj < N; jj = jj + 1)
          if (ii + jj < K) e4 = e4 + (((ii == N-1) != (jj == N-1)) ? 3 : 1) * (1 << (ii + jj));
      compval = (e4 + 2) / 4;
    end
  endfunction
  localparam integer CV = (COMP != 0 && K != 0 && MODE == 0) ? compval(0) : 0;
  reg [2*N-1:0] row, acc, orl;
  reg t;
  integer i, j;
  always @* begin
    acc = 0; orl = 0;
    for (j = 0; j < N; j = j + 1) begin
      row = 0;
      for (i = 0; i < N; i = i + 1) begin
        t = a[i] & b[j];
        if ((i == N-1) != (j == N-1)) t = ~t;
        row[i + j] = t;
      end
      acc = acc + (row & KEEP);
      orl = orl | (row & LOWM);
    end
    acc = acc + (1 << N) + (1 << (2*N-1)) + CV;
    y = (MODE == 1) ? (acc | orl) : acc;
  end
endmodule
