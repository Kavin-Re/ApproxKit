// ApproxKit v2: sign-magnitude wrapper around a v1 UNSIGNED multiplier.
// Matches approxkit_v2.signed.SignMag.  UNIT: 0 exact('*'), 1 trunc_mult(K=P1),
// 2 drum_mult(K=P1), 3 mitchell_mult, 4 booth_mult(K=P1), 5 bam_mult(VBL=P1,HBL=P2).
// SAT=0: unit width N.  SAT=1: unit width N-1, |x| saturated to 2^(N-1)-1.
// Needs v1 rtl/{trunc_mult,bam_mult,drum_mult,mitchell_mult,booth_mult}.v on the compile line.
module smag_mult #(parameter N = 8, parameter UNIT = 0, parameter P1 = 0, parameter P2 = 0,
                   parameter SAT = 0) (input [N-1:0] a, b, output [2*N-1:0] y);
  localparam W = SAT ? N - 1 : N;
  wire sa = a[N-1], sb = b[N-1];
  wire [N-1:0] ma = sa ? (~a + 1'b1) : a;       // |a| <= 2^(N-1): fits N bits
  wire [N-1:0] mb = sb ? (~b + 1'b1) : b;
  wire [W-1:0] ua, ub;
  generate if (SAT) begin : g_sat
    assign ua = ma[N-1] ? {W{1'b1}} : ma[W-1:0];
    assign ub = mb[N-1] ? {W{1'b1}} : mb[W-1:0];
  end else begin : g_nosat
    assign ua = ma;
    assign ub = mb;
  end endgenerate
  wire [2*W-1:0] p;
  generate
    if (UNIT == 0)      begin : u0 assign p = ua * ub; end
    else if (UNIT == 1) begin : u1 trunc_mult   #(.N(W), .K(P1))             u (.a(ua), .b(ub), .y(p)); end
    else if (UNIT == 2) begin : u2 drum_mult    #(.N(W), .K(P1))             u (.a(ua), .b(ub), .y(p)); end
    else if (UNIT == 3) begin : u3 mitchell_mult #(.N(W))                    u (.a(ua), .b(ub), .y(p)); end
    else if (UNIT == 4) begin : u4 booth_mult   #(.N(W), .K(P1))             u (.a(ua), .b(ub), .y(p)); end
    else                begin : u5 bam_mult     #(.N(W), .VBL(P1), .HBL(P2)) u (.a(ua), .b(ub), .y(p)); end
  endgenerate
  wire [2*N-1:0] pz  = p;                                   // zero-extended
  wire [2*N-1:0] mag = pz[2*N-1] ? {1'b0, {(2*N-1){1'b1}}} : pz;   // clamp to 2^(2N-1)-1
  assign y = (sa ^ sb) ? (~mag + 1'b1) : mag;
endmodule
