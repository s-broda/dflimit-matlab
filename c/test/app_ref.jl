# reference values of the application functions for the C engine test (c/dfl_test_app.c): c/test/app_ref.txt.  The ADFmax
# statistic (adfmax_cdf, adfmax_pdf), the region {U <= g(V)} (ltratio_region), the linear form U - aV (ltratio_lincdf) and
# the power envelope (power_envelope) of julia/DFExact.jl.
# Columns: "adfmax model stat c x F f", "region type p1 p2 v0 c P", "lincdf model a c cbar k b_1 P_1 ... b_k P_k" (model 1
# case DK, 4 case DG), "vcdf c k b_1 P_1 ... b_k P_k" (P(V <= b), case DK) and "envelope model alpha k c_1 P_1 ... c_k P_k".
using Printf
using DFLimit; const DFExact = DFLimit.DFExact
out = open(ARGS[1], "w")
for (cas, c, x) in (("DU", 0.0, -2.5), ("DU", -10.0, -2.0), ("DU", 0.0, 0.5), ("DU", -30.0, -4.5), ("DT", 0.0, -3.2), ("DT", -10.0, -3.0), ("DT", -5.0, -0.2),
                    ("SU", -10.0, -2.5), ("ST", -10.0, -3.0), ("SU", -5.0, 0.3), ("ST", -25.0, -4.0))
    model = cas[2] == 'U' ? 2 : 3
    @printf(out, "adfmax %d %d %.17g %.17g %.17g %.17g\n", model, cas[1] == 'S' ? 1 : 0, c, x, adfmax_cdf(x, cas, c), adfmax_pdf(x, cas, c))
end
for (c, a, b) in ((0.0, -3.5, -0.2639), (-5.0, -3.5, 1.0))
    @printf(out, "region 1 %.17g %.17g 0 %.17g %.17g\n", a, b, c, ltratio_region((type = "line", a = a, b = b), "DK", c))
end
for (c, x, c0) in ((-7.0, -1.94087, 0.0), (-5.0, -2.0, -5.0))
    @printf(out, "region 2 %.17g %.17g 0 %.17g %.17g\n", x, c0, c, ltratio_region((type = "parabola", x = x, c0 = c0), "DK", c))
end
for (c, k, om) in ((0.0, 3.6658, 10.0), (-10.0, 3.6658, 10.0), (-5.0, 3.6357, 5.0), (0.0, 12.0615, 10.0), (-20.0, 2.3871, 20.0))
    @printf(out, "region 3 %.17g %.17g 0 %.17g %.17g\n", k, om, c, ltratio_region((type = "wap", k = k, omega = om), "DK", c))
end
let (c, a, b, v0) = (0.0, -3.5, -0.2639, 0.05)                      # a boundary that is +Inf below v0
    @printf(out, "region 1 %.17g %.17g %.17g %.17g %.17g\n", a, b, v0, c, ltratio_region((type = "line", a = a, b = b, v0 = v0), "DK", c))   # (P(V < v0) is added)
end
# the linear form U - aV (several b at one a, c, cbar per line) and V
for (cas, a, c, cbar, bs) in (("DK", -3.5, 0.0, -13.5, [-0.2639, 0.5, 2.0]), ("DK", -3.5, -7.0, -13.5, [-0.4, -0.2639, -0.1]), ("DK", 0.0, 0.0, -13.5, [-0.4, 0.0, 1.0]),
                              ("DK", 2.0, -5.0, -13.5, [-0.8, -0.6]), ("DK", 1.0, 2.0, -13.5, [-0.5, 1.5]), ("DK", -10.0, -20.0, -13.5, [-0.3, -0.1]),
                              ("DG", -6.75, 0.0, -13.5, [-0.3, 0.0, 0.5]), ("DG", -6.75, -13.5, -13.5, [-0.3, -0.1]), ("DG", -10.0, -5.0, -20.0, [0.0, 0.5]),
                              ("DG", 1.0, -2.0, -13.5, [-0.8, -0.5]), ("DG", -2.5, 0.0, -5.0, [-0.3]),
                              # where the inversion of Gil-Pelaez is not trusted (small |a|, c > 0): the conditional integrals
                              ("DK", 0.0, -5.0, -13.5, [-0.6, -0.5, -0.4999, -0.3]), ("DK", 0.1, -1.0, -13.5, [-0.55, -0.45, 0.2]),
                              ("DK", -0.01, 0.0, -13.5, [-0.5, -0.49, 0.5]), ("DK", 1.0, 5.0, -13.5, [-0.5, 3.0]),
                              ("DG", 0.0, 0.0, -13.5, [-0.6, -0.5, -0.499, -0.45]), ("DG", -0.3, -5.0, -13.5, [-0.48, -0.4]),
                              ("DG", 0.05, 2.0, -13.5, [-0.5, -0.3]), ("DG", -0.001, 0.0, -1.0, [-0.499, -0.49]))
    P = ltratio_lincdf(bs, a, cas, c; cbar = cbar)
    @printf(out, "lincdf %d %.17g %.17g %.17g %d", cas == "DK" ? 1 : 4, a, c, cbar, length(bs))
    for (b, p) in zip(bs, P); @printf(out, " %.17g %.17g", b, p); end
    @printf(out, "\n")
end
for (c, bs) in ((0.0, [0.1, 0.5, 2.0]), (-10.0, [0.05]), (3.0, [1.0]), (-200.0, [0.0025, 0.003]), (5.0, [2.0]))
    P = ltratio_lincdf(bs, NaN, "DK", c)
    @printf(out, "vcdf %.17g %d", c, length(bs))
    for (b, p) in zip(bs, P); @printf(out, " %.17g %.17g", b, p); end
    @printf(out, "\n")
end
# the power envelope (several c per line)
for (cas, alpha, cs) in (("DK", 0.05, [-1.0, -7.0, -20.0, 0.0]), ("DK", 0.01, [-10.0]), ("DG", 0.05, [-5.0, -13.5, 1.0]), ("DG", 0.1, [-25.0]),
                         ("DK", 0.05, [-0.01, -1e-4]), ("DG", 0.05, [-1e-3, -0.1]))   # (near c = 0: the conditional integrals)
    P = power_envelope(cs, cas; alpha = alpha)
    @printf(out, "envelope %d %.17g %d", cas == "DK" ? 1 : 4, alpha, length(cs))
    for (c, p) in zip(cs, P); @printf(out, " %.17g %.17g", c, p); end
    @printf(out, "\n")
end
close(out)
