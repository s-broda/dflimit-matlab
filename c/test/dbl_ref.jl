# reference values of the double integral (cdf and density) for the C engine test: c/test/dbl_ref.txt (the checked
# computation, verified = true: dfl_test_dbl compares the checked computation with them and one evaluation, by the rule
# of one evaluation since November 2026, at a looser tolerance)
using Printf
using DFLimit; const DFExact = DFLimit.DFExact; const D = DFExact
out = open(ARGS[1], "w")
P = ((1, 0.0, 0.0, 0.0, false, [-2.5, -0.5, 1.5]), (2, 0.0, 0.0, 0.0, false, [-3.4, -1.6, 0.2]), (3, -5.0, 0.0, 0.0, false, [-3.5, -2.5]),
     (1, 2.0, 2.0, 0.0, false, [-1.5, 0.0, 1.5]), (2, 2.0, 1.0, 0.0, false, [-3.0, -0.5, 1.3]), (3, -2.5, -5.0, 0.0, false, [-3.6, -1.1]), (2, 5.0, -5.0, 0.0, false, [57.11]),
     (1, -5.0, 0.0, 0.0, true, [-2.5, -1.0]), (2, -5.0, 0.0, 0.0, true, [-3.0, -1.2]), (3, -10.0, 0.0, 0.0, true, [-3.5, -2.2]),
     (1, -5.0, 0.0, 1.5, false, [-3.0, -2.0, 0.0]), (2, -5.0, 0.0, 1.5, false, [-3.0, -2.0]), (3, -5.0, 0.0, 1.5, false, [-3.5, -2.5]), (2, 2.0, 0.0, 1.0, false, [-2.0, 0.0]), (3, -10.0, -3.0, 0.8, false, [-3.0, -1.0]), (2, -2.5, 1.0, 0.0, true, [-2.0, 0.0]))
for (model, c, c0, g0, st, xs) in P
    F, _ = D.tr_limdist(xs, model, c; c0, g0, start = st ? :stationary : :zero, verified = true); f, _ = D.tr_limdist(xs, model, c; c0, g0, start = st ? :stationary : :zero, pdf = true, verified = true)
    for (k, x) in enumerate(xs); @printf(out, "%d %.17g %.17g %.17g %d %.17g %.17g %.17g\n", model, c, c0, g0, st ? 1 : 0, x, F[k], f[k]); end
end
# far from the unit root: cases DU, DT with a zero start at c0 = 0 (the lower end of the outer integral for c < -20)
for (model, c, xs) in ((2, -200.0, [-10.84, -10.08, -9.07]), (2, -480.0, [-16.54, -15.54, -14.54]), (3, -60.0, [-6.73, -5.71, -4.67]), (3, -200.0, [-10.89, -10.13, -9.12]))
    F, _ = D.tr_limdist(xs, model, c; verified = true); f, _ = D.tr_limdist(xs, model, c; pdf = true, verified = true)
    for (k, x) in enumerate(xs); @printf(out, "%d %.17g %.17g %.17g %d %.17g %.17g %.17g\n", model, c, 0.0, 0.0, 0, x, F[k], f[k]); end
end
close(out)
