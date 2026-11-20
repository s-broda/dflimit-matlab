# reference values of the single integral with a fixed start g0 != 0 (the noncentral kernel): cdf and density, for the C
# engine test c/test/fixed_ref.txt (columns: model c c0 g0 stat x F f)
using Printf
using DFLimit; const DFExact = DFLimit.DFExact; const D = DFExact
out = open(ARGS[1], "w")
P = ((1, -5.0, 0.0, 1.0, [-2.0, 0.5]), (1, 0.0, 0.0, 1.5, [-1.0]), (1, 2.0, 0.0, 1.0, [1.0, -1.5]), (1, -5.0, -5.0, 1.0, [-2.0]), (1, 2.0, 2.0, 1.0, [-1.0, -2.5]),
     (2, -5.0, 0.0, 1.0, [-2.86, -1.0, 0.5]), (2, -2.0, 0.0, 1.5, [-2.0]), (2, 2.0, 0.0, 1.0, [0.0, -3.0]), (2, -5.0, -5.0, 1.0, [-2.0]), (2, 1.0, 1.0, 1.0, [0.0]),
     (3, -5.0, 0.0, 2.0, [-3.0, 0.0]), (3, 2.0, 0.0, 1.0, [-2.0]), (3, -2.5, -5.0, 1.0, [-1.1]), (3, 0.0, 2.0, 0.5, [-2.3, -1.9]),
     (2, -20.0, 0.0, 1.0, [-4.0, -1.0]), (2, -5.0, 0.0, 3.0, [-3.0]), (3, -10.0, 0.0, 3.0, [-3.0]), (2, 5.0, 0.0, 1.0, [-3.0, 2.0]), (1, 5.0, 0.0, 1.0, [0.0]),
     (2, 3.0, 0.0, 1.0, [0.0]), (2, -20.0, 0.0, 1.5, [-4.0]),
     # November 2026: (c g0)^2/2 > 709 in case DU (the line integral returned NaN), and case DK where P[J(1)^2 <= 1 + g0^2] is of
     # the size e^{-lam} (the sum of the two erf cancelled: F off by 0.33 at c = 0, g0 = 3)
     (2, -18.0, 0.0, 3.0, [-9.0806812422194216]), (2, -26.0, 0.0, 1.5, [-6.4421804207241742]), (1, 0.0, 0.0, 3.0, [-0.16692]), (1, 2.0, 0.0, 3.0, [22.09087]))
for (model, c, c0, g0, xs) in P
    for x in xs
        F = D.cdf(x, model; c, c0, g0); f = D.density(x, model; c, c0, g0)
        @printf(out, "%d %.17g %.17g %.17g %d %.17g %.17g %.17g\n", model, c, c0, g0, 0, x, F, f)
    end
end
close(out)
