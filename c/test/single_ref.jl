# reference values of the single integral at c0 = 0 (cdf and density) for the C engine test: c/test/single_ref.txt
using Printf
using DFLimit; const DFExact = DFLimit.DFExact; const D = DFExact
out = open(ARGS[1], "w")
for model in 1:4, (c, st) in ((0.0, false), (-5.0, false), (-20.0, false), (0.5, false), (3.0, false), (10.0, false), (-5.0, true), (-10.0, true))
    model == 4 && (st || c > 2) && continue                 # (model 4, the GLS-detrended trend case: zero start; c <= 2 here)
    xs = model == 1 ? [-3.0, -1.94, -0.72, 0.0, 0.5, 2.0, 8.0] : [-4.0, -2.86, -1.5, -0.5, 0.3, 2.0, 8.0]
    for x in xs
        c > 2 && x < 0 && continue
        F = D.cdf(x, model; c, start = st ? :stationary : :zero); f = D.density(x, model; c, start = st ? :stationary : :zero)
        @printf(out, "%d %.17g %d %.17g %.17g %.17g\n", model, c, st ? 1 : 0, x, F, f)
    end
end
# far from the unit root: the points of test/engine/largec_ref.txt of DFLimit.jl with c0 = 0 (the values stored there, from the line integral)
for l in eachline(joinpath(pkgdir(DFLimit), "test", "engine", "largec_ref.txt"))
    p = split(l); (p[1] == "val" && p[4] == "0") || continue
    @printf(out, "%d %s %d %s %s %s\n", findfirst(==(p[2][2]), "KUT"), p[3], p[2][1] == 'S' ? 1 : 0, p[6], p[7], p[8])
end
# case DG (model 4) far from the unit root, at the median (November 2026: stable to 2e-12 under node doubling, the adaptive
# rules and a second line at c = -100, -480, the density integrating to one to 1e-13)
for (c, x) in ((-100.0, -7.128497349301127), (-480.0, -15.446799628855823))
    @printf(out, "%d %.17g %d %.17g %.17g %.17g\n", 4, c, 0, x, D.cdf(x, 4; c), D.density(x, 4; c))
end
close(out)
