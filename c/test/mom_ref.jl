# reference moments for the C engine test: c/test/mom_ref.txt (model c c0 g0 stat m1..m4)
using Printf
using DFLimit; const DFExact = DFLimit.DFExact; const D = DFExact
out = open(ARGS[1], "w")
for (model, c, c0, g0, st) in ((1, 0.0, 0.0, 0.0, false), (2, 0.0, 0.0, 0.0, false), (3, 0.0, 0.0, 0.0, false), (2, -5.0, 0.0, 0.0, false), (3, -10.0, 0.0, 0.0, false),
                                (1, 5.0, 0.0, 0.0, false), (2, 5.0, 0.0, 0.0, false), (2, 10.0, 10.0, 0.0, false), (3, -5.0, -5.0, 0.0, false), (1, 2.0, 2.0, 0.0, false),
                                (2, -5.0, 0.0, 0.0, true), (1, -2.5, 0.0, 0.0, true), (3, -10.0, 1.0, 0.0, true),
                                (1, -5.0, 0.0, 1.5, false), (2, -5.0, 0.0, 1.5, false), (3, -5.0, 0.0, 0.8, false), (2, 2.0, 1.0, 1.0, false))
    m = D.moments(model; c, c0, g0, start = st ? :stationary : :zero, pmax = 4)
    @printf(out, "%d %.17g %.17g %.17g %d %.17g %.17g %.17g %.17g\n", model, c, c0, g0, st ? 1 : 0, m...)
end
# far from the unit root: the moments of test/engine/largec_ref.txt of DFLimit.jl (c = -60 ... -480, the six cases) with c0 = 0 and, with a zero
# start, c0 = c down to c = -200 (with c0 = c the third and fourth moments converge slowly in the step in log s, to 1e-10
# with a zero start and 1e-9 to 1e-6 with a stationary one, and the C and Julia values differ at that level)
for l in eachline(joinpath(pkgdir(DFLimit), "test", "engine", "largec_ref.txt"))
    p = split(l); (p[1] == "mom" && (p[4] == "0" || (p[2][1] == 'D' && p[3] != "-480"))) || continue
    @printf(out, "%d %s %s 0 %d %s %s %s %s\n", findfirst(==(p[2][2]), "KUT"), p[3], p[4], p[2][1] == 'S' ? 1 : 0, p[5], p[6], p[7], p[8])
end
close(out)
