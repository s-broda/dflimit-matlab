# Reference values of the general engine (DFLimit.DFDet of DFLimit.jl) for the C test c/dfl_test_det.c: writes
# c/test/det_ref.txt.  Run from the repository root with the Julia package's environment:
#   JLCAP_GB=4 jlcap --project=/path/to/DFLimit.jl c/test/det_ref.jl
# Each line starts with a tag.  'design id spec...' defines a design (trend p | perron M q lam... | ltrend p | fdesign
# name | matrix base n | first base | raw np deg K loc br... C...), the other lines refer to it by id: 'prep' (the sizes
# and the nonzero eigenvalues of A), 'dval' (the determinant at complex (s, t)), 'roots' (log D_V and the reciprocal roots
# at t), 'tstar', and the user-level values (see the end of the file).  Values are printed with 17 digits; every value line
# ends with its tolerance.
using DFLimit, Printf, Logging
const DG = DFLimit.DFDet
const OUT = joinpath(@__DIR__, "det_ref.txt")
io = open(OUT, "w")
fmt(x::Real) = @sprintf("%.17g", x)
fmt(z::Complex) = @sprintf("%.17g %.17g", real(z), imag(z))
line(args...) = println(io, join([a isa AbstractString ? a : a isa Symbol ? String(a) : fmt(a) for a in args], " "))
quiet(f) = with_logger(NullLogger()) do; f(); end

# ---------------------------------------------------------------------------------------------------------------------
# the designs
designs = Dict{Int,Any}()
function defd(id, spec::AbstractString, D)
    designs[id] = D
    println(io, "design $id $spec")
end
defd(1, "trend -1", DG.trend(-1)); defd(2, "trend 0", DG.trend(0)); defd(3, "trend 1", DG.trend(1)); defd(4, "trend 2", DG.trend(2))
defd(5, "trend 3", DG.trend(3)); defd(6, "trend 4", DG.trend(4))
defd(7, "perron A 1 0.3", DG.perron(:A, 0.3)); defd(8, "perron B 1 0.5", DG.perron(:B, 0.5)); defd(9, "perron C 1 0.3", DG.perron(:C, 0.3))
defd(10, "perron C 2 0.3 0.7", DG.perron(:C, [0.3, 0.7])); defd(11, "perron A 3 0.25 0.5 0.75", DG.perron(:A, [0.25, 0.5, 0.75]))
defd(12, "perron C 3 0.25 0.5 0.75", DG.perron(:C, [0.25, 0.5, 0.75]))
defd(13, "ltrend 3", DG.ltrend(3))
defd(14, "fdesign sincos", DG.fdesign(r -> [sin(2pi * r), cos(2pi * r)]))
defd(15, "fdesign perronC04", DG.fdesign(r -> [1.0, r, r > 0.4 ? 1.0 : 0.0, max(r - 0.4, 0.0)], [0.4]))
defd(16, "matrix 9 500", DG.design_from_matrix(DG.regmat(DG.perron(:C, 0.3), 500)))
defd(17, "perron B 1 0.3", DG.perron(:B, 0.3))
defd(18, "first 17", DG.detrend_first(DG.perron(:B, 0.3)))
D19 = DG.Determ([0.0, 0.4, 1.0], [[1.0 1.0; 1.0 0.0], [1.0 0.0; 1.0 0.0]], "1 + r, 1{r < 0.4}")
println(io, "design 19 raw 2 1 2 0 0 0.4 1 ", join(fmt.(vcat(vec(D19.C[1]), vec(D19.C[2]))), " "))
designs[19] = D19
defd(20, "perron A 1 0.4", DG.perron(:A, 0.4))
defd(21, "fdesign expr", DG.fdesign(r -> [exp(r)]))

# ---------------------------------------------------------------------------------------------------------------------
# the prepared designs: N p nj r rp rm, the eigenvalues of A (ascending), and the model of each start: a, kappa
for (id, cb) in [(i, NaN) for i in 1:21]
    P = DG.prep(designs[id])
    line("prep", string(id), cb, string(P.N), string(P.p), string(P.nj), string(P.r), string(P.rp), string(P.rm), P.Dv..., 1e-10)
end
for (id, cb) in ((3, -13.5), (20, -13.5), (2, -13.5), (1, -13.5), (9, -7.0), (5, -10.0))
    P = DG.prep_gls(designs[id], cb)
    line("prep", string(id), cb, string(P.N), string(P.p), string(P.nj), string(P.r), string(P.rp), string(P.rm), P.Dv..., 1e-10)
end
for id in (1, 2, 3, 7, 9, 14, 19, 21)
    Mo = DG.model(DG.prep(designs[id]), :fixed, 1.5)
    line("model", string(id), NaN, "1", 1.5, string(Mo.N), string(Mo.r), string(Mo.rp), string(Mo.rm), Mo.a, Mo.kappa, 1e-12)
    Ms = DG.model(DG.prep(designs[id]), :stationary)
    line("model", string(id), NaN, "2", 0.0, string(Ms.N), string(Ms.r), string(Ms.rp), string(Ms.rm), Ms.a, Ms.kappa, 1e-12)
end

# the determinant: Dval(model, c, s, t) at complex points (zero, fixed: that of the zero start; stationary: extended)
pts = ((0.3 + 2.0im, 0.5 - 7.0im), (-0.8 - 5.0im, -3.0 + 40.0im), (1.1 + 0.0im, 0.9 + 0.0im))
for (id, cb, st, c) in ((1, NaN, 0, -5.0), (2, NaN, 0, 0.0), (3, NaN, 0, -5.0), (4, NaN, 0, 2.0), (6, NaN, 0, -20.0), (9, NaN, 0, -5.0),
                        (10, NaN, 0, 0.0), (12, NaN, 0, -13.5), (13, NaN, 0, -5.0), (14, NaN, 0, 0.0), (15, NaN, 0, -5.0), (16, NaN, 0, -5.0),
                        (18, NaN, 0, -5.0), (19, NaN, 0, 0.0), (3, -13.5, 0, -5.0), (9, -7.0, 0, 0.0), (3, NaN, 2, -5.0), (12, NaN, 2, -8.0),
                        (7, NaN, 2, -0.5), (21, NaN, 0, 0.5))
    D = designs[id]; P = isnan(cb) ? DG.prep(D) : DG.prep_gls(D, cb)
    Mo = DG.model(P, st == 2 ? :stationary : :zero)
    for (s, t) in pts
        d = DG.Dval(Mo, c, s, t)
        line("dval", string(id), cb, string(st), 0.0, c, s, t, d, 1e-11 * abs(d))
    end
end
# the roots at t: log D_V (real part; the imaginary part modulo 2 pi) and the reciprocal roots mu (any order); with a
# fixed start also cf (the exponent's constant)
for (id, cb, st, g0, c) in ((3, NaN, 0, 0.0, -5.0), (5, NaN, 0, 0.0, 0.0), (9, NaN, 0, 0.0, -5.0), (10, NaN, 0, 0.0, 0.0), (12, NaN, 0, 0.0, -13.5),
                            (12, NaN, 2, 0.0, -13.5), (14, NaN, 0, 0.0, 2.0), (18, NaN, 0, 0.0, -5.0), (3, -13.5, 0, 0.0, -5.0), (9, NaN, 1, 1.5, -5.0),
                            (2, NaN, 1, 2.0, 0.5), (19, NaN, 1, 1.0, 0.0), (7, NaN, 2, 0.0, -2.0), (16, NaN, 0, 0.0, -50.0))
    D = designs[id]; P = isnan(cb) ? DG.prep(D) : DG.prep_gls(D, cb)
    Mo = DG.model(P, st == 0 ? :zero : st == 1 ? :fixed : :stationary, g0)
    W = DG.MWork(Mo)
    for t in (0.5 + 0.0im, 2.0 - 7.0im, -3.0 + 40.0im, 300.0 + 800.0im)
        th, ld, mu, fx = DG.mgfparts(Mo, c, t, W)
        cf = fx === nothing ? 0.0im : fx.cf
        line("roots", string(id), cb, string(st), g0, c, t, real(ld), cf, string(length(mu)), mu..., 1e-10)
    end
end
# t_*
for (id, cb, st, c) in ((1, NaN, 0, 0.0), (2, NaN, 0, -5.0), (3, NaN, 0, 0.0), (3, NaN, 0, 3.0), (3, NaN, 0, -150.0), (6, NaN, 0, -10.0),
                        (9, NaN, 0, -5.0), (12, NaN, 0, 0.0), (14, NaN, 0, -5.0), (16, NaN, 0, 0.0), (18, NaN, 0, -5.0), (3, -13.5, 0, -5.0),
                        (2, NaN, 2, -5.0), (12, NaN, 2, -2.0), (5, -10.0, 0, 0.0), (19, NaN, 0, 0.0))
    D = designs[id]; P = isnan(cb) ? DG.prep(D) : DG.prep_gls(D, cb)
    Mo = DG.model(P, st == 2 ? :stationary : :zero)
    ts, m = DG.tstar(Mo, c)
    line("tstar", string(id), cb, string(st), 0.0, c, ts, string(m), 1e-11 * ts)
end
# the double zero of model C at lambda = 1/2 (8 pi^2, multiplicity 2): its t_* by the golden section, to sqrt(eps)
P = DG.prep(DG.perron(:C, 0.5)); ts, m = DG.tstar(P, 0.0)
println(io, "design 22 perron C 1 0.5")
line("tstar", "22", NaN, "0", 0.0, 0.0, ts, string(m), 1e-7 * ts)
close(io)
println("written ", OUT)
