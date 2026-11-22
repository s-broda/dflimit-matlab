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

# ---------------------------------------------------------------------------------------------------------------------
# the outer paths: tau t_* a b h phi K line minpos maxneg, and the roots and log P at three nodes
startsym(st) = st == 0 ? :zero : st == 1 ? :fixed : :stationary
function getpath(id, cb, st, g0, c)
    D = designs[id]
    DG.cpath_cached(DG.model(DG.prep_cached(D, isnan(cb) ? nothing : cb), startsym(st), g0), c)
end
for (id, cb, st, g0, c) in ((1, NaN, 0, 0.0, -5.0), (3, NaN, 0, 0.0, 0.0), (3, NaN, 0, 0.0, -13.5), (9, NaN, 0, 0.0, -5.0), (10, NaN, 0, 0.0, -20.0),
                            (12, NaN, 0, 0.0, 0.0), (3, NaN, 0, 0.0, 2.0), (2, NaN, 0, 0.0, -150.0), (3, NaN, 1, 2.0, -5.0), (1, NaN, 1, 1.5, 0.0),
                            (9, NaN, 1, 1.0, -5.0), (2, NaN, 2, 0.0, -5.0), (12, NaN, 2, 0.0, -5.0), (3, -13.5, 0, 0.0, -5.0), (3, -5.0, 0, 0.0, 0.0),
                            (18, NaN, 0, 0.0, -5.0), (14, NaN, 0, 0.0, -5.0), (16, NaN, 0, 0.0, 0.0), (19, NaN, 0, 0.0, 0.0), (5, NaN, 0, 0.0, -50.0))
    CP = quiet(() -> getpath(id, cb, st, g0, c))
    K = length(CP.t); ks = unique([1, cld(K, 2), K])
    line("path", string(id), cb, string(st), g0, c, CP.tau, CP.tstar, CP.a, CP.b, CP.h, CP.phi, string(K), CP.minpos, CP.maxneg,
         string(length(ks)), [x for k in ks for x in (Float64(k - 1), real(CP.lP[k]), imag(CP.lP[k]), (reim(z)[i] for z in CP.R[:, k] for i in 1:2)...)]..., 1e-9)
end

DG.empty_caches!(); GC.gc()
# (GLS detrending at cbar = -2 with c = 0 is left out: its path exhausts the memory cap of 4 GB in the Julia engine, as the
# README of DFLimit.jl notes for the envelope near the unit root)
# ---------------------------------------------------------------------------------------------------------------------
# the user-level values: 'val func id cbar start g0 c c0 verified par n x... values... tol'; the points of a line are
# evaluated together (as one call).  The tolerance is 1e-9 for distribution functions and densities, 1e-7 for quantiles,
# and where the Julia engine warns that one evaluation is brittle (its rounding estimate or its known regimes) the
# tolerance is 1e-6 and the line says so ('brittle' after the tolerance).
struct WarnCount <: Logging.AbstractLogger; n::Base.RefValue{Int}; end
Logging.min_enabled_level(::WarnCount) = Logging.Warn
Logging.shouldlog(::WarnCount, args...) = true
Logging.handle_message(w::WarnCount, args...; kw...) = (w.n[] += 1; nothing)
function warned(f)
    w = WarnCount(Ref(0))
    v = with_logger(f, w)
    v, w.n[] > 0
end
function val(func, id, cb, st, g0, c, c0, ver, par, xs, f; tol = 1e-9)
    v, w = warned(f)
    v = vcat(v)
    DG.empty_caches!()
    t = w ? max(tol, 1e-6) : tol
    println(io, join(["val", func, string(id), fmt(cb), string(st), fmt(g0), fmt(c), fmt(c0), string(ver), fmt(par), string(length(xs)),
                      fmt.(xs)..., string(length(v)), fmt.(v)..., fmt(t)], " ") * (w ? " brittle" : ""))
end
kw(cb, st, g0, c0, ver) = (c0 = c0, start = startsym(st), gamma0 = g0, gls = isnan(cb) ? nothing : cb, verified = ver)
kw4(cb, st, g0, c0) = (c0 = c0, start = startsym(st), gamma0 = g0, gls = isnan(cb) ? nothing : cb)
kwl(cb, st, g0) = (start = startsym(st), gamma0 = g0, gls = isnan(cb) ? nothing : cb)
function tcell(id, cb, st, g0, c, c0, xs; ver = false, pdf = true)
    D = designs[id]; k = kw(cb, st, g0, c0, ver)
    val("cdf", id, cb, st, g0, c, c0, ver ? 1 : 0, 0.0, xs, () -> DG.dist_cdf(D, xs, c; k...))
    pdf && val("pdf", id, cb, st, g0, c, c0, ver ? 1 : 0, 0.0, xs, () -> DG.dist_pdf(D, xs, c; k...))
end
# every constructor, breaks, degrees 0 to 4, c in {0, -5, -13.5, -20, -50, 2}
xs5 = [-4.0, -3.0, -2.2, -1.5, -0.5, 0.8]
for id in (1, 2, 3, 4, 5, 6), c in (0.0, -5.0)
    tcell(id, NaN, 0, 0.0, c, 0.0, xs5)
end
for (id, c) in ((3, -13.5), (3, -20.0), (3, -50.0), (3, 2.0), (2, -13.5), (1, 2.0), (2, -50.0), (6, -20.0), (4, 2.0))
    tcell(id, NaN, 0, 0.0, c, 0.0, xs5)
end
# (design 19, d = (1 + r, 1{r < 0.4}) without a proper intercept, is left out of the values: its outer path is the vertical
# line, and its distribution function took more than the memory cap of 4 GB in the Julia engine; it is in the tests of
# the layers and of the paths above)
for id in (7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21), c in (0.0, -13.5)
    tcell(id, NaN, 0, 0.0, c, 0.0, [-5.0, -3.8, -3.0, -2.0, -0.5])
end
# the far tails and the bulk at one point each
for (id, c) in ((1, 0.0), (3, 0.0), (9, -5.0))
    tcell(id, NaN, 0, 0.0, c, 0.0, [-8.0]); tcell(id, NaN, 0, 0.0, c, 0.0, [-6.5]); tcell(id, NaN, 0, 0.0, c, 0.0, [3.5]); tcell(id, NaN, 0, 0.0, c, 0.0, [5.0])
end
# the starts
for (id, g0, c) in ((1, 1.5, -5.0), (2, 1.0, -5.0), (3, 2.0, -5.0), (3, 1.0, 0.0), (1, 1.0, 2.0), (9, 1.0, -5.0), (7, 1.5, -13.5), (2, 1.5, -20.0), (14, 1.0, -5.0))
    tcell(id, NaN, 1, g0, c, 0.0, [-4.0, -2.5, -1.5, 0.0])
end
for (id, c) in ((1, -5.0), (2, -5.0), (3, -10.0), (2, -0.5), (9, -13.5), (12, -5.0), (3, -50.0), (13, -5.0))
    tcell(id, NaN, 2, 0.0, c, 0.0, [-4.0, -2.5, -1.5, 0.0])
end
# GLS detrending
for (id, cb, c) in ((3, -13.5, 0.0), (3, -13.5, -5.0), (2, -7.0, -5.0), (20, -13.5, -5.0), (9, -13.5, -5.0), (10, -13.5, 0.0), (5, -10.0, -13.5), (3, -13.5, 2.0), (3, -2.0, -5.0))
    tcell(id, cb, 0, 0.0, c, 0.0, [-3.5, -2.5, -1.8, -1.0, 0.3])
end
# a local null value
for (id, cb, st, g0, c, c0) in ((3, NaN, 0, 0.0, -5.0, 3.0), (3, NaN, 0, 0.0, -5.0, -3.0), (2, NaN, 0, 0.0, -5.0, -5.0), (1, NaN, 0, 0.0, -10.0, 2.0),
                                (9, NaN, 0, 0.0, -13.5, -13.5), (3, NaN, 2, 0.0, -5.0, -5.0), (1, NaN, 1, 1.0, -10.0, -2.0), (3, -13.5, 0, 0.0, -5.0, 3.0),
                                (10, NaN, 0, 0.0, -20.0, -20.0), (2, NaN, 0, 0.0, 0.0, 2.0))
    tcell(id, cb, st, g0, c, c0, [-4.5, -3.0, -1.8, -0.5, 1.5])
end
# the tangency of case DK with c0 > 0 (x* = -sqrt(2 c0)): just above and below it, and at it (an infinite density)
for x in (-2.0 + 1e-3, -2.0 - 1e-3, -2.0)
    tcell(1, NaN, 0, 0.0, -5.0, 2.0, [x])
end
# the statistic that detrends first; c > 0 on the vertical line (Filon's rule with c0 != 0)
tcell(18, NaN, 0, 0.0, 0.0, 0.0, [-4.5, -3.8, -3.0, -2.0])
tcell(1, NaN, 0, 0.0, 5.0, -3.0, [-1.0, 1.0, 4.0]; pdf = true)
tcell(2, NaN, 0, 0.0, 4.0, 0.0, [-1.0, 1.0, 3.0])
# the checked computation
tcell(3, NaN, 0, 0.0, -5.0, 0.0, [-3.0, -1.5]; ver = true)
tcell(2, NaN, 1, 1.5, -5.0, 0.0, [-3.0, -1.5]; ver = true)
# quantiles at several levels
for (id, cb, st, g0, c, c0) in ((1, NaN, 0, 0.0, 0.0, 0.0), (3, NaN, 0, 0.0, 0.0, 0.0), (9, NaN, 0, 0.0, -5.0, 0.0), (12, NaN, 0, 0.0, 0.0, 0.0), (3, -13.5, 0, 0.0, 0.0, 0.0),
                                (2, NaN, 1, 1.0, -5.0, 0.0), (2, NaN, 2, 0.0, -5.0, 0.0), (3, NaN, 0, 0.0, -5.0, -5.0), (18, NaN, 0, 0.0, 0.0, 0.0), (14, NaN, 0, 0.0, -13.5, 0.0))
    D = designs[id]
    for p in (0.01, 0.05, 0.1, 0.5, 0.95)
        val("quant", id, cb, st, g0, c, c0, 0, -2.0, [p], () -> DG.dist_quantile(D, p, c; x0 = -2.0, kw4(cb, st, g0, c0)...); tol = 1e-7)
    end
end
# the coefficient statistic
for (id, cb, st, g0, c, c0) in ((3, NaN, 0, 0.0, -5.0, 0.0), (2, NaN, 0, 0.0, 0.0, 0.0), (9, NaN, 0, 0.0, -5.0, 0.0), (3, -13.5, 0, 0.0, -5.0, 0.0), (2, NaN, 2, 0.0, -5.0, -5.0), (1, NaN, 1, 1.0, -5.0, 0.0))
    D = designs[id]; ys = [-25.0, -12.0, -6.0, -1.0]; k = kw(cb, st, g0, c0, false)
    val("coefcdf", id, cb, st, g0, c, c0, 0, 0.0, ys, () -> DG.coef_cdf(D, ys, c; k...))
    val("coefpdf", id, cb, st, g0, c, c0, 0, 0.0, ys, () -> DG.coef_pdf(D, ys, c; k...))
    val("coefq", id, cb, st, g0, c, c0, 0, -10.0, [0.05], () -> DG.coef_quantile(D, 0.05, c; kw4(cb, st, g0, c0)...); tol = 1e-7)
end
for (id, c) in ((2, 0.0), (9, -5.0))
    ys = [-12.0, -3.0]
    val("coefgp", id, NaN, 0, 0.0, c, 0.0, 0, 0.0, ys, () -> DG.coef_cdf_gp(designs[id], ys, c))
end
# moments
for (id, cb, st, g0, c, c0) in ((3, NaN, 0, 0.0, -5.0, 0.0), (2, NaN, 0, 0.0, -5.0, -5.0), (9, NaN, 0, 0.0, 0.0, 0.0), (3, NaN, 2, 0.0, -5.0, 0.0), (2, NaN, 1, 1.0, -5.0, 0.0),
                                (3, -13.5, 0, 0.0, 0.0, 0.0), (12, NaN, 0, 0.0, -13.5, 0.0), (1, NaN, 0, 0.0, -50.0, 0.0))
    D = designs[id]
    val("mom", id, cb, st, g0, c, c0, 0, 4.0, [1.0, 2.0, 3.0, 4.0], () -> DG.dist_moments(D, c; p = 4, kw4(cb, st, g0, c0)...)[1])
end
# local power
val("lpow", 3, NaN, 0, 0.0, 0.0, 0.0, 0, 0.0, [-5.0, -13.5], () -> DG.local_power(designs[3], [-5.0, -13.5]); tol = 1e-8)
val("lpow", 9, NaN, 0, 0.0, 0.0, 0.0, 0, 1.0, [-10.0], () -> DG.local_power(designs[9], [-10.0]; statistic = :coef); tol = 1e-8)
# the linear form, the law of V, the edge of the support of U (b near -1/2)
for (id, cb, st, g0, c, a) in ((1, NaN, 0, 0.0, 0.0, -3.5), (1, NaN, 0, 0.0, -7.0, -1.0), (3, -13.5, 0, 0.0, 0.0, -7.3), (9, NaN, 0, 0.0, -5.0, 2.0), (2, NaN, 2, 0.0, -5.0, -2.0),
                               (3, NaN, 1, 1.0, -5.0, -2.0), (1, NaN, 0, 0.0, 0.0, 0.0), (3, -5.0, 0, 0.0, 0.0, -2.0))
    D = designs[id]; bs = [-0.6, -0.2639, 0.0, 0.5]; k = kw(cb, st, g0, 0.0, false)
    val("lincdf", id, cb, st, g0, c, 0.0, 0, a, bs, () -> DG.lin_cdf(D, bs, a, c; kwl(cb, st, g0)...))
    val("linpdf", id, cb, st, g0, c, 0.0, 0, a, bs, () -> DG.lin_pdf(D, bs, a, c; kwl(cb, st, g0)...))
end
val("lincdf", 1, NaN, 0, 0.0, 0.0, 0.0, 0, 0.0, [-0.499, -0.4999, -0.49999], () -> DG.lin_cdf(designs[1], [-0.499, -0.4999, -0.49999], 0.0, 0.0))
val("linq", 1, NaN, 0, 0.0, 0.0, 0.0, 0, -3.5, [0.05], () -> DG.lin_quantile(designs[1], 0.05, -3.5, 0.0); tol = 1e-7)
val("linq", 3, -13.5, 0, 0.0, 0.0, 0.0, 0, -7.3, [0.05], () -> DG.lin_quantile(designs[3], 0.05, -7.3, 0.0; gls = -13.5); tol = 1e-7)
for (id, cb, st, c) in ((1, NaN, 0, -5.0), (1, NaN, 0, 5.0), (9, NaN, 0, 0.0), (3, -13.5, 0, -5.0), (2, NaN, 2, -5.0))
    D = designs[id]; ys = [0.02, 0.06, 0.2, 0.6]
    val("lincdf", id, cb, st, 0.0, c, 0.0, 0, NaN, ys, () -> DG.v_cdf(D, ys, c; gls = isnan(cb) ? nothing : cb, start = startsym(st)))
    val("linpdf", id, cb, st, 0.0, c, 0.0, 0, NaN, ys, () -> DG.v_pdf(D, ys, c; gls = isnan(cb) ? nothing : cb, start = startsym(st)))
end
# ADFmax
for (id, st, c) in ((2, 0, 0.0), (3, 0, -10.0), (2, 2, -10.0), (9, 0, -5.0))
    D = designs[id]; xs = [-3.0, -2.4, -1.5]
    val("adfcdf", id, NaN, st, 0.0, c, 0.0, 0, 0.0, xs, () -> DG.adfmax_cdf(D, xs, c; start = startsym(st)))
    val("adfpdf", id, NaN, st, 0.0, c, 0.0, 0, 0.0, xs, () -> DG.adfmax_pdf(D, xs, c; start = startsym(st)))
end
val("adfq", 2, NaN, 0, 0.0, 0.0, 0.0, 0, -2.0, [0.05], () -> DG.adfmax_quantile(designs[2], 0.05, 0.0); tol = 1e-7)
# the power envelope
val("env", 1, NaN, 0, 0.0, 0.0, 0.0, 0, 0.05, [-5.0, -10.0], () -> DG.power_envelope(designs[1], [-5.0, -10.0]); tol = 1e-8)
val("env", 3, NaN, 0, 0.0, 0.0, 0.0, 0, 0.05, [-5.0, -10.0], () -> DG.power_envelope(designs[3], [-5.0, -10.0]); tol = 1e-8)
# the companions: the point-optimal statistic and the envelope, the LM statistic, MZ_alpha, MSB
for (id, cb, c) in ((3, -13.5, 0.0), (3, -13.5, -10.0), (2, -7.0, 0.0), (1, -7.0, -5.0), (20, -10.0, -5.0), (9, -10.0, 0.0))
    qs = [2.0, 4.0, 7.0, 12.0]
    val("ptcdf", id, cb, 0, 0.0, c, 0.0, 0, cb, qs, () -> DG.pt_cdf(designs[id], qs, cb, c))
    val("ptpdf", id, cb, 0, 0.0, c, 0.0, 0, cb, qs, () -> DG.pt_pdf(designs[id], qs, cb, c))
end
val("ptq", 3, -13.5, 0, 0.0, 0.0, 0.0, 0, -13.5, [0.05], () -> DG.pt_quantile(designs[3], 0.05, -13.5); tol = 1e-7)
val("ptpow", 20, -10.0, 0, 0.0, -10.0, 0.0, 0, -10.0, [0.05], () -> collect(DG.pt_power(designs[20], -10.0, -10.0)); tol = 1e-7)
val("ptenv", 9, NaN, 0, 0.0, 0.0, 0.0, 0, 0.05, [-5.0, -10.0, 1.0], () -> DG.pt_envelope(designs[9], [-5.0, -10.0, 1.0]); tol = 1e-8)
for (id, c) in ((3, 0.0), (3, -10.0), (8, -5.0), (9, 0.0), (10, -13.5))
    xs = [-3.0, -2.2, -1.6, -1.0, 0.5]
    val("lmcdf", id, NaN, 0, 0.0, c, 0.0, 0, 0.0, xs, () -> DG.lm_cdf(designs[id], xs, c))
end
val("lmq", 3, NaN, 0, 0.0, 0.0, 0.0, 0, 0.0, [0.05], () -> DG.lm_quantile(designs[3], 0.05); tol = 1e-7)
val("lmq", 9, NaN, 0, 0.0, -5.0, 0.0, 0, 0.0, [0.1], () -> DG.lm_quantile(designs[9], 0.1; c = -5.0); tol = 1e-7)
val("lmpow", 8, NaN, 0, 0.0, 0.0, 0.0, 0, 0.05, [-5.0, -10.0], () -> (r = DG.lm_power(designs[8], [-5.0, -10.0]); [r[1]; r[2]]); tol = 1e-7)
for (id, cb, c) in ((3, -13.5, 0.0), (3, -13.5, -10.0), (20, -13.5, -5.0), (5, -10.0, 0.0))
    ys = [-25.0, -12.0, -5.0, -1.0]
    val("mzacdf", id, cb, 0, 0.0, c, 0.0, 0, cb, ys, () -> DG.mza_cdf(designs[id], ys, c; gls = cb))
    xs = [0.08, 0.15, 0.25, 0.4]
    val("msbcdf", id, cb, 0, 0.0, c, 0.0, 0, cb, xs, () -> DG.msb_cdf(designs[id], xs, c; gls = cb))
end
val("mzaq", 3, -13.5, 0, 0.0, 0.0, 0.0, 0, -13.5, [0.05], () -> DG.mza_quantile(designs[3], 0.05, 0.0; gls = -13.5); tol = 1e-7)
val("msbq", 3, -13.5, 0, 0.0, 0.0, 0.0, 0, -13.5, [0.05], () -> DG.msb_quantile(designs[3], 0.05, 0.0; gls = -13.5); tol = 1e-7)
# ---------------------------------------------------------------------------------------------------------------------
# the router (merged.jl): 'route func engine spec a b c d c c0 par n x... nv v... tol', spec 'case' (a = model, b = start,
# c = g0) or 'design' (a = id, b = start, c = g0, d = gls); engine 0 auto, 1 classical, 2 general
const ENG = (:auto, :classical, :general)
casestr(m, st) = m == 4 ? "DG" : (st == 2 ? "S" : "D") * "KUT"[m]
function rline(func, e, spec, a, b, cc, d, c, c0, par, xs, f; tol = 1e-9)
    v, w = warned(f); v = vcat(v); DG.empty_caches!()
    println(io, join(["route", func, string(e), spec, string(a), string(b), fmt(cc), fmt(d), fmt(c), fmt(c0), fmt(par), string(length(xs)), fmt.(xs)...,
                      string(length(v)), fmt.(v)..., fmt(w ? max(tol, 1e-6) : tol)], " ") * (w ? " brittle" : ""))
end
for e in 0:2
    eng = ENG[e+1]
    for (m, st, g0, c, c0) in ((2, 0, 0.0, -5.0, 0.0), (1, 0, 0.0, -5.0, 2.0), (2, 2, 0.0, -10.0, -10.0), (1, 1, 1.0, -5.0, 0.0), (3, 0, 0.0, 0.0, 0.0), (4, 0, 0.0, -5.0, 0.0))
        xs = [-3.5, -2.5, -1.6, -0.5, 1.0]; cs = casestr(m, st)
        rline("cdf", e, "case", m, st, g0, NaN, c, c0, 0.0, xs, () -> DFLimit.ltratio_cdf(xs, cs, c; c0, g0, engine = eng))
        rline("pdf", e, "case", m, st, g0, NaN, c, c0, 0.0, xs, () -> DFLimit.ltratio_pdf(xs, cs, c; c0, g0, engine = eng))
    end
    for (id, st, g0, gl, c, c0) in ((3, 0, 0.0, NaN, -5.0, 0.0), (9, 0, 0.0, NaN, -5.0, 0.0), (3, 0, 0.0, -13.5, -5.0, 0.0), (2, 2, 0.0, NaN, -5.0, -5.0))
        e == 1 && id == 9 && continue
        xs = [-3.5, -2.5, -1.6, -0.5, 1.0]
        rline("cdf", e, "design", id, st, g0, gl, c, c0, 0.0, xs, () -> DFLimit.ltratio_cdf(xs, designs[id], c; c0, g0, start = startsym(st), gls = isnan(gl) ? nothing : gl, engine = eng))
    end
    for (m, st, g0, c) in ((3, 0, 0.0, -5.0), (1, 0, 0.0, 0.0), (2, 0, 1.0, -5.0))
        ps = [0.05, 0.5]; cs = casestr(m, st)
        rline("icdf", e, "case", m, st, g0, NaN, c, 0.0, 0.0, ps, () -> DFLimit.ltratio_icdf(ps, cs, c; g0, engine = eng); tol = 1e-7)
    end
    for (m, st, c, c0) in ((2, 0, -5.0, 0.0), (1, 0, -5.0, -5.0), (3, 2, -5.0, 0.0))
        rline("mom", e, "case", m, st, 0.0, NaN, c, c0, 0.0, [4.0], () -> DFLimit.ltratio_mom(4, casestr(m, st), c; c0, engine = eng); tol = 1e-8)
    end
    ys = [-20.0, -10.0, -5.0]
    rline("coef", e, "case", 2, 0, 0.0, NaN, 0.0, 0.0, 0.0, ys, () -> DFLimit.coef_cdf(ys, "DU", 0.0; engine = eng))
    rline("coefq", e, "case", 2, 0, 0.0, NaN, 0.0, 0.0, 0.0, [0.05], () -> DFLimit.coef_icdf(0.05, "DU"; engine = eng); tol = 1e-7)
    bs = [-0.6, -0.2639, 0.3]
    rline("lin", e, "case", 1, 0, 0.0, NaN, 0.0, 0.0, -3.5, bs, () -> DFLimit.ltratio_lincdf(bs, -3.5, "DK", 0.0; engine = eng))
    vs = [0.02, 0.06, 0.2]
    rline("lin", e, "case", 1, 0, 0.0, NaN, -5.0, 0.0, NaN, vs, () -> DFLimit.ltratio_lincdf(vs, NaN, "DK", -5.0; engine = eng))
    rline("lin", e, "case", 4, 0, 0.0, -13.5, 0.0, 0.0, -7.3, [-0.9, -0.5, 0.0], () -> DFLimit.ltratio_lincdf([-0.9, -0.5, 0.0], -7.3, "DG", 0.0; engine = eng))
    rline("env", e, "case", 1, 0, 0.0, NaN, 0.0, 0.0, 0.05, [-5.0, -10.0], () -> DFLimit.power_envelope([-5.0, -10.0], "DK"; engine = eng); tol = 1e-8)
    rline("env", e, "case", 4, 0, 0.0, NaN, 0.0, 0.0, 0.05, [-5.0], () -> DFLimit.power_envelope([-5.0], "DG"; engine = eng); tol = 1e-8)
    for (m, st, c) in ((2, 0, 0.0), (2, 2, -10.0))
        xs = [-3.0, -2.4, -1.5]; cs = casestr(m, st)
        rline("adf", e, "case", m, st, 0.0, NaN, c, 0.0, 0.0, xs, () -> DFLimit.adfmax_cdf(xs, cs, c; engine = eng))
        rline("adfpdf", e, "case", m, st, 0.0, NaN, c, 0.0, 0.0, xs, () -> DFLimit.adfmax_pdf(xs, cs, c; engine = eng))
    end
    rline("region1", e, "case", 1, 0, 0.0, NaN, -5.0, 0.0, -3.5, [-0.2639], () -> DFLimit.ltratio_region((type = "line", a = -3.5, b = -0.2639), "DK", -5.0; engine = eng))
    rline("region2", e, "case", 1, 0, 0.0, NaN, -7.0, 0.0, 0.0, [-1.94087], () -> DFLimit.ltratio_region((type = "parabola", x = -1.94087), "DK", -7.0; engine = eng))
end
# the rule: 'auto f model stat g0 gls c c0 x engine' with Julia's autoengine (0 classical, 1 general), f 1 dist, 2 icdf,
# 3 mom, 4 lincdf (x = a), 5 adfmax, 6 envelope
AF = (:dist, :icdf, :mom, :lincdf, :adfmax, :envelope)
function aline(fi, m, st, g0, gl, c, c0, x)
    cs = casestr(m, st)
    S = fi == 4 ? DFLimit.lin_spec(cs, isnan(gl) ? -13.5 : gl) : DFLimit.case_spec(cs, g0)
    e = DFLimit.autoengine(AF[fi], S, Float64(c), Float64(c0), Float64(x))
    println(io, join(["auto", string(fi), string(m), string(st), fmt(g0), fmt(gl), fmt(c), fmt(c0), fmt(x), e === :general ? "1" : "0"], " "))
end
cgrid = (-300.0, -100.0, -50.0, -30.0, -20.0, -10.0, -5.0, -1.0, -0.5, -1e-3, 0.0, 1.0, 2.0, 3.0, 5.0)
for m in 1:4, st in (0, 2), c in cgrid
    (m == 4 && st == 2) && continue
    st == 2 && c >= 0 && continue
    for c0 in unique((0.0, -5.0, 2.0, c))
        aline(1, m, st, 0.0, NaN, c, c0, -2.0)
        m == 1 && st == 0 && c0 > 0 && aline(1, m, st, 0.0, NaN, c, c0, -sqrt(2c0) + 5e-4)
        aline(2, m, st, 0.0, NaN, c, c0, 0.0)
        m != 4 && aline(3, m, st, 0.0, NaN, c, c0, 0.0)
    end
    if st == 0 && m != 4
        for g0 in (0.5, 1.0, -1.0, 3.0)
            aline(1, m, st, g0, NaN, c, 0.0, -2.0); aline(2, m, st, g0, NaN, c, 0.0, 0.0); aline(3, m, st, g0, NaN, c, 0.0, 0.0)
        end
    end
    m in (2, 3) && aline(5, m, st, 0.0, NaN, c, 0.0, 0.0)
end
for c in cgrid, a in (NaN, -3.5, 0.0); aline(4, 1, 0, 0.0, NaN, c, 0.0, a); end
for c in cgrid, gl in (-13.5, -5.0, -2.0, -1.0), a in (-7.3, 0.0); aline(4, 4, 0, 0.0, gl, c, 0.0, a); end
for m in (1, 4), c in (-10.0, -1.0, -0.5, -0.01, -1e-4); aline(6, m, 0, 0.0, NaN, c, 0.0, 0.0); end
close(io)
println("written ", OUT)
