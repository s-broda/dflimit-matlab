# reference values of the devices of Bao and Kan for the C engine test (c/dfl_test_kan.c): c/test/kan_ref.txt, from the module
# DFLimit.DFKan of DFLimit.jl (src/DFKan.jl, src/DFKan_ratio.jl, src/DFKan_det.jl) and the package's user-level functions
# (src/merged.jl).  Run from the repository root: julia --project=c/test c/test/kan_ref.jl c/test/kan_ref.txt
#
# Lines (model 1, 2, 3 = K, U, T; stat 1 a stationary start; g0 the fixed start):
#   lcf model stat c g0 t1r t1i t2r t2i Lr Li                      log E exp(i t1 A + i t2 B) (DFKan.logcf)
#   jet model stat c g0 n t1 (n complex) t2 (n complex) L (n complex)  the same in Taylor-mode arithmetic
#   scale model stat c g0 Eu Ev sdv Vu                             ratio_scale
#   coef model stat c g0 c0 y P e                                  coef_cdf1 (tol 1e-12)
#   coefpdf model stat c g0 c0 y f e                               coef_pdf1 (tol 1e-12)
#   coefq model stat c g0 c0 p y                                   coef_icdf1 (xtol 1e-9, tol 1e-12)
#   lin model stat c g0 a b P e                                    lin_cdf1 (tol 1e-12)
#   env alpha c P e                                                envelope1 (tol 1e-12)
#   tmom model stat c g0 c0 k m_1 .. m_k                           tratio_mom1 (tol 1e-10)
#   cmom model stat c g0 c0 k m_1 .. m_k                           coef_mom1 (tol 1e-10)
#   ray model stat c g0 P                                          ray_cdf
# and the user-level functions with their defaults (the router, dfl_engine_*; start 0 zero, 1 fixed, 2 stationary):
#   rcoef model start g0 c c0 n y_1 .. y_n P_1 .. P_n               coef_cdf(y, cas, c; c0, g0)
#   rcoefpdf model start g0 c c0 n y.. f..                         coef_pdf
#   rcoefq model start g0 c c0 n p.. y..                           coef_icdf (tol 1e-9)
#   rcoefmom model start g0 c c0 k m..                             coef_mom
#   rmom model start g0 c c0 k m..                                 ltratio_mom
#   rlin p start g0 c a n b.. P..                                  ltratio_lincdf(b, a, trend(p), c; g0, start) (p = -1 with a
#                                                                  zero start: also the case code DK, model 1)
#   renv alpha n c.. P..                                           power_envelope(c, "DK"; alpha)
# and the general engine (a design: D kind np par.., kind 0 trend(par), 1, 2, 3 Perron's models A, B, C with np breaks):
#   dcf D.. start gls c t1 t2 Lr Li                               DFDet.DetCF: logcf with the branch from the pivots
#   ptcdf D.. cbar c n q.. P..                                     pt_cdf (the default: the single inversion up to c = 4)
#   ptpdf D.. cbar c n q.. f..                                     pt_pdf
#   ptq D.. cbar p q                                               pt_quantile
#   ptpow D.. cbar c level b power                                 pt_power
#   ptenv D.. level n c.. P..                                      pt_envelope
#   gp D.. start gls c c0 n y.. P..                                coef_cdf_gp (the default: the single inversion up to c = 5)
using Printf
using DFLimit
const K = DFLimit.DFKan
out = open(ARGS[1], "w")
fmt(x::Real) = @sprintf("%.17g", x)
fmt(z::Complex) = fmt(real(z)) * " " * fmt(imag(z))
CASES = (("DK", 1, 0), ("DU", 2, 0), ("DT", 3, 0), ("SK", 1, 1), ("SU", 2, 1), ("ST", 3, 1))
casecode(m, st) = (st == 1 ? "S" : "D") * "KUT"[m]
line(xs...) = println(out, join(xs, " "))

# ---- the transform at points, scalars and jets ----
pts = ((0.3, -0.2), (4.0, 30.0), (-25.0, 700.0), (60.0, -2000.0), (0.7, -1.3), (1e-3, 2e-3), (400.0, 1200.0), (-0.05, 1.5), (3.0, -0.5im), (-0.2im, 1.5im), (2.0 - 0.1im, 0.5 + 0.3im))
for (cas, m, st) in CASES, c in (1.5, 0.0, -0.7, -6.0, -60.0, 5.0, 20.0), g0 in (0.0, 1.0)
    (st == 1 && (c >= 0 || g0 != 0)) && continue
    F = K.limcf(cas, c, g0)
    for (t1, t2) in pts
        L = K.logcf(F, complex(t1), complex(t2))
        line("lcf", m, st, fmt(c), fmt(g0), fmt(complex(t1)), fmt(complex(t2)), fmt(L))
    end
    for (n, a1, a0, b1, b0) in ((9, 1.0, 37.0, -2.5, -2.5 * 37.0), (3, 0.0, 0.0, 1.0, 0.0), (5, -1im, 0.0, 0.5im, 3.0im), (2, 0.0, 12.0, 1.0, -30.0))
        t1 = K.Jet{n}(K.SVector{n,ComplexF64}(ntuple(k -> k == 1 ? ComplexF64(a0) : k == 2 ? ComplexF64(a1) : 0.0im, n)))
        t2 = K.Jet{n}(K.SVector{n,ComplexF64}(ntuple(k -> k == 1 ? ComplexF64(b0) : k == 2 ? ComplexF64(b1) : 0.0im, n)))
        L = K.logcf(F, t1, t2)
        line("jet", m, st, fmt(c), fmt(g0), n, join(fmt.(t1.c), " "), join(fmt.(t2.c), " "), join(fmt.(L.c), " "))
    end
    line("scale", m, st, fmt(c), fmt(g0), join(fmt.(collect(K.ratio_scale(F))), " "))
end

# ---- the single inversion: the coefficient statistic, its density and quantiles, the linear form ----
for (cas, m, st) in CASES, c in (-60.0, -20.0, -5.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, 20.0), g0 in (0.0, 1.0)
    (st == 1 && (c >= 0 || g0 != 0)) && continue
    (g0 != 0 && c > 2) && continue
    F = K.limcf(cas, c, g0)
    for (i, p) in enumerate((0.01, 0.5, 0.99))
        r = K.ratio_icdf(F, p)
        c0 = i == 2 ? 2.0 : 0.0
        y = r - (c0 - c)
        P, e = K.coef_cdf1(F, y, c0, 1e-12)
        line("coef", m, st, fmt(c), fmt(g0), fmt(c0), fmt(y), fmt(P), fmt(e))
        f, ef, _ = K.coef_pdf1(F, y, c0, 1e-12)
        line("coefpdf", m, st, fmt(c), fmt(g0), fmt(c0), fmt(y), fmt(f), fmt(ef))
    end
    for p in (0.05, 0.9)
        y, _ = K.coef_icdf1(F, p, 0.0, 1e-9, 1e-12)
        line("coefq", m, st, fmt(c), fmt(g0), fmt(0.0), fmt(p), fmt(y))
    end
    if -100 <= c <= 4
        Eu, Ev, sdv, Vu = K.ratio_scale(F)
        for a in (-3.5, 1.0)
            bq = c - a; mu = Eu + bq * Ev; sd = sqrt(Vu + bq^2 * sdv^2)
            for p in (0.01, 0.5, 0.99)
                b = K.illinois(x -> K.lin_cdf1(F, x, a, 1e-12)[1] - p, mu - sd, mu + sd, 1e-9 * sd)
                P, e = K.lin_cdf1(F, b, a, 1e-12)
                line("lin", m, st, fmt(c), fmt(g0), fmt(a), fmt(b), fmt(P), fmt(e))
            end
        end
    end
    for c0 in unique((0.0, c, -2.0, 2.0))
        mm, _ = K.tratio_mom1(F, 4, c0, 1e-10)
        line("tmom", m, st, fmt(c), fmt(g0), fmt(c0), 4, join(fmt.(mm), " "))
    end
    mc, _ = K.coef_mom1(F, 2, 0.0, 1e-10)
    line("cmom", m, st, fmt(c), fmt(g0), fmt(0.0), 2, join(fmt.(mc), " "))
end
for c in (-1e-4, -0.01, -1.0, -7.0, -13.5, -30.0, -100.0, -500.0, 0.0), alpha in (0.05, 0.01)
    P, e = K.envelope1(c, alpha, 1e-12)
    line("env", fmt(alpha), fmt(c), fmt(P), fmt(e))
end
for (cas, c, g0) in (("DK", 2.0, 0.0), ("DK", 2.0, 1.0), ("DK", 0.0, 0.0), ("DK", 0.0, 1.0), ("DK", -0.1, 1.0), ("DK", -5.0, 0.0), ("DK", -5.0, 1.0), ("DK", -50.0, 0.0),
                     ("DU", 2.0, 0.0), ("DU", 2.0, 1.0), ("DU", 0.0, 0.0), ("DU", -0.1, 1.0), ("DU", -5.0, 0.0), ("DU", -5.0, 1.0), ("DU", -50.0, 0.0), ("DU", -50.0, 1.0),
                     ("SK", -0.1, 0.0), ("SK", -5.0, 0.0), ("SK", -50.0, 0.0), ("SU", -0.1, 0.0), ("SU", -5.0, 0.0), ("SU", -50.0, 0.0))
    line("ray", cas[2] == 'K' ? 1 : 2, cas[1] == 'S' ? 1 : 0, fmt(c), fmt(g0), fmt(K.ray_cdf(cas, c; g0)))
end

# ---- the user-level functions (src/merged.jl) with their defaults ----
vec(xs) = join(fmt.(xs), " ")
ST = (zero = 0, fixed = 1, stationary = 2)
for (cas, c, g0, ys, c0) in (("DK", -5.0, 0.0, [-9.0, -3.0], 0.0), ("DU", -5.0, 1.0, [-12.0, -6.0], 0.0), ("DT", 1.0, 0.0, [-10.0, 0.5], 0.0), ("SK", -1.0, 0.0, [-3.0], 0.0),
                             ("SU", -10.0, 0.0, [-20.0, -8.0], 0.0), ("ST", -5.0, 0.0, [-25.0], 0.0), ("DU", -5.0, 0.0, [-14.0], 0.0), ("DK", 0.0, 0.0, [-8.039, 1.0, -Inf, NaN, Inf], 0.0),
                             ("DT", -2.0, 0.0, [-4.0], 0.0), ("DU", -60.0, 1.5, [-30.0, -10.0], 0.0), ("DK", 20.0, 0.0, [-0.5, 0.5], 20.0), ("DU", -20.0, 0.0, [-12.0], -3.0),
                             ("DT", 2.0, 0.5, [-15.0], 2.0), ("SU", -100.0, 0.0, [-15.0], -100.0))
    m = findfirst(==(cas[2]), "KUT"); st = cas[1] == 'S' ? 2 : g0 != 0 ? 1 : 0
    line("rcoef", m, st, fmt(g0), fmt(c), fmt(c0), length(ys), vec(ys), vec(coef_cdf(ys, cas, c; c0, g0)))
    yf = filter(isfinite, ys)
    line("rcoefpdf", m, st, fmt(g0), fmt(c), fmt(c0), length(yf), vec(yf), vec(coef_pdf(yf, cas, c; c0, g0)))
    ps = [0.05, 0.5]
    line("rcoefq", m, st, fmt(g0), fmt(c), fmt(c0), 2, vec(ps), vec(coef_icdf(ps, cas, c; c0, g0)))
    line("rcoefmom", m, st, fmt(g0), fmt(c), fmt(c0), 3, vec(coef_mom(3, cas, c; c0, g0)))
end
for (cas, c, c0, g0) in (("DK", 0.0, 0.0, 0.0), ("DU", -5.0, 0.0, 0.0), ("DT", -5.0, -5.0, 0.0), ("SU", -3.0, 0.0, 0.0), ("DK", 2.0, -2.0, 1.0), ("DK", 20.0, 20.0, 0.0),
                     ("DU", -100.0, -100.0, 0.0), ("ST", -20.0, 2.0, 0.0), ("DT", 1.0, 0.0, 1.0), ("SK", -0.5, -0.5, 0.0), ("DU", 10.0, 10.0, 0.0))
    m = findfirst(==(cas[2]), "KUT"); st = cas[1] == 'S' ? 2 : g0 != 0 ? 1 : 0
    line("rmom", m, st, fmt(g0), fmt(c), fmt(c0), 4, vec(ltratio_mom(4, cas, c; c0, g0)))
end
for (p, st, g0, c, a, bs) in ((-1, :zero, 0.0, 0.0, -3.5, [-0.2639, 0.5]), (-1, :zero, 0.0, -5.0, -2.5, [0.3]), (-1, :zero, 0.0, -13.5, 1.0, [-0.1]), (-1, :zero, 0.0, 2.0, -3.5, [5.0]),
                              (-1, :zero, 0.0, 4.0, -1.0, [2.0, 20.0]), (-1, :zero, 0.0, -100.0, -50.0, [-0.3]), (0, :zero, 0.0, -5.0, -3.5, [-0.2639]), (1, :zero, 0.0, -5.0, -3.5, [-0.2639]),
                              (0, :stationary, 0.0, -5.0, -3.5, [-0.2639, 0.4]), (-1, :fixed, 1.5, -5.0, -3.5, [-0.2639]), (-1, :fixed, 1.0, 1.0, -7.0, [1.0, 30.0]), (0, :fixed, 1.0, 2.0, -3.5, [2.0]),
                              (1, :stationary, 0.0, -20.0, -10.0, [-0.3, 0.2]), (-1, :stationary, 0.0, -1.0, 2.0, [-0.5]))
    P = ltratio_lincdf(bs, a, trend(p), c; g0, start = st)
    line("rlin", p, ST[st], fmt(g0), fmt(c), fmt(a), length(bs), vec(bs), vec(P))
end
for (alpha, cs) in ((0.05, [-1e-4, -0.01, -1.0, -7.0, -13.5, -30.0, -500.0, 0.0, 1.0]), (0.01, [-10.0]), (0.1, [-2.0, -25.0]))
    line("renv", fmt(alpha), length(cs), vec(cs), vec(power_envelope(cs, "DK"; alpha)))
end
# ---- the general engine: the transform with the branch from the pivots, the point-optimal statistic, coef_cdf_gp ----
G = DFLimit.DFDet
dspec(kind, par) = kind == 0 ? ("D 0 1 " * fmt(par[1]), G.trend(Int(par[1]))) :
    ("D $kind $(length(par)) " * join(fmt.(par), " "), G.perron((:A, :B, :C)[kind], length(par) == 1 ? par[1] : collect(par)))
for (kind, par, st, cs) in ((0, (-1,), :zero, (-4.0, 0.0, 2.0)), (0, (0,), :zero, (-4.0, 0.0)), (0, (1,), :zero, (-4.0, 0.0, 3.0)), (0, (0,), :stationary, (-4.0, -20.0)),
                            (0, (1,), :stationary, (-4.0,)), (0, (2,), :zero, (-5.0,)), (3, (0.5,), :zero, (-5.0, 1.0)), (2, (0.3,), :stationary, (-10.0,)))
    ds, D = dspec(kind, par)
    for c in cs
        F = G.DetCF(G.startmodel(D, c, st, 0.0, nothing), c)
        for (t1, t2) in ((0.3, -0.2), (4.0, 30.0), (-25.0, 700.0), (60.0, -2000.0), (1e-3, 0.5))
            line("dcf", ds, ST[st], fmt(NaN), fmt(c), fmt(t1), fmt(t2), fmt(K.logcf(F, t1, t2)))
        end
    end
end
let (ds, D) = dspec(0, (1,)), c = -2.0
    F = G.DetCF(G.startmodel(D, c, :zero, 0.0, -13.5), c)
    for (t1, t2) in ((0.3, -0.2), (4.0, 30.0), (-25.0, 700.0))
        line("dcf", ds, 0, fmt(-13.5), fmt(c), fmt(t1), fmt(t2), fmt(K.logcf(F, t1, t2)))
    end
end
for (kind, par, cbar, c, qs) in ((0, (0,), -7.0, 0.0, [18.1014815, 3.0]), (0, (1,), -13.5, 4.0, [969.3008382]), (2, (0.5,), -17.8, 2.0, [16.5697]), (0, (0,), -7.0, 4.0, [2223.714162]),
                                 (0, (1,), -13.5, -30.0, [1.5, 4.0]), (0, (1,), -13.5, 0.0, [5.0, 10.0, 2000.0]), (0, (2,), -20.0, -5.0, [6.0, 12.0]), (3, (0.3, 0.7), -20.0, -10.0, [5.0]))
    ds, D = dspec(kind, par)
    line("ptcdf", ds, fmt(cbar), fmt(c), length(qs), vec(qs), vec(pt_cdf(D, qs, cbar, c)))
    line("ptpdf", ds, fmt(cbar), fmt(c), length(qs), vec(qs), vec(pt_pdf(D, qs, cbar, c)))
end
for (kind, par, cbar, p) in ((0, (1,), -13.5, 0.05), (0, (0,), -7.0, 0.05), (2, (0.5,), -17.8, 0.1), (0, (2,), -20.0, 0.5))
    ds, D = dspec(kind, par)
    line("ptq", ds, fmt(cbar), fmt(p), fmt(pt_quantile(D, p, cbar)))
end
for (kind, par, cbar, c, lev) in ((0, (1,), -13.5, -13.5, 0.05), (0, (0,), -7.0, -3.0, 0.05), (3, (0.5,), -20.0, -20.0, 0.1))
    ds, D = dspec(kind, par)
    b, P = pt_power(D, cbar, c; level = lev)
    line("ptpow", ds, fmt(cbar), fmt(c), fmt(lev), fmt(b), fmt(P))
end
for (kind, par, lev, cs) in ((0, (0,), 0.05, [-7.0, NaN, 0.5]), (2, (0.5,), 0.05, [-17.8, -5.0]), (0, (2,), 0.1, [-20.0]))
    ds, D = dspec(kind, par)
    line("ptenv", ds, fmt(lev), length(cs), vec(cs), vec(pt_envelope(D, cs; level = lev)))
end
for (kind, par, st, gls, c, c0, ys) in ((0, (0,), :zero, nothing, -5.0, 0.0, [-13.83853904, -5.0]), (0, (2,), :zero, nothing, 0.0, 0.0, [-20.0]), (3, (0.3, 0.7), :zero, nothing, -13.5, 0.0, [-15.0]),
                                       (0, (1,), :stationary, nothing, -13.5, 0.0, [-20.0]), (0, (0,), :zero, nothing, 5.0, 0.0, [9.942147816]), (3, (0.5,), :zero, nothing, -5.0, 0.0, [-15.0]),
                                       (0, (1,), :zero, -13.5, 0.0, 0.0, [-20.0]), (0, (-1,), :zero, nothing, -5.0, 2.0, [-10.0, Inf, NaN]))
    ds, D = dspec(kind, par)
    P = coef_cdf_gp(D, ys, c; c0, start = st, gls)
    line("gp", ds, ST[st], fmt(gls === nothing ? NaN : gls), fmt(c), fmt(c0), length(ys), vec(ys), vec(P))
end
close(out)
