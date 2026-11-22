# The reference values of the general engine for the MATLAB programs (dfdet_test.m): the Julia engine DFLimit.DFDet and
# the merged interface of DFLimit.jl (branch engine-merge) on the designs, starts and functions of the MATLAB port.
# Writes det_ref.txt beside this file, one value set per line:
#
#     <func> <design> <c> <options> | <arguments> | <values>
#
# <design> a token of dfdet_test.m (t-1 .. t4 the polynomial trends, pA0.4, pB0.3, pC0.3, pC0.3,0.7, pA0.25,0.5,0.75,
# pC0.25,0.5,0.75 Perron's models, fourier, lt3, lt8, dfm, dfmb the designs from a matrix, df the statistic that detrends
# first, lev, fdC the function designs; DK .. DG the case strings of the router), <options> key=value pairs separated by
# commas or "-", the numbers in full precision.  Run from the repository root as
#
#     JLCAP_GB=4 jlcap --project=/path/to/DFLimit.jl matlab/validation/det_ref.jl
#
# with an environment in which DFLimit.jl (branch engine-merge) is instantiated (or one that develops it).
using DFLimit, Printf, LinearAlgebra, Logging
const G = DFLimit.DFDet
const OUT = joinpath(@__DIR__, "det_ref.txt")
fm(x::Real) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Inf" : "-Inf") : @sprintf("%.17g", x)
fm(z::Complex) = fm(real(z)) * " " * fm(imag(z))
quiet(f) = with_logger(NullLogger()) do; f(); end

designs = Dict{String,Any}(
    "t-1" => G.trend(-1), "t0" => G.trend(0), "t1" => G.trend(1), "t2" => G.trend(2), "t3" => G.trend(3), "t4" => G.trend(4),
    "pA0.4" => G.perron(:A, 0.4), "pB0.3" => G.perron(:B, 0.3), "pC0.3" => G.perron(:C, 0.3), "pB0.5" => G.perron(:B, 0.5),
    "pB0.1" => G.perron(:B, 0.1), "pC0.5" => G.perron(:C, 0.5),
    "pC0.3,0.7" => G.perron(:C, [0.3, 0.7]), "pA0.25,0.5,0.75" => G.perron(:A, [0.25, 0.5, 0.75]), "pC0.25,0.5,0.75" => G.perron(:C, [0.25, 0.5, 0.75]),
    "fourier" => G.fdesign(r -> [1.0, sin(2pi * r), cos(2pi * r)]), "lt3" => G.ltrend(3), "lt8" => G.ltrend(8),
    "dfm" => G.design_from_matrix(G.regmat(G.perron(:C, 0.3), 500)),
    "dfmb" => G.design_from_matrix(hcat(ones(400), (1:400) ./ 400, [t > 120 ? 1.0 : 0.0 for t in 1:400]); breaks = [120]),
    "df" => G.detrend_first(G.perron(:B, 0.5)),
    "lev" => G.Determ([0.0, 0.4, 1.0], [[1.0 1.0; 1.0 0.0], [1.0 0.0; 1.0 0.0]], "1 + r, 1{r < 0.4}"),
    "fdC" => G.fdesign(r -> [1.0, r, r > 0.4 ? 1.0 : 0.0, max(r - 0.4, 0.0)], [0.4]))
io = open(OUT, "w")
function okw(kw)
    isempty(kw) && return "-"
    join([string(k) * "=" * (v isa Symbol ? string(v) : fm(Float64(v))) for (k, v) in pairs(kw)], ",")
end
line(func, d, c, kw, args, vals) = println(io, join([func, d, fm(Float64(c)), okw(kw), "|", join(fm.(args), " "), "|", join(fm.(vals), " ")], " "))
gkw(kw) = Dict{Symbol,Any}((k == :g0 ? :gamma0 : k) => v for (k, v) in pairs(kw))

# ---------- the designs and their preparation ----------
for d in ("t-1", "t0", "t1", "t2", "t3", "t4", "pA0.4", "pB0.3", "pC0.3", "pC0.3,0.7", "pA0.25,0.5,0.75", "pC0.25,0.5,0.75", "fourier", "lt3", "lt8",
          "dfm", "dfmb", "df", "lev", "fdC")
    P = G.prep(designs[d])
    line("prep", d, 0, (;), Float64[], Float64[P.N, P.r, P.rp, P.rm, P.nj, P.p])
    line("dv", d, 0, (;), Float64[], sort(P.Dv))
end
for (d, cb) in (("t1", -13.5), ("pA0.4", -13.5), ("t0", -7.0), ("pB0.3", -10.0))
    P = G.prep_gls(designs[d], cb)
    line("prep", d, 0, (gls = cb,), Float64[], Float64[P.N, P.r, P.rp, P.rm, P.nj, P.p])
end
# ---------- the determinant, the roots, t_* ----------
for (d, c, kw) in (("t1", -5.0, (;)), ("pC0.3", 0.0, (;)), ("pC0.3,0.7", -13.5, (;)), ("t2", 2.0, (;)), ("t1", -5.0, (start = :stationary,)),
                   ("t0", -5.0, (start = :fixed, g0 = 1.5)), ("t1", -5.0, (gls = -13.5,)), ("lt3", -20.0, (;)), ("pA0.25,0.5,0.75", -50.0, (;)))
    Mo = G.startmodel(designs[d], c, get(kw, :start, :zero), get(kw, :g0, 0.0), get(kw, :gls, nothing))
    for (s, t) in ((0.3 + 2.0im, 0.5 - 7.0im), (-0.8 - 5.0im, -3.0 + 40.0im))
        dv = G.Dval(Mo, c, s, t)
        line("dval", d, c, kw, [real(s), imag(s), real(t), imag(t)], [real(dv), imag(dv)])
    end
    for t in (1.0 + 3.0im, -20.0 + 50.0im)
        _, ld, μ, _ = G.mgfparts(Mo, c, t, G.MWork(Mo))
        rts = sort(1 ./ μ; by = z -> (real(z), imag(z)))
        line("roots", d, c, kw, [real(t), imag(t)], reduce(vcat, [[real(z), imag(z)] for z in rts]))
    end
    line("tstar", d, c, kw, Float64[], [G.tstar(Mo, c)[1]])
end
# ---------- distribution function and density ----------
cdfcases = Any[
    ("t-1", 0.0, (;), [-2.5, -1.95, 0.0, 1.5]), ("t-1", -13.5, (;), [-4.0, -2.0, 0.0]),
    ("t0", 0.0, (;), [-6.0, -3.5, -2.86, -1.0, 0.5, 3.0]), ("t1", -5.0, (;), [-6.5, -3.5, -2.0, -0.7, 0.4, 2.5]),
    ("t2", -20.0, (;), [-6.0, -5.0, -4.0]), ("t3", 0.0, (;), [-5.0, -4.0, -3.0]), ("t4", -5.0, (;), [-5.5, -4.5, -3.5]),
    ("t1", -13.5, (;), [-5.0, -3.5, -2.0]), ("t1", -50.0, (;), [-10.0, -8.0, -6.0]), ("t1", 2.0, (;), [-3.0, -1.0, 1.0]),
    ("pA0.4", -5.0, (;), [-5.0, -4.0, -3.0]), ("pB0.3", 0.0, (;), [-5.0, -4.0, -3.0]),
    ("pC0.3", 0.0, (;), [-5.0, -4.0, -3.0, -2.0]), ("pC0.3", -5.0, (;), [-5.0, -4.0, -2.0]), ("pC0.3", -13.5, (;), [-6.0, -4.5, -3.0]),
    ("pC0.3", -20.0, (;), [-6.5, -5.0, -3.5]), ("pC0.3", -50.0, (;), [-11.0, -9.0, -7.0]), ("pC0.3", 2.0, (;), [-4.0, -2.0, 0.0]),
    ("pC0.3,0.7", 0.0, (;), [-6.0, -5.0, -4.0]), ("pA0.25,0.5,0.75", -5.0, (;), [-5.0, -4.0]), ("pC0.25,0.5,0.75", -5.0, (;), [-6.0, -5.0, -4.0]),
    ("fourier", 0.0, (;), [-5.0, -4.0, -3.0]), ("lt3", 0.0, (;), [-4.0, -3.0, -2.0]), ("lt8", -5.0, (;), [-6.0, -5.0]),
    ("dfm", 0.0, (;), [-4.0, -3.0, -2.0]), ("dfmb", -5.0, (;), [-4.5, -3.5]), ("df", -5.0, (;), [-4.0, -3.0]),
    ("fdC", -5.0, (;), [-4.0, -3.0]),
    # the starts
    ("t1", -10.0, (start = :stationary,), [-3.5, -2.0, -0.7]), ("pC0.3", -5.0, (start = :stationary,), [-5.0, -4.0, -3.0]),
    ("t0", -5.0, (start = :fixed, g0 = 1.5), [-3.5, -2.0, 0.4]), ("t-1", 0.0, (start = :fixed, g0 = 1.0), [-2.0, -1.0, 0.0, 1.0]),
    ("pC0.3", -5.0, (start = :fixed, g0 = 1.0), [-5.0, -4.0, -3.0]), ("t1", -20.0, (start = :fixed, g0 = 2.0), [-5.0, -4.0, -3.0]),
    # GLS detrending
    ("t1", -5.0, (gls = -13.5,), [-3.0, -2.0, -1.0, 0.5]), ("t1", 0.0, (gls = -13.5,), [-3.5, -2.8, -2.0]), ("pB0.3", -5.0, (gls = -10.0,), [-4.0, -3.0]),
    # a local null value
    ("t1", -5.0, (c0 = 3.0,), [-3.5, -2.0, 0.4, 2.5]), ("t0", -5.0, (c0 = -3.0,), [-3.5, -0.7, 2.5]), ("pC0.3", -5.0, (c0 = -5.0,), [-2.0, 0.0, 2.0]),
    ("t-1", 0.0, (c0 = 2.0,), [-2.5, -2.0 - 1e-3, -1.5]), ("t-1", 5.0, (c0 = -3.0,), [22.91, 94.61]),
    # far from the unit root, far tails
    ("t0", -300.0, (;), [-25.0, -24.0, -23.0]), ("t1", 0.0, (;), [-7.0, 1.5, 3.0])]
for (d, c, kw, xs) in cdfcases
    line("cdf", d, c, kw, xs, quiet(() -> G.dist_cdf(designs[d], xs, c; gkw(kw)...)))
    line("pdf", d, c, kw, xs, quiet(() -> G.dist_pdf(designs[d], xs, c; gkw(kw)...)))
end
# the design (1 + r, 1{r < 0.4}) without a proper intercept: on the vertical line a positive root runs off, and its U-shapes
# would need more than 1e6 nodes at every c (the guard of DFLimit.jl ac3b4d8; before it, 0.2064 and 0.7434 at c = -5 by
# U-shapes of up to 1e6 nodes): an error ('err', fn the function, txt the text of the message, '_' for ' ')
line("err", "lev", -5.0, (fn = :cdf, txt = Symbol("U-shapes_of_the_zero_start")), [-3.0, -2.0], Float64[])
line("err", "lev", 1.0, (fn = :cdf, txt = Symbol("U-shapes_of_the_zero_start")), [-2.0], Float64[])
# the infinite density at the tangency of case DK with c0 > 0
line("pdf", "t-1", 0.0, (c0 = 2.0,), [-2.0], quiet(() -> G.dist_pdf(designs["t-1"], [-2.0], 0.0; c0 = 2.0)))
# ---------- quantiles ----------
for (d, c, kw) in (("t0", 0.0, (;)), ("pC0.3", -5.0, (;)), ("df", 0.0, (;)), ("t1", 0.0, (gls = -13.5,)), ("t-1", -5.0, (start = :fixed, g0 = 1.0)),
                   ("pA0.25,0.5,0.75", 0.0, (;)), ("t2", -13.5, (;)), ("t1", -10.0, (start = :stationary,)), ("pC0.3", 0.0, (c0 = -2.0,)))
    ps = [0.01, 0.05, 0.10, 0.5, 0.95]
    line("icdf", d, c, kw, ps, [quiet(() -> G.dist_quantile(designs[d], p, c; x0 = -2.0, tol = 1e-12, gkw(kw)...)) for p in ps])
end
# ---------- moments ----------
for (d, c, kw) in (("t0", -5.0, (;)), ("pC0.3", -5.0, (;)), ("t1", -5.0, (start = :stationary,)), ("t0", -5.0, (start = :fixed, g0 = 1.0)),
                   ("t0", -5.0, (c0 = -5.0,)), ("t1", 0.0, (gls = -13.5,)), ("t3", 0.0, (;)), ("pA0.25,0.5,0.75", -20.0, (;)))
    line("mom", d, c, kw, Float64[], quiet(() -> G.dist_moments(designs[d], c; p = 4, gkw(kw)...)[1]))
end
# ---------- the related statistics ----------
ys = [-20.0, -8.0, -1.0]
for (d, c, kw) in (("t1", -5.0, (;)), ("pC0.3", 0.0, (;)), ("t1", -5.0, (gls = -13.5,)))
    line("coefcdf", d, c, kw, ys, G.coef_cdf(designs[d], ys, c; gkw(kw)...))
    line("coefpdf", d, c, kw, ys, G.coef_pdf(designs[d], ys, c; gkw(kw)...))
end
line("coeficdf", "t-1", 0.0, (;), [0.05], [G.coef_quantile(designs["t-1"], 0.05; tol = 1e-12)])
line("coeficdf", "pC0.3", 0.0, (;), [0.05], [G.coef_quantile(designs["pC0.3"], 0.05; y0 = -25.0, tol = 1e-12)])
bs = [-2.0, -0.3, 0.0, 0.5, 2.0]
for (d, a, c, kw) in (("t-1", -3.5, -7.0, (;)), ("t1", 1.0, 0.0, (gls = -13.5,)), ("pC0.3", -2.0, -5.0, (;)), ("t-1", -1.0, 5.0, (;)))
    line("lincdf", d, c, merge(kw, (a = a,)), bs, quiet(() -> G.lin_cdf(designs[d], bs, a, c; gkw(kw)...)))
    line("linpdf", d, c, merge(kw, (a = a,)), bs, quiet(() -> G.lin_pdf(designs[d], bs, a, c; gkw(kw)...)))
end
# the law of U at the edge of its support (a = 0, b next to -1/2)
be = [-0.5 + 1e-7, -0.49999, -0.4999, -0.499, -0.4, 0.0, 1.0]
for c in (0.0, -5.0, 2.0)
    line("lincdf", "t-1", c, (a = 0.0,), be, G.lin_cdf(designs["t-1"], be, 0.0, c))
end
line("linicdf", "t-1", 0.0, (a = -3.5,), [0.05], [G.lin_quantile(designs["t-1"], 0.05, -3.5, 0.0; tol = 1e-12)])
for (d, c, kw, yv) in (("t-1", 0.0, (;), [0.05, 0.3, 2.0]), ("pC0.3", -5.0, (;), [0.01, 0.05, 0.2]), ("t-1", 4.0, (;), [0.5, 2.0, 5.0]), ("t1", -5.0, (gls = -13.5,), [0.01, 0.05, 0.2]))
    line("vcdf", d, c, kw, yv, G.v_cdf(designs[d], yv, c; gkw(kw)...))
    line("vpdf", d, c, kw, yv, G.v_pdf(designs[d], yv, c; gkw(kw)...))
end
xa = [-3.0, -2.0, -1.0, 0.0, 1.0]
for (d, c, kw) in (("t0", 0.0, (;)), ("t1", -10.0, (start = :stationary,)), ("pC0.3", -5.0, (;)))
    line("adfcdf", d, c, kw, xa, G.adfmax_cdf(designs[d], xa, c; gkw(kw)...))
    line("adfpdf", d, c, kw, xa, G.adfmax_pdf(designs[d], xa, c; gkw(kw)...))
end
line("adficdf", "t0", 0.0, (;), [0.05], [G.adfmax_quantile(designs["t0"], 0.05; tol = 1e-12)])
line("env", "t-1", 0.0, (;), [-7.0, -15.0], G.power_envelope(designs["t-1"], [-7.0, -15.0]))
line("env", "t1", 0.0, (;), [-7.0], G.power_envelope(designs["t1"], [-7.0]))
line("lp", "pC0.3", 0.0, (;), [-5.0, -10.0], G.local_power(designs["pC0.3"], [-5.0, -10.0]))
line("region", "t-1", -7.0, (a = -3.5, b = -0.26), Float64[], [G.region_prob(designs["t-1"], (type = "line", a = -3.5, b = -0.26), -7.0)])
line("region", "pC0.3", -5.0, (x = -3.0,), Float64[], [G.region_prob(designs["pC0.3"], (type = "parabola", x = -3.0), -5.0)])
# ---------- the companions ----------
line("ptq", "pB0.1", 0.0, (cbar = -16.0,), [0.05], [G.pt_quantile(designs["pB0.1"], 0.05, -16.0)])
line("ptcdf", "pB0.5", 0.0, (cbar = -17.8,), [6.789, 10.0], G.pt_cdf(designs["pB0.5"], [6.789, 10.0], -17.8))
line("ptpdf", "pB0.5", -5.0, (cbar = -17.8,), [6.789, 10.0], G.pt_pdf(designs["pB0.5"], [6.789, 10.0], -17.8, -5.0))
line("ptenv", "pC0.5", 0.0, (;), [-10.0], G.pt_envelope(designs["pC0.5"], [-10.0]))
line("ptpow", "t1", 0.0, (cbar = -10.0,), [-5.0], [G.pt_power(designs["t1"], -10.0, -5.0)...])
line("lmq", "t1", 0.0, (;), [0.01, 0.05, 0.10], [G.lm_quantile(designs["t1"], p) for p in (0.01, 0.05, 0.10)])
line("lmq", "pC0.5", 0.0, (;), [0.01, 0.05, 0.10], [G.lm_quantile(designs["pC0.5"], p) for p in (0.01, 0.05, 0.10)])
line("lmp", "t1", 0.0, (;), [-5.0, -10.0], G.lm_power(designs["t1"], [-5.0, -10.0])[2])
line("lmcdf", "t2", -5.0, (;), [-4.0, -3.5, -3.0], G.lm_cdf(designs["t2"], [-4.0, -3.5, -3.0], -5.0))
line("mzacdf", "pB0.5", 0.0, (gls = -17.8,), [-23.659, -10.0], G.mza_cdf(designs["pB0.5"], [-23.659, -10.0]; gls = -17.8))
line("mzaq", "pB0.5", 0.0, (gls = -17.8,), [0.05], [G.mza_quantile(designs["pB0.5"], 0.05; gls = -17.8, y0 = -25.0, tol = 1e-12)])
line("msbcdf", "pB0.5", 0.0, (gls = -17.8,), [0.144, 0.2], G.msb_cdf(designs["pB0.5"], [0.144, 0.2]; gls = -17.8))
line("msbq", "pB0.5", 0.0, (gls = -17.8,), [0.05], [G.msb_quantile(designs["pB0.5"], 0.05; gls = -17.8)])
# ---------- the router: the merged interface with engine = auto, classical, general ----------
xr = [-3.0, -2.0, -1.0]
for (cas, c, kw) in (("DU", -5.0, (;)), ("DT", -300.0, (c0 = -300.0,)), ("SU", -10.0, (c0 = -10.0,)), ("DU", -5.0, (c0 = 2.0,)), ("DK", 0.0, (g0 = 1.0,)),
                     ("DG", -10.0, (;)), ("DK", -5.0, (c0 = 2.0,)))
    for eng in (:auto, :classical, :general)
        line("rcdf", cas, c, merge(kw, (engine = eng,)), xr, quiet(() -> DFLimit.ltratio_cdf(xr, cas, c; engine = eng, kw...)))
        line("rpdf", cas, c, merge(kw, (engine = eng,)), xr, quiet(() -> DFLimit.ltratio_pdf(xr, cas, c; engine = eng, kw...)))
    end
end
for (d, c, kw) in (("pC0.3", -5.0, (;)), ("t0", -5.0, (;)), ("t1", -10.0, (start = :stationary,)), ("t1", -10.0, (gls = -13.5,)), ("df", -5.0, (;)))
    line("rcdf", d, c, kw, xr, quiet(() -> DFLimit.ltratio_cdf(xr, designs[d], c; kw...)))
end
for (cas, c) in (("DT", -5.0), ("DK", 0.0), ("DU", 0.0))
    for eng in (:auto, :classical, :general)
        line("ricdf", cas, c, (engine = eng,), [0.05, 0.5], quiet(() -> DFLimit.ltratio_icdf([0.05, 0.5], cas, c; engine = eng, tol = 1e-12)))
    end
end
line("ricdf", "pC0.3", 0.0, (;), [0.05], [quiet(() -> DFLimit.ltratio_icdf(0.05, designs["pC0.3"], 0.0; tol = 1e-12))])
for (cas, c) in (("DT", -5.0), ("DU", -100.0))
    for eng in (:auto, :classical, :general)
        line("rmom", cas, c, (engine = eng,), Float64[], quiet(() -> DFLimit.ltratio_mom(4, cas, c; engine = eng)))
    end
end
for eng in (:auto, :classical, :general)
    line("rcoefcdf", "DU", 0.0, (engine = eng,), [-14.094], [quiet(() -> DFLimit.coef_cdf(-14.094, "DU", 0.0; engine = eng))])
    line("rlincdf", "DK", -5.0, (a = NaN, engine = eng), [0.1, 0.3], quiet(() -> DFLimit.ltratio_lincdf([0.1, 0.3], NaN, "DK", -5.0; engine = eng)))
    line("rlincdf", "DK", 0.0, (a = -3.5, engine = eng), [-0.6, -0.2639, 0.3], quiet(() -> DFLimit.ltratio_lincdf([-0.6, -0.2639, 0.3], -3.5, "DK", 0.0; engine = eng)))
    line("radfcdf", "SU", -10.0, (engine = eng,), [-3.0, -2.4], quiet(() -> DFLimit.adfmax_cdf([-3.0, -2.4], "SU", -10.0; engine = eng)))
    line("renv", "DK", 0.0, (engine = eng,), [-10.0], [quiet(() -> DFLimit.power_envelope(-10.0, "DK"; engine = eng))])
    line("rregion", "DK", -5.0, (a = -3.5, b = -0.2639, engine = eng), Float64[], [quiet(() -> DFLimit.ltratio_region((type = "line", a = -3.5, b = -0.2639), "DK", -5.0; engine = eng))])
end
line("renv", "pC0.5", 0.0, (;), [-10.0], [quiet(() -> DFLimit.power_envelope(-10.0, designs["pC0.5"]))])
# ---------- the repairs of November 2026 (DFLimit.jl 6e3109e to 30fd050) ----------
# the LM statistic: the edge t_* of its transform from the quadratic form of its V; the level breaks have the law of the
# linear trend at c = 3, 4; c = 5, 7 on the path of V alone (dfdet_test checks the refusal above c = 10)
designs["pA0.3"] = G.perron(:A, 0.3); designs["pA0.5"] = G.perron(:A, 0.5)
for d in ("t1", "t2", "pA0.3", "pB0.5", "pC0.3,0.7"), c in (-20.0, 2.0, 5.0)
    line("lmtstar", d, c, (;), Float64[], [G.tstar(G.model(G.lmprep(designs[d])), c)[1]])
end
xl = [-0.5, -0.2, -0.1, -0.05]
for d in ("t1", "pA0.3", "pA0.5"), c in (3.0, 4.0)
    line("lmcdf", d, c, (;), xl, G.lm_cdf(designs[d], xl, c))
end
for c in (2.5, 5.0, 7.0); line("lmcdf", "t1", c, (;), xl, G.lm_cdf(designs["t1"], xl, c)); end
# the point-optimal statistic just below the first c that it refuses (dfdet_test checks the refusals), at non-finite q
for (d, cb, c) in (("t1", -13.5, 3.6), ("t0", -7.0, 3.0), ("pB0.5", -17.8, 4.3), ("pC0.5", -20.0, 4.2))
    line("ptcdf", d, c, (cbar = cb,), [4.0, 10.0, 30.0], G.pt_cdf(designs[d], [4.0, 10.0, 30.0], cb, c))
end
line("ptcdf", "t1", 0.0, (cbar = -13.5,), [1e3, Inf, -Inf, NaN], G.pt_cdf(designs["t1"], [1e3, Inf, -Inf, NaN], -13.5, 0.0))
# the law of V off its support, MSB at NaN and Inf, the envelopes at NaN (by hand: the Julia functions do not give these
# in this form; at c = NaN they return the level)
line("vcdf", "t-1", -5.0, (;), [-0.1, 0.0, NaN, Inf, -Inf], [0.0, 0.0, NaN, 1.0, 0.0])
line("vpdf", "t-1", -5.0, (;), [-0.1, 0.0, NaN, Inf, -Inf], [0.0, 0.0, NaN, 0.0, 0.0])
line("msbcdf", "pB0.5", 0.0, (gls = -17.8,), [NaN, Inf, -1.0], [NaN, 1.0, 0.0])
line("env", "t-1", 0.0, (;), [NaN, -7.0], [NaN; G.power_envelope(designs["t-1"], [-7.0])])
line("ptenv", "pC0.5", 0.0, (;), [NaN, 1.0], [NaN, 0.05])
# the line integral with the explicit kernel for the designs of rank one or two (method = :line of the Julia engine, the
# default at c0 = 0), which the MATLAB programs do not have: their inversion on paths (opts.method 'paths', and 'auto')
# against its values, polynomial trends of order 2 to 4, the intercept with the first sine and cosine, the statistic that
# detrends first, zero and stationary start; opts.method = 'line' refused (DFDET_METHOD)
for d in ("t2", "t3", "t4", "fourier", "df"), c in (0.0, -5.0, -13.5, -20.0, -50.0, 2.0)
    D = designs[d]; Dm = d == "df" ? G.perron(:B, 0.5) : D
    m = G.dist_moments(Dm, c; p = 2)[1]; sd = sqrt(m[2] - m[1]^2); xs = m[1] .+ sd .* [-2.0, 0.0, 2.0]
    line("cdf", d, c, (method = :paths,), xs, G.dist_cdf(D, xs, c; method = :line))
    line("pdf", d, c, (method = :auto,), xs, G.dist_pdf(D, xs, c; method = :line))
end
for d in ("t2", "fourier"), c in (-5.0, -13.5)
    D = designs[d]
    m = G.dist_moments(D, c; p = 2, start = :stationary)[1]; sd = sqrt(m[2] - m[1]^2); xs = m[1] .+ sd .* [-2.0, 0.0, 2.0]
    line("cdf", d, c, (start = :stationary, method = :paths), xs, G.dist_cdf(D, xs, c; start = :stationary, method = :line))
end
line("err", "t2", 0.0, (fn = :cdf, method = :line, txt = Symbol("not_available_in_the_MATLAB_programs")), [-2.0], Float64[])
line("rcdf", "t2", -5.0, (method = :paths,), [-3.5, -2.5], DFLimit.ltratio_cdf([-3.5, -2.5], designs["t2"], -5.0; method = :line))
# the router: the linear form and the regions of no deterministic terms with a fixed start (the general engine; the
# classical programs dropped the start), ADFmax with a fixed start (the general engine)
for eng in (:auto, :general)
    line("rlincdf", "t-1", -5.0, (a = -3.5, g0 = 1.5, engine = eng), [-0.2639], [DFLimit.ltratio_lincdf(-0.2639, -3.5, designs["t-1"], -5.0; g0 = 1.5, engine = eng)])
    line("rregion", "t-1", -5.0, (a = -3.5, b = -0.2639, g0 = 1.5, engine = eng), Float64[], [DFLimit.ltratio_region((type = "line", a = -3.5, b = -0.2639), designs["t-1"], -5.0; g0 = 1.5, engine = eng)])
    line("rregion", "t-1", -5.0, (x = -2.0, g0 = 1.5, engine = eng), Float64[], [DFLimit.ltratio_region((type = "parabola", x = -2.0), designs["t-1"], -5.0; g0 = 1.5, engine = eng)])
    line("radfcdf", "DU", -5.0, (g0 = 1.5, engine = eng), [-3.0, -2.4, -1.5], G.adfmax_cdf(designs["t0"], [-3.0, -2.4, -1.5], -5.0; start = :fixed, gamma0 = 1.5))
    line("radfcdf", "t1", -5.0, (g0 = 1.0, engine = eng), [-3.0, -2.4, -1.5], G.adfmax_cdf(designs["t1"], [-3.0, -2.4, -1.5], -5.0; start = :fixed, gamma0 = 1.0))
end
close(io)
println("wrote ", OUT)
