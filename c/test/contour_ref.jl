# reference values of the contour integral representation (a fixed start g0 != 0 in cases DU, DT; DFExact.tr_contour) for
# the C engine test c/test/contour_ref.txt (columns: model c c0 g0 x F f err_F err_f), read by dfl_test_fixed; the same
# values in matlab/ltratio_test.m.  The points: the grid of the experiment (c in {-10, -5, -2, 1}, c0 in {0, c}), two roots
# of X(v) with c0 > 0, c0 < 0 far from c, near the double root x = -sqrt(2 c0), x = 0, c = 0 and c = 2, c = 1.9 in the upper
# tail, far from the unit root (c = -94.4, -70.5), large noncentralities ((c g0)^2/2 = 364 at c = -18, g0 = 1.5, and 7200 at
# c = -40, g0 = 3 in case DT, the corner of the region) and x -> 0- at c0 = 0 (the root of X at 1/(2|x|) = 5e7).
#   julia --project=c/test -t auto c/test/contour_ref.jl c/test/contour_ref.txt
using Printf
using DFLimit; const DFExact = DFLimit.DFExact; const D = DFExact
const PTS = [  # model, c, c0, g0, x
    (2, -5.0, 0.0, 1 / sqrt(10), -2.86), (2, -5.0, 0.0, 1 / sqrt(10), 1.0), (2, -5.0, 0.0, 1.5, -2.0), (2, -5.0, 0.0, 1.0, 0.0),
    (2, -2.0, -2.0, 1.5, 0.0), (2, -10.0, 0.0, 3 / sqrt(20), -4.0), (2, -10.0, -10.0, 3 / sqrt(20), -1.5), (2, 1.0, 1.0, 2.0, -2.86),
    (2, 1.0, 1.0, 0.5, -1.5), (2, 1.9, 0.0, 1.298, 3.304), (2, -94.41, 0.0, 0.051, -7.929), (2, -17.81, 2.19, 0.249, -4.144),
    (2, -15.66, -2.02, 0.557, -1.943), (2, 0.0, 0.0, 1.0, -1.5), (2, 2.0, 0.0, 1.0, 0.5), (2, -18.0, 0.0, 1.5, -9.0806812422194216),
    (3, -5.0, 0.0, 2.0, -3.0), (3, -10.0, -10.0, 3 / sqrt(20), -1.5), (3, -70.45, -4.25, 0.268, -4.3), (3, -7.53, 1.84, 0.9, -1.92),
    (3, 1.26, 1.95, -1.582, -1.901), (3, 1.0, 0.0, 2.0, -2.86), (3, -2.0, -2.0, 1.5, 1.0), (3, -52.54, 1.56, 0.172, -3.734),
    (3, 2.0, 0.0, 1.0, -2.0), (3, 0.0, 0.0, 3.0, -2.5), (3, -2.0, 0.0, 1.5, -4.0), (3, -20.0, 0.0, 1.5, -4.5),
    (2, -5.0, 0.0, 1.0, -1e-8), (3, -40.0, 0.0, 3.0, -13.0)]
if abspath(PROGRAM_FILE) == @__FILE__
    out = open(ARGS[1], "w")
    for (model, c, c0, g0, x) in PTS
        # (the checked computation, verified = true, with its error estimates; one evaluation, by the rule of one
        # evaluation since November 2026, is compared with these values at a looser tolerance)
        F, iF = D.tr_contour([x], model, c; c0, g0, pdf = false, warn = false, verified = true)
        f, i_f = D.tr_contour([x], model, c; c0, g0, pdf = true, warn = false, verified = true)
        (iF.status[1] == 0 && i_f.status[1] == 0) || @warn "flagged reference point" model c c0 g0 x
        @printf(out, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.3e %.3e\n", model, c, c0, g0, x, F[1], f[1], iF.err[1], i_f.err[1])
    end
    close(out)
end
