# News

## Unreleased

**The general engine in C.** `c/dfl_det_*.c` is a port of the general engine of DFLimit.jl (`src/DFDet.jl` and
`src/DFDet_paper2.jl`, branch engine-merge): the limiting law of the Dickey-Fuller t-ratio for any piecewise smooth
deterministic design (polynomial trends of any order, breaks in level and slope at several dates, designs given by a
function or by the regressor matrix of a sample, the statistic that detrends first), with a zero, fixed or stationary
start, a local null value and GLS detrending, by the inversion formula on paths around the roots of the joint moment
generating function. The functions are those of the Julia engine: the distribution function, density and quantiles of
the t-ratio and of the coefficient statistic, moments, local power, the linear form U - aV and the law of V, ADFmax, the
power envelope, the regions below a line or a parabola, and the companions of the second paper (the point-optimal
statistic and its envelope for any design, the LM statistic, MZ_alpha, MSB). A design crosses the interface as plain
arrays (`dfl_design`), with constructors for the standard designs; the declarations are in `c/dflimit.h`. The engine is
self-contained like the classical one (no LAPACK: the small eigenvalue problems, solves and null spaces are written out),
multithreaded on the same pool, and keeps its prepared designs and outer paths in caches behind one lock. The stored
outer paths of cases DK, DU, DT at the unit root are compiled in (`c/dfl_det_tables.c`, written by `c/gen_det_tables.jl`).
On one thread a value takes 4 to 17 percent longer than in the Julia engine, on four threads 0.85 to 1.17 times its
time, on ten 0.78 to 1.0 times, and the outer paths are built 1.1 to 2 times faster (twelve designs and starts, warm
calls).

**One interface over the two engines.** `dfl_engine_*` (`c/dfl_det_route.c`, a port of `src/merged.jl`) takes a classical
case or a design and the engine (auto, classical, general); `dfl_autoengine` is the rule of the Julia package, with the
cells where the timings of the two C engines point elsewhere moved (`c/test/engine_table.c`, `c/test/engine_table_c.txt`):
the t-ratio at c0 = 0 beyond c = -20 and with c0 < 0 in cases DT (-30 <= c <= 0) and DU (-10 <= c <= 0) and most
quantiles at c0 = 0 go to the general engine, the moments with c0 != 0, the power envelope and the linear form of case DG
below c = 2 to the classical engine.

**In MATLAB.** The MEX file reaches the general engine through designs: `ltratio('design', ...)` makes one,
`ltratio('cdf', x, D, c, opts)` (and `'pdf'`, `'icdf'`, `'mom'`) takes the router, `ltratio('det', fn, D, ...)` the other
functions; `c/mex/det_check.m` checks a built MEX file.

**Tests.** `c/dfl_test_det` compares the general engine and the router with reference values of the Julia engine
(`c/test/det_ref.jl`, `c/test/det_ref.txt`): every constructor, one to three breaks, polynomial degrees 0 to 4, the three
starts, GLS detrending, a local null value, the tails and the edge of the support, every companion, and the three settings
of the router; it is built by `c/build.sh` and run by the workflow of the C engine.

**MATLAB.** The general engine of DFLimit.jl (`DFLimit.DFDet`: any deterministic terms, polynomial trends of any
order, breaks in level and slope at several known dates, smooth functions of t/n, by the inversion on paths around the
roots of the joint moment generating function) is ported to plain MATLAB: the constructors of the designs
(`dfdet_trend`, `dfdet_perron`, `dfdet_fdesign`, `dfdet_ltrend`, `dfdet_design_from_matrix`, `dfdet_detrend_first`,
`dfdet_determ`), the distribution function, density, quantiles and moments of the t-ratio (`dfdet_cdf`, `dfdet_pdf`,
`dfdet_icdf`, `dfdet_mom`, with a local null value, a zero, fixed or stationary start, GLS detrending and the
statistic that detrends first), the coefficient statistic, the linear form U - aV, the law of V, ADFmax, the envelope,
local power and regions (`dfdet_*`), and the companions `pt_*`, `lm_*`, `mza_*`, `msb_*`; the user-level routines of
the classical programs take `opts.engine` ('auto', the default, by the rule `ltratio_engine` set from the timings of
the two MATLAB engines; 'classical'; 'general') and a design in place of the case. `dfdet_test` checks the port
against 1113 values of the Julia engine (within 1.3e-11 in distribution functions and densities, 2.2e-11 in
quantiles); the values of the classical programs are unchanged where they compute.

**Repairs of the general engine (November 2026), in C and in MATLAB, as in DFLimit.jl 6e3109e to f824375.** The LM
statistic's edge t_* of its transform is estimated from the quadratic form of its own denominator (the estimate of the
orthogonal projection was up to 2.3 times too large, and above the unit root the LM functions failed, Perron's model A at
c = 3 to 5); above the unit root the LM statistic takes the path of its denominator alone, a hyperbola, in place of the
vertical line of the t-ratio (P(V <= y) was off by 6e-3 at c = 5 with a linear trend), and it refuses c > 10, where the
contours disagree. The point-optimal statistic refuses where the grid of its inversion does not resolve the transform
(far above the unit root: with a linear trend and cbar = -13.5 from c = 3.7, where the distribution function had come out
between -0.1 and 3.4 at c = 5; the error names the grid size that serves) and at |q| beyond 0.15/h (1000 with the default
grid); 0, 1 and NaN at -Inf, Inf and NaN. The U-shapes of the zero and the stationary start stop with an error beyond 1e6
nodes, as those of the fixed start did (the design (1 + r, 1{r < 0.4}), whose positive root runs off on the vertical
line, took gigabytes within seconds; it is now refused at every c). The refusal comes before any U-shape is built: the
outer path carries the roots at its nodes, and the count is checked at every node where the transform is not negligible
before the outer sum (`ushape_guard`, as in DFLimit.jl), so that such a design fails within 0.1 s and a few MB in C, where
the guard in the U-shapes refused after 1.4 s and 1.2 GB, and in MATLAB within the 0.8 s its outer path takes, where it
refused after 11 s. The checked computation (`verified`) falls back to the
variant at the default angle after any error of the steeper variant, except an exhausted memory. Values at c <= 0
unchanged (bit for bit in the reference values of both tests).

**The rules far from the unit root from c = -60 on (C and MATLAB, as in DFLimit.jl).** On -100 <= c < -60 the outer
contour of the general engine was the 45-degree hyperbola, on which the terms grow far from the unit root: the
distribution function was off by 1e-11 at c = -70 and by up to 1.4e-9 at c = -100 against the line integral (cases DU,
DT, SK, SU at c0 = 0, case DU at c0 = 2), the linear form of the envelope by up to 9e-10. The angle of `fixed_phi` (with
every start) and the caps 20/min(A, 40/|sigma|) of the segments of the U-shapes now apply from c = -60 on
(`DET_FAR_CMIN`), and the contour ends at the larger of 1e5 and 20 c^2 for every c: within 5e-13 there (case DK 3e-11, as
below c = -100), values elsewhere unchanged bit for bit, and on that range the C general engine at c0 = 0 takes 0.6 to 0.9
times the classical engine's time in cases DU, DT, SU (1.4 to 1.7 in case DK, as below c = -100; in DFLimit.jl the
general engine takes 0.2 to 0.95 times its former time there in cases DU, DT and 1.3 to 1.5 times in cases DK, SU).
The MATLAB rule (`ltratio_engine`) no longer keeps the
classical programs on that range; `dfl_test_det` and `dfdet_test` check it (the general engine within 1e-10 of the line
integral at c = -70, -85, -100).

**Non-finite arguments (C and MATLAB, as in DFLimit.jl).** Every public function of the general engine and of the router
gives at +-Inf and NaN in each argument either the limit or an error that names the argument. `dfl_det_coef_cdf_gp` did not
return at y = +-Inf (the full sum on all threads; now 1 and 0 without the sum, NaN at NaN); the searches from a start
value of +-Inf (`dfl_det_quantile`, `dfl_det_coef_quantile`, `dfl_det_lin_quantile`, `dfl_det_adfmax_quantile`,
`dfl_det_mza_quantile`) did not return or returned +-Inf, and a tolerance of +Inf stopped them after one step (now
errors; NaN keeps the default); the point-optimal statistic at c = +-Inf or cbar = -Inf and with umax = +Inf, and
`dfl_det_coef_cdf_gp` with umax = +Inf, did not return (now errors); the linear form at a = +-Inf gave NaN, through the
router's classical code 0, also at cbar = +-Inf in case DG (now errors; a = NaN is the law of V as before); the region
below a line gave NaN at a = +-Inf and the law of V at a = NaN, through the classical code 0 at b = Inf and at a = +-Inf
(now an error for the slope, and 1, 0, NaN at b = Inf, -Inf, NaN); the envelopes failed at c = -Inf (now 1, the limit);
GLS detrending at cbar = -Inf and MSB at a non-finite c gave errors that named another quantity (now named). The
sizes of the general engine are bounded: the outer hyperbola refuses more than 1e5 nodes (a factor hfac of its step far
below one, or not positive; in C the count overflowed an int) and the moments a grid of more than 1e6 points (c beyond
2.5e4). `dfl_det_local_power` checks its own level and c0, and the router's message of the level is that of the general
engine, "the level must lie in (0, 1)".

**Defects of the C interface (the R and Python wrappers).** `dfl_warning()` and `dfl_det_warning()` are cleared at the
start of every `dfl_engine_*` function and of `dfl_design_from_matrix` (a warning of an earlier call was read again), and
the new public functions `dfl_warning_clear()` and `dfl_det_warning_clear()` clear each; the router's quantile search
kept its hold on the general engine's warnings when it refused p outside (0, 1), after which that thread's warnings were
never cleared again (every return releases it now). ADFmax with a fixed start through the router computed the zero start of the classical ADFmax:
it takes the general engine; ADFmax with GLS detrending or c0 != 0 (dropped without notice) and a case code with gls
other than NaN are refused. The quantile functions of the general engine (`dfl_det_lin_quantile`,
`dfl_det_adfmax_quantile`, `dfl_det_pt_quantile`) and the levels of the power functions refuse p outside (0, 1) instead
of growing the bracket (`dfl_det_lin_quantile` at p = 1.5 had exhausted the memory). The law of V is 0 at b <= 0, the
envelopes of the general engine are NaN at c = NaN, as the classical envelope, and MSB is NaN at NaN.

**The MATLAB router.** The regions and the linear form with no deterministic terms and a fixed start take the general
engine (the classical programs have the zero start only, and the start was dropped: 0.285 for 0.975), ADFmax with a
fixed start too; `opts.start` or `opts.gls` with a case string, a zero start named with `g0 ~= 0`, `opts.g0` with the
cases of the linear form, the regions and the envelope, and `opts.c0` where there is no local null value are errors
(they were dropped without notice). With the folder of the MEX file ahead of `matlab/` on the path, the user-level
routines hand their calls with a design, or with a case string and `opts.engine`, without options of the MATLAB
engines to the C router (`ltratio_croute`): its rule and both engines in C; before, the classical engine was taken
throughout and the designs went to the general engine in MATLAB (a case string without `opts.engine` keeps the
classical programs, as `ltratio.m` does).
The MEX file's `route` operation has ADFmax, the linear form, the envelope and the regions, and gives the warning
`ltratio:brittle` where it takes the classical engine in a brittle regime.

## 0.5.0

Four errors are repaired, in the C engine (and with it the MEX file and the R and Python packages) and in the MATLAB
programs, as in DFLimit.jl 0.5.0. Values outside the regions named below are unchanged.

**Large x above the unit root.** In the bulk of the distribution at c = 7 and beyond (x in the hundreds) the inner
integrals of the line integral need a Gauss-Legendre rule whose node count grows like x^2. Beyond 4096 nodes the C
engine used a rule that it does not hold and crashed: in case DT at c = 7, x = 751 a segmentation fault after 20
seconds, which ends MATLAB with the MEX file, and R or Python with their packages. The MATLAB programs built every rule
by Golub-Welsch, with its full matrix, and passed 2 GB at the same point. The MATLAB programs now compute rules above
1024 nodes by Newton's method, as the C engine does, and both refuse a rule above 4096 nodes before anything is
evaluated: the automatic method takes the conditional integral at such a point (0.990823749760806 at the point above,
in 0.2 seconds with the C engine).

**A negative local null value near zero.** With c0 < 0 close to zero and x > 0, in every case and with every start, a
root of X(v) = c0 v^2 + x v + 1/2 lies far out, near x/|c0|, and the paths of the line integral's inner integrals grow
like (x/c0)^2; the C engine stored every node, and a call took more than 2 GB within two seconds (case SK at c = c0 =
-1e-3, x = 0.05 and at c = c0 = -1e-2, x = 1; the same in the MATLAB programs). The C engine now sums the paths as it
builds them, and both refuse a path of more than two million nodes before it is evaluated. The automatic method takes
the conditional integral at such a point, with its checks, and warns where such a value is flagged (`dfl_warning`; the
warning `ltratio:bound` in MATLAB); with the line integral requested, and in case DG, which has the line integral only,
an error names the arguments (the return code 10 + the bound in C, the identifier `ltratio:bound` in MATLAB). Below the
bound the footprint of such a point fell from 1 GB to 5 MB (case SK, c = c0 = -0.1, x = 1.5).

**A piece of the path that climbed.** Where the exponent of the line integral's integrand rises into a segment at one
of its ends, the piece of steepest descent from that end could climb by many orders of magnitude and cancel against the
next piece. With a fixed start (the line integral of the C engine: case DK by default, cases DU, DT with the line
integral requested) and c0 of -0.01 or -0.001 at x > 0 and c >= 0, the values were wrong, from errors of 1e-7 and
2.4e-3 that look plausible to 6e74: in case DK at c = 0, c0 = -0.01, g0 = 1, x = 1.5 the engine returned 6.03e74 for
0.960992, at c = 2 0.08440 for 0.08676. The same fault gave -3.5e86 and 1.2e8 for 0.48 with c0 = 1 at c = -20 in cases
DU, DT with the line integral requested, and, with a zero start, 1.1e13 (the MATLAB programs 1.7e13) for 0.68267 in
case DK at c = 0, c0 = -1e-4, x = 0. Such a piece is now cut back where it climbs above the larger end of the segment.

Of 6945 values of the C engine across its functions (the six cases and DG, the three starts, c from -50 to 10, c0 of
both signs, the three methods with and without `verified`, quantiles, moments, the coefficient statistic, ADFmax, the
regions, the linear form and the power envelope) 6932 are unchanged bit for bit; the others are four of the wrong values
above, seven that move by at most 3e-13 where the descent piece is now cut back, and two in case DK at c = 5,
c0 = -0.1, x = 14.17, whose paths now exceed the bound and which the conditional integral gives within 3e-15 of
before. At the repaired points the C engine agrees with DFLimit.jl 0.5.0 within 1e-13. Of 1398 values of the MATLAB
programs across the same functions all are unchanged bit for bit, and `ltratio_test` passes, also on the MEX file. At
the flagged points near c0 = 0 the conditional integral of the MATLAB programs differs from that of the C engine by up
to 4e-9 (the density of case SK at c = c0 = -1e-3, x = 0.05, which is off by 1.2e-3 in both). The C test program
`dfl_test_bounds` checks the repaired points. The MEX file gives the warning of `dfl_warning` and passes the refusals
of the line integral to the MATLAB programs as `dfx:bound:rule` and `dfx:bound:path`, so that `ltratio_dist` routes
them as it does those of the MATLAB engine.

**The linear form for small |a|.** `ltratio_lincdf` (`dfl_lincdf`) inverted the characteristic function of U - aV by
the formula of Gil-Pelaez on a grid that ends at u = 4000. At a = 0 (the law of U, whose density is infinite at the edge
b = -1/2 of its support) the transform decays only like u^(-1/2), and the values were off by up to 7e-2 within 1e-3 of
the edge, 1e-2 within 1e-2 of it and 1e-4 in the bulk, in cases DK and DG alike, as low as -1.5e-2 below the support
(where the value is 0) and above 1 by up to 5e-5 (2380 of 12541 values of a map over a, b and c lay outside [0, 1]).
The error stays at that size for |a| up to about 3e-3, falls to 1e-3 in case DK (4e-5 in case DG) at |a| = 0.3 and to
3e-8 at |a| = 1, and is larger where U and V are concentrated (1.3e-3 at c = -50 for small |a|; at c = -200 and -500
even for |a| of 3 to 10) and where U - aV is wide (up to 1.2e-2 in case DK at c = 5). P(V <= b) (a = NaN) was off by
8e-8 at c = -50, 4e-2 at c = -200, 8e-2 at c = -500 and 1.2e-5 at c = 5. The C engine and the MATLAB programs agreed
with DFLimit.jl to 1e-14 in these errors. Where the transform has not fallen below 5e-9 on the last 101 points of the
grid, or c > 0, the value now comes from a conditional integral: in case DK that of the region {U <= aV + b}
(`dfl_region`, `tr_region`), in case DG the conditional integral given R = U + 1/2, whose law is a gamma law of shape 1/2
under every tilt in V (`lc_x` in `c/dfl_gilpelaez.c`, `tr_lincond`), with outer panels that are halved where 24 and 12
Gauss-Legendre points disagree: to about 1e-10, within [0, 1] and 0 below the support (against the closed form at a = 0
within 3e-15, against a converged inversion within 3e-13). P(V <= b) takes a longer grid for c <= -50 and a finer one
for c > 2 (`dfl_vgrid`). `power_envelope` near c = 0, where its weight a is small, takes the conditional integral with
a bisection in x = b + 1/2 (in case DG it fell below alpha there; now alpha + O(c^2), 0.0500000034 at c = -1e-3; at
c = -1 0.0531924225 for 0.0531924279). `ltratio_lincdf` requires c >= -550, case DG c up to about 3.5 at cbar = -13.5
(less for cbar nearer 0) and a = NaN c up to about 7, with an error beyond; b = -Inf, Inf give 0, 1 in the MATLAB
programs too. `ltratio_region` with v0 > 0 changes with P(V < v0) at c <= -50 and c > 2. In the C engine and the
MATLAB programs the values agree with DFLimit.jl within 6e-14 (`dfl_test_app`, which now includes such points, and the
MATLAB programs on the same points within 2e-14). Of 6849 values of the C engine across its functions 6764 are
unchanged bit for bit; the 85 others are values of the linear form and of the envelope, of which 84 are closer to a
reference (the closed form at a = 0, an inversion on 1.2e6 nodes up to u = 1e5 elsewhere) than before, within 3e-13
of it, and one, 1 - 2.2e-14 before, moves by 1.5e-14. The MATLAB programs other than these four are not touched, and
`ltratio_test` (now 268 checks) passes, also on the MEX file.

## 0.4.1

Documentation only; the programs and the values are those of 0.4.0. The README gives the paper's address on SSRN;
`CITATION.cff` names the paper as the suggested citation, and `.zenodo.json` describes the package for its archive on
Zenodo; the pointers to the reproduction of the paper name the paper's replication package; the workflows run again on a
push that changes `c/`.

## 0.4.0

The first release as a repository of its own. The MATLAB programs, their manual and validation material, the C
engine with its tests and build script, the MEX gateway and the workflows that test the C engine and build the MEX
file were developed in the paper's repository (DF-Limit, the directories `public/` and `c/`); the code here is that of
DF-Limit commit b484d4f, with the MATLAB programs in `matlab/` (formerly `public/`). The values are unchanged.

The changes are those of the move. The paths that named `public/` name `matlab/`; `dfx_test` reads its reference
values from `matlab/validation/`; the reference scripts `c/test/*_ref.jl` and `c/gen_tables.jl` load the Julia engine
from the package DFLimit.jl (the environment `c/test/Project.toml`); the reproduction of the paper (`ltratio_paper`,
`mackinnon94_p`) stays in the paper's repository; the logs of earlier validation runs stay in its history. The workflows
run when started by hand and on a tag `v*`, no longer on every push. `Contents.m` and the manual state the version.
The R package dflimit and the Python package dflimit-py now copy the C engine from this repository.
