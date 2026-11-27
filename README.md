# dflimit-matlab

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24404382.svg)](https://doi.org/10.5281/zenodo.24404382)

The exact limiting distribution of the Dickey-Fuller t-ratio in MATLAB, and the C engine of the R and Python packages.

The programs compute the distribution function, density, quantiles and moments of the limits of the Dickey-Fuller
t-ratios under local alternatives, phi = 1 + c/n, by the exact representations of two papers. The lead paper,

> Bao, Y., Broda, S. A. and Kan, R., "Exact Limiting Distributions of Unit Root t-Ratios", working paper, 2026,
> available at <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7745921>,

gives the inversion formula on paths around the singularities of the joint moment generating function and with it the
limit for any piecewise smooth deterministic terms (the general engine, `dfdet_*` and `c/dfl_det_*.c`); its companion,

> Bao, Y., Broda, S. A. and Kan, R., "Low-Rank Representations of the Classical Dickey-Fuller Limits, with an Extension
> to Predictive Regressions", working paper, 2026, available at <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619>,

gives the representations of the classical cases, in which the numerator has rank one or two (the line integral and the
conditional integral), and the probabilities of regions of the plane of the numerator and the denominator.

The programs cover the regression without deterministic terms, with an intercept and with an intercept and a trend; a
zero, a fixed or a stationary initial observation; the t-ratio centred at the unit root or at a local null value c0; and
the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996). They also have the ADFmax test of Leybourne (1995),
the coefficient statistic, the probabilities of regions {U <= g(V)} of the limit experiment (point-optimal and
weighted-average-power tests), the distribution of the linear form U - aV and the asymptotic power envelope, and, by the
general engine, the same limits for any deterministic terms ("The general engine in MATLAB" below).

```matlab
addpath /path/to/dflimit-matlab/matlab
ltratio_icdf([0.01 0.05 0.10], 'DU')         % asymptotic critical values of tau_mu
ltratio_cdf(-2.86157, 'DU', -5)               % limiting power of the 5% test at c = -5
ltratio_cdf(-2.86157, 'SU', -5)               % the same with a stationary start
ltratio_cdf(0, 'DK', 2, struct('c0', 2))      % t-ratio centred at the true value
ltratio_cdf(-2, 'DU', -5, struct('g0', 1.5))  % fixed start
ltratio_mom(4, 'DT')                          % first four moments of tau_tau
```

## Two engines

The repository holds two implementations of the same functions with the same arguments.

The MATLAB programs in `matlab/` are plain MATLAB (R2019b or later; they also run, more slowly, on GNU Octave). They
take 0.007 to 0.13 seconds per value on one thread by the default (up to 2 seconds with a representation forced, the
conditional integral with a fixed start; `c/mex/README.md`, November 27, 2026) and are the readable form of the method: `Contents.m` lists them and
the manual `matlab/doc/ltratio_manual.pdf` describes the user-level routines, the engines (the line integral `dfx_*`,
the conditional integral `tr_limdist`, the contour integral `tr_contour`), accuracy and cost.

The C engine in `c/` is a port of the Julia engine `DFExact.jl` of the package DFLimit.jl, multithreaded and 5 to 90
times faster (the same representation, ten threads against the MATLAB programs' one; `c/mex/README.md`). In MATLAB it is reached through the MEX file `ltratio_mex`, built from `c/` with the gateway
`c/mex/ltratio_mex.c`. The folder of the MEX file has `ltratio.m`, which offers the functions of the MATLAB programs
with the same arguments (`ltratio('cdf', -2.86157, 'DU', -5)`), and drop-in replacements of the engines of the MATLAB
programs (`dfx_core.m`, `tr_limdist.m`, `tr_contour.m`, `tr_limmom.m`, `tr_adfmax.m`, `tr_region.m`,
`tr_gilpelaez.m`): with that folder ahead of `matlab/` on the path, `ltratio_cdf` and the other user-level routines run
unchanged on the C engine. The C engine also has the line integral with a fixed start, which the MATLAB programs lack
(they take the conditional integral there), and its line integral takes the inner integrals on fixed node layouts where
they serve (c <= 0 and c0 = 0, as the Julia engine: `c/dfl_batch.c`), 15 to 20 times faster than the adaptive rules in
case DK on one thread (8 to 13 on four) and 1.9 to 3.2 times in the other cases (1.7 to 2.9 on four; F and f at five
points, warm calls on mains power, November 27, 2026); the MATLAB programs keep the adaptive rules there (they are the readable
form, and their values agree with the layouts' to about 1e-12).

Use the MATLAB programs when nothing but MATLAB may be installed or when the code is to be read, and the MEX file for
speed: power curves, confidence belts, anything with many values. Both give the same values to well below 1e-5.

`c/` also holds the general engine (`c/dfl_det_*.c`), a port of the general engine of DFLimit.jl (`DFDet.jl`): any
piecewise smooth deterministic design, given as the break points and the coefficient matrices of its pieces
(`dfl_design`, with constructors for polynomial trends, Perron's models with several breaks, designs given by a function
or by the regressor matrix of a sample, and the statistic that detrends first), with a zero, fixed or stationary start,
a local null value and GLS detrending, by the inversion formula on paths around the roots of the joint moment generating
function; and the companions of the lead paper (the point-optimal statistic and its envelope,
the LM statistic, MZ_alpha, MSB). For the designs whose numerator has rank one or two (polynomial trends of any order, an
intercept with sines and cosines, the statistic that detrends first for a continuous design that does not vanish at the
origin; zero or stationary start, least squares detrending) the distribution function, the density and the quantiles of
the t-ratio also come from the line integral with the explicit kernel: the classical line integral with P(t) and the two
roots supplied by the general engine's transform (`c/dfl_det_line.c`), which `dfl_det_opts.method` chooses
(`DFL_DET_AUTO`: the line integral at c0 = 0 for c >= -20, except just above the unit root with deterministic terms, 0 <
c <= 0.5, where the paths are 1.5 to 2.3 times faster, and with c0 > 0 at the points x < 0 where the outer contour of the
paths is the vertical line, where the paths do not resolve the two roots of the parabola of such a point and their one
evaluation was off by up to 1e-3 in the lower tail (the quantiles keep the paths there); the paths otherwise; `DFL_DET_LINE`, `DFL_DET_PATHS`). Measured on one thread and on
four (`c/test/line_timing.c` on mains power, November 27, 2026, and interleaved timings, November 26): below the unit root the line
integral is 1.0 to 1.6 times faster per value for -20 <= c <= 0 (at c = -30 and -50 the paths up to 2 times faster);
above it the paths are 1.3 to 2.2 times faster up to c = 2 for some designs with deterministic terms (at c = 2 1.4 to 1.6
times per value and 1.8 to 2.2 per quantile), but below 1.5 for others beyond c = 0.5 (the Julia rule, whose costs differ, takes the paths up
to c = 1.5 with one regressor and c = 2 with more). `dfl_engine_*` (`c/dfl_det_route.c`) is one interface over the two C engines, with
the automatic choice of the Julia package (`dfl_autoengine`), re-derived on November 26, 2026 from interleaved timings of
the two C engines (a cell moved where the other engine was at least 1.5 times faster) and differing from the Julia rule
where the C costs differ (the comments of `dfl_autoengine` and NEWS list the cells). It agrees with the Julia engine to 1e-9 in distribution
functions and densities and to 1e-7 in quantiles (`c/dfl_test_det`); in MATLAB the MEX file reaches it through designs
(`ltratio('design', ...)`, `c/mex/README.md`). Above the unit root the general engine computes up to c = 16 with
deterministic terms (10 with them and a fixed start) and c = 20 without them (`DFL_DET_CMAX`, `DFL_DET_CMAX_FIXED`,
`DFL_DET_CMAX0`, `dfl_det_cmax`): beyond, one root of its
transform is of the order 2c e^{-2c} and the others are no longer resolved in double precision, and every function of the
general engine stops with an error that names the statistic, while `dfl_autoengine` takes the classical engine. Far
above the unit root, for c > 5.5, the inversion on paths by one evaluation may be off by more than 1e-5 in the upper tail
for a design beyond the classical cases and for a fixed start with deterministic terms (Perron's model C at c = 7 by 1.5e-5
in F, a linear trend with a fixed start at c = 7 by 8.7e-5), and the general engine warns there (`dfl_det_warning`; the
MATLAB programs `dfdet:brittle`; `verified` checks it).

Above the unit root the t-ratio with c0 != c has its upper tail at x of the order (c - c0) e^c/(2c). There the classical
engine takes the conditional integral (beyond the line integral's bound, from x ~ 700, in case DK the checked value
whatever its flag; for c > 4 in the three cases also at x > 300, and in case DK with c0 > c at x < -20, `far_tail` in
`c/dfl_api.c` and `matlab/ltratio_dist.m`), within 1.3e-8 of an independent computation in case DK and 5e-7 of the large-c
law in cases DU, DT up to c = 14, and refuses these points above c = 14 with an error (before November 23, 2026 case DK
returned one evaluation of the conditional integral there, values outside [0, 1] among them, and cases DU, DT the line
integral, which exhausted the memory at c = 16). Every distribution function of both engines ends with a check of its values: outside
[0, 1] by more than 1e-8 an error that names the call (`dfl_cdf_range_check`, `ltratio_rangecheck`).

## The general engine in MATLAB

The folder `matlab/` also holds a port of the general engine of DFLimit.jl (`DFLimit.DFDet`, the engine of the
lead paper on general deterministic terms) to plain MATLAB, the routines `dfdet_*`: the limiting t-ratio, the
coefficient statistic, the linear form U - aV, the denominator V and ADFmax for any deterministic terms, polynomial
trends of any order, breaks in level and slope at several known dates, smooth functions of t/n, with a zero, fixed or
stationary start, GLS detrending, a local null value and the statistic that detrends first; and the companions of
that paper, the point-optimal statistic and the Gaussian power envelope for any design (`pt_*`), the LM statistic
(`lm_*`; for trend functions whose span contains the constant and a linear trend, an error otherwise), MZ_alpha and MSB (`mza_*`, `msb_*`). A design is built by a constructor (`dfdet_trend`, `dfdet_perron`,
`dfdet_fdesign`, `dfdet_ltrend`, `dfdet_design_from_matrix`, `dfdet_detrend_first`, `dfdet_determ`); the joint moment
generating function of the numerator and the denominator is a small determinant, its roots are followed along an
outer contour, and the distribution function and the density are the inversion on paths around the roots, a double
integral of functions in closed form (`dfdet_contour`).

```matlab
D = dfdet_perron('C', 0.3);                       % a break in level and slope at 0.3
dfdet_icdf([0.01 0.05 0.10], D)                   % critical values
dfdet_cdf(-4.2, D, -10)                           % local power at c = -10
dfdet_cdf(-2.9, dfdet_trend(1), -5, struct('gls', -13.5))      % DF-GLS
dfdet_icdf(0.05, dfdet_detrend_first(dfdet_perron('B', 0.5)))  % a slope break, detrended first
ltratio_cdf(-2.86157, 'DU', -5, struct('engine', 'general'))   % a classical case by the general engine
```

The user-level routines of the classical programs (`ltratio_cdf`, `ltratio_pdf`, `ltratio_icdf`, `ltratio_mom`,
`coef_cdf`, `coef_icdf`, `adfmax_*`, `ltratio_lincdf`, `power_envelope`, `ltratio_region`) take `opts.engine`:
`'classical'` (the programs above), `'general'`, or `'auto'` (the default), which takes in each cell the engine of
`ltratio_engine`, the rule of the Julia package set by the timings of the two MATLAB engines
(`matlab/validation/engine_table_matlab.txt`), with the cells that the C rule moved on November 26, 2026 moved as in C; they also take a design in place of the case, a classical design being
its case and any other going to the general engine. With the folder of the MEX file ahead of `matlab/` on the path
these routines hand a call with a design or a case string, without options of the MATLAB
engines (at most `c0`, `g0`, `start`, `gls`, `engine` `'auto'` or `'general'`, `method` `'auto'`, `verified` false,
and `tol`, `cbar`, `alpha` where the routine has them) to the C router (`dfl_engine_*`, `ltratio_croute`): the choice of engine is then that of the C rule, set by the timings
of the two C engines, and both engines compute in C, with the same checks of the arguments and the warnings of the C
engine (`opts.engine = 'auto'` is the default, for a case string as for a design; where the rule takes the general
engine for a case string the two engines agree to about 1e-9), and with `engine = 'classical'`, a representation named, `verified = true` or further options they
compute as before, the classical programs on the C engine through the drop-ins and the general engine in MATLAB. `dfdet_test` compares the general engine and the router with reference values of the Julia engine
(`matlab/validation/det_ref.txt`, written by `matlab/validation/det_ref.jl`).

The general engine in MATLAB has the range of the C engine (`dfdet_cmax`: c up to 16 with deterministic terms, 10
with them and a fixed start, 20 without them; beyond, the routines refuse with an error that names the statistic, and `ltratio_engine` takes the classical
programs).

The MATLAB programs have the inversion on paths only: `dfdet_cdf`, `dfdet_pdf` and `dfdet_icdf` take `opts.method`
`'auto'` or `'paths'` and refuse `'line'`, the line integral with the explicit kernel of the C engine and of the Julia,
R and Python packages, with a message that points to the MEX file (`dfdet_method`). Their line integral (`dfx_core`, with
the adaptive rules) is slower than their inversion on paths in the rank-two cases DU, DT and SU at c0 = 0 even with its
closed forms of P(t) and the roots: 3.0 to 3.9 times at one point and 3.0 to 3.6 times at five points in one
call, whose paths are shared (c = 0 and -10, warm calls on one thread on mains power, November 27, 2026; 1.04 to 1.64
times before the thinner rules of the paths), so that it would not be faster for any design. With the folder of the MEX file ahead on the path, `ltratio_cdf(x, D, c, struct('method', 'line'))` (and
`ltratio_pdf`, `ltratio_icdf`, `ltratio('cdf', x, D, c, struct('method', 'line'))`) take the C engine's line integral,
`'paths'` the C engine's paths, and `'auto'` the C rule.

## The extensions of the limit

Three extensions of the limit of the t-ratio, ported from the Julia package (`src/DFDet_vol.jl`, `DFDet_cov.jl`,
`DFDet_seas.jl` of DFLimit.jl) to the C engine (`c/dfl_det_ext.c`, the volatility transform in `c/dfl_det_mgf.c`, its
paths in `c/dfl_det_path.c`) and to the plain MATLAB programs, with the same names and arguments:

- **Volatility and persistence regimes.** `volpath(breaks, omega)` makes a path: the standard deviation of the
  innovations piecewise constant between known dates (`breaks`, strictly increasing in (0, 1); `omega` of any positive
  scale, only the ratios matter). With `opts.vol = volpath(...)`, `dfdet_cdf`, `dfdet_pdf`, `dfdet_icdf` and
  `ltratio_cdf`, `ltratio_pdf`, `ltratio_icdf` (a case string or a design) give the limit of the t-ratio when the
  process is dJ = c(r) J dr + omega(r) dW, with `c` a number or a vector of one value per regime (changing
  persistence). The transform is the general engine's determinant with a diagonal tilt at the regime boundaries and the
  covariance of the functionals under the regimes; the outer contour, the inner integrals, the rules of one evaluation,
  the checked computation and the quantile search are the engine's. Zero start, least squares detrending, the one-step
  regression, every `c` in [-480, 0], `opts.method` `'auto'` or `'paths'`; one regime with a scalar `c` is the engine
  itself.
- **A normal component.** `opts.rho`: the limit rho tau_c0 + sqrt(1 - rho^2) Z, Z standard normal and independent of
  the t-ratio's limit tau_c0, which is the limit of the covariate-augmented t-ratio of Hansen (1995) at
  `c0 = cadf_c0(c, rho, R)` and of the t-ratio of the slope in a predictive regression with a regressor local to unity
  at `c0 = predictive_c0(c, rho, b)`. `rho` in [-1, 1] (1 the t-ratio itself, 0 the normal law, a negative `rho` by the
  reflection); `opts.method` `'auto'`, `'lines'` (the inversion on two lines with the normal component in the kernel, a
  zero or stationary start, c <= 0) or `'mixture'` (the normal average of the t-ratio's distribution function, graded
  toward its kink when c0 > 0); `'auto'` takes the lines for 0.1 <= rho^2 <= 0.99 with a zero or stationary start at
  c <= 0 and the mixture otherwise. Least squares detrending and the one-step regression; not together with `vol`.
- **The first seasonal t-ratio** of the regression of Hylleberg, Engle, Granger and Yoo at a harmonic frequency (pi/2
  for quarterly data, the t-ratio of pi_3; each pair of harmonic frequencies of monthly data), the one of Dickey-Fuller
  shape: `hegy_t3_cdf`, `hegy_t3_pdf`, `hegy_t3_quantile` for any design (`dfdet_trend(-1)` no seasonal terms,
  `dfdet_trend(0)` seasonal dummies, `dfdet_trend(1)` seasonal dummies and seasonal trends, breaks), with `opts.c0` and
  `opts.vol`, c <= 0. The moment generating function of its numerator and denominator has poles in its first argument,
  and the inner integral of the inversion on paths is the sum of the residues at the roots on the engine's outer
  contour. The second t-ratio and the F statistic of the regression are not part of this release (nor of the Julia
  package's).

The costs of the extensions in the C engine (warm calls on four threads of an M4 laptop on mains power, November 27,
2026, on the cells of the documentation pages of DFLimit.jl; F and f are two calls): with `vol`, F and f at one point
3.5 to 9.6 ms and a 5% quantile 9 to 40 ms (10 to 52 ms with the outer contour built; no deterministic terms, an
intercept, a linear trend with the variance tripling at mid-sample, a break in level and slope with the variance falling);
with `rho` (rho^2 = 0.5 and 0.9; no deterministic terms, an intercept, a linear trend, a level break), F and f at one
point by the lines 28 to 80 ms and by the mixture 70 to 530 ms at c0 = 0 (6.4 s with an intercept at c = -5, c0 = 5),
a 5% quantile 69 to 381 ms; `hegy_t3_*` (single-threaded) F and f 0.09 to 0.15 ms on a cached contour and a 5% quantile
0.2 to 0.7 ms. On one thread these times are 2 to 4 times longer (`hegy_t3_*` the same).

```matlab
V = volpath(0.3, [1 3]);                                        % the standard deviation triples at 30% of the sample
dfdet_cdf(-2.86, dfdet_trend(0), 0, struct('vol', V))
dfdet_icdf(0.05, dfdet_trend(1), [-3 0 -10], struct('vol', volpath([0.2 0.7], [2 1 0.5])))
ltratio_cdf(1.645, 'DU', 0, struct('rho', -0.95, 'c0', predictive_c0(0, -0.95)))   % the predictive regression
dfdet_cdf(-2.5, dfdet_trend(1), -5, struct('rho', sqrt(0.6), 'c0', cadf_c0(-5, sqrt(0.6), sqrt(0.4))))
hegy_t3_quantile(0.05, dfdet_trend(0))                          % the 5% critical value with seasonal dummies
```

The MEX file takes them too: `ltratio('cdf', x, cas, c, opts)` (and `'pdf'`, `'icdf'`) with `opts.vol` (a structure
with the fields `breaks` and `omega`, as `volpath` makes it) and `opts.rho`, `ltratio('det', 't3cdf', D, x, c, opts)`
(and `'t3pdf'`, `'t3icdf'`); with its folder ahead on the path the user-level routines hand such a call to the C
router. `ext_test` checks the MATLAB programs against the reference values of the Julia engine
(`matlab/validation/ext_ref.txt`) and, where the MEX file is found, against the C engine; `c/dfl_test_ext.c` checks the
C engine against the same file.

## The rules of one evaluation

One evaluation (the default, `verified` false) takes thinner rules than before November 2026 where they hold, in C, in
the MEX file and in the MATLAB programs alike, with the rules, regions and constants of DFLimit.jl (`PATH_RULES`,
`path_region`, `ctr_rule`, `cond_rule`; its README, "The rules of one evaluation", has the sweeps that set them).

The general engine's inversion on paths (`dfl_det_dist`, `dfl_det_quantile`, and the paths that the C line integral
leaves to them; `dfdet_cdf`, `dfdet_pdf`, `dfdet_icdf`): the outer hyperbola with 1.5 times the step and the U-shapes
with 16 nodes per panel and 32 on the tails of their rays for 24 and 40 (32 and 48 for 48 and 64 with a fixed start),
by regions of the start, the deterministic terms, c and c0 (`dfl_det_path_region`, `dfdet_rules`):

| region | rule |
|---|---|
| c > -1, c0 < 0, zero or stationary start | the thin inner rule on the present hyperbola (`inner`) |
| c > -1, c0 < 0, fixed start, deterministic terms | the thin inner rule on the present hyperbola with half its step (`fine`) |
| c > -1, c0 < 0, fixed start, no deterministic terms | the present inner rule on the present hyperbola, the step min(1/2, 2.5/(g0 \|c0\|)) times its own where it is flatter than 70 degrees (`finestep`; the present step was off by up to 0.86 in the upper tail at c = 0.5, g0 = 1, c0 = -20) |
| c > -1, c0 = 0, fixed start, no deterministic terms | `inner` |
| c > 2 with a fixed start, c > 5; c < -60 with c0 other than 0 and c | the present rule (`present`) |
| elsewhere | the thin rule (`thin`) |

`dfl_det_opts.rule` (`DFL_RULE_AUTO` and the six rules), `opts.rule` of `dfdet_cdf` and of `ltratio('det', ...)`
name a rule (`present` the rule of before, which `verified` takes); `hfac`, `n`, `ntail` override it (`hfac` is now
NaN, by the rule, unless given). The guard on the U-shapes counts the nodes of the present rule.

The classical contour integral (a fixed start in cases DU, DT; `dfl_ctr_rule`, `ltratio_rules('contour')`): 32 nodes
per panel and ray for 64 and the outer step 1.5 times as long where they hold (32 on the present step above the unit
root, 64 on the longer step with c0 > 0 below c = -2.5 and everywhere below c = -10, the present rule where |c g0| > 20);
`dfl_contour_rule` and `tr_contour` take the rule `present`, `n1` and `hfac`. The classical conditional integral
(`dfl_cond_rule`, `ltratio_rules('conditional')`): 48 nodes for 96 with a stationary start for c <= -2 and with a zero
start in cases DU, DT for c >= -5 and in case DK for c <= -13.5, 64 with a stationary start for c > -2, 96 elsewhere,
with a fixed start and above c = 2; `Nv` and `Nl` override it (in C they are 0, the rule, by default).

The rule of the engines: the fixed start of cases DU, DT goes to the general engine with c0 > 0 below c = -20 (below
c = -10 where |c g0| > 20) and with any c0 below c = -20 where |c g0| > 60, where the contour integral's one evaluation
was off by up to 2.2e-2 of the peak of the density and 3.9e-5 in F (`dfl_autoengine`, `ltratio_engine`; the warning
`ltratio:brittle`, `dfl_brittle`, names that regime where the classical engine computes); case DK with a fixed start
and c0 < 0 near the unit root (c > -1) takes the classical engine in C (the MATLAB rule took it there already). The
placement of the U-shapes of a fixed start counts the growth of the v-kernels of the upper tail at every c (in case DK
near the unit root one evaluation was off by up to 1e39 in the far upper tail). The quantile searches start from the
moments from jets of Bao and Kan where they hold (`dfl_mom_start`, `mom_start` in `ltratio_icdf.m`), and where the rule
gives a Newton step's point to the general engine its paths return F and f from one pass (`dfl_det_dist_both`,
`dfdet_dist(..., 'both')`), the values of the two evaluations bit for bit; `DFL_MOM_START`, `DFL_DIST2`,
`ltratio_switch` and `ltratio('switch', ...)` restore the computation of before.

Where C and Julia take the same rule their values agree to 3.7e-14 (paths), 2.3e-14 (contour integral) and 7.1e-13
(conditional integral), the quantiles to 5.6e-11 (the search tolerance is 1e-10; `c/test/det_ref.txt`); the MATLAB
programs to 2e-14, 1.6e-14 and 2.1e-12, the quantiles to 1.9e-11 (`matlab/validation/det_ref.txt`). Warm time per call
new/old in C (minima of seven calls, ten threads): the general engine's paths 0.47 to 0.71 (1.4 to 1.7 in the cells of
the finer steps), the contour integral 0.39 to 0.69, the conditional integral 0.55 to 0.84, the router's quantiles 0.26
to 0.31 with c0 != 0 (0.80 to 0.98 where the classical search computes). In the MATLAB programs: warm calls (minima of
five after a first, one MATLAB session each): the paths 0.48 to 0.78 (1.6 in the cells of the finer steps), the contour
integral 0.44 to 0.83, the conditional integral 0.87 to 0.98, the quantiles of `ltratio_icdf` 0.34 to 0.48 where the
general engine computes some of their points, 1.01 with the classical engine.

## The devices of Bao and Kan

The half-plane statistics of the six classical cases and the moments of the t-ratio are computed by default with the
numerical devices of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares
Estimator in First-Order Autoregressions: O(1) Algorithms", as in the Julia package (its module `DFKan`): in C
`c/dfl_kan.c` and `c/dfl_det_kan.c`, in MATLAB `kan_core`, `dfdet_kan` and `kan_route`. The joint characteristic function
of the limit pair comes in closed forms that are free of cancellation on both sides of the unit root, with the branch of
the square root fixed at each point from factors whose arguments are known, and in Taylor-mode arithmetic (jets), so that
its derivatives come out of the same assembly. A statistic that rejects on a half-plane is one inversion along a ray, by
adaptive quadrature over decades (the algorithm of QuadGK.jl) with the tail beyond S = 400/|w| completed by repeated
integration by parts from a jet of order eight; the moments of the t-ratio are one-dimensional integrals of the j-th
derivative of the moment generating function, which one jet of order j gives at each node.

| function (MATLAB; C router) | the default | warm time in C against `method = engine` |
|---|---|---|
| `coef_cdf`, `coef_icdf`, `coef_pdf`, `coef_mom`; `dfl_engine_coef*` | the cases DK, DU, DT (zero or fixed start), SK, SU, ST, c <= 20 | 0.005 to 0.18 (quantiles 0.02; the density 0.03 against the general engine's) |
| `ltratio_lincdf`; `dfl_engine_lincdf` | case DK, a != 0, -100 <= c <= 4; the six cases through their designs, -100 <= c <= 2 | 0.07 to 0.24 |
| `power_envelope`; `dfl_engine_envelope` | case DK (DU), c <= -1e-4 | 0.13 to 0.16 |
| `ltratio_mom`; `dfl_engine_mom` | the six cases, any c0 and fixed start, c <= 20 | 0.013 to 0.23 |
| `pt_cdf`, `pt_pdf`, `pt_quantile`, `pt_power`, `pt_envelope`; `dfl_det_pt*` | any design, c <= 4 | 0.08 to 0.15 (quantiles 0.3 to 0.4, envelopes 0.2) |
| `dfl_det_coef_cdf_gp` (C only: the MATLAB programs have no `coef_cdf_gp`) | zero or stationary start, c <= 5 | 0.05 to 0.1 |

`opts.method` (`dfl_inv_opts.method` in C): `'auto'`, `'inversion'` (`'jets'` for the moments) to force them, `'engine'` for
the computation of before, which reproduces those values, `'grid'` for the fixed grid of the point-optimal statistic and
of `coef_cdf_gp`; `opts.tol` (1e-12, for the moments 1e-10 relative), `opts.verified` (recomputes by the engines the values
whose error estimate exceeds the tolerance), and every value comes with its error estimate (a second output in MATLAB and
the MEX file, `errest` in C; a warning names the values whose estimate exceeds the tolerance). At a = 0 (the law of U
alone) and a = NaN the linear form stays with the engines' closed forms. The regions, the tolerances and the warnings
are those of the Julia package, whose README lists the comparisons behind them and the cases where the engines' one
evaluation was wrong (the coefficient statistic of case DU with a fixed start far from the unit root, the linear form and
the moments with a fixed start above the unit root). The same devices fix the branches inside the engines: P(t), M_0(t)
and the labels of the roots of the classical line and contour integrals (`dfl_tparts`, `dfx_core('tparts')`) and the
argument of log P(t) in the provider of the general engine's line integral are determined at each point instead of being
continued along the grids, to 4e-13 the values of before, at the same cost. C and MATLAB agree with DFLimit.jl to 6e-15 in
the single inversions and the moments (`c/dfl_test_kan`, `matlab/kan_test`).

## Installation

Clone the repository (or download it as a zip) and add the folder `matlab` to the MATLAB path:
`addpath /path/to/dflimit-matlab/matlab` (and `savepath` to keep it).

The MEX file needs nothing else: no compiler, no runtime, one file of about 2.5 MB. There are two ways to get it.

1. Download a built one. The workflow "MATLAB MEX" (`.github/workflows/mex.yml`), which runs on a push that changes `c/`, on a tag `v*`
   and by hand from the Actions tab ("Run workflow"), builds one zip per platform (`ltratio-linux-x86_64.zip`, `ltratio-macos.zip`
   with the Apple silicon and the Intel MEX files, `ltratio-windows-x86_64.zip`) and attaches them to the run as
   artifacts; on a tag it also attaches them to a GitHub release of the repository
   (https://github.com/s-broda/dflimit-matlab/releases). Unzip, and add the folder to the
   path ahead of `matlab/`: `addpath ltratio-macos`.
2. Build it yourself, on macOS or Linux from the repository root with a C and a C++ compiler, or on Windows with the
   Visual Studio tools; `c/mex/README.md` has the commands, the timings and the details of the drop-ins.

## Tests

The C engine: `sh c/build.sh` compiles the engine and its test programs in `c/`; each test compares the engine with the
reference values of the Julia engine in `c/test/`:

```sh
cd c
for t in prim single gen fixed mom dbl app det kan ext; do ./dfl_test_$t test/${t}_ref.txt; done
./dfl_test_batch test/single_ref.txt    # the fixed node layouts against the adaptive rules
./dfl_test_vmath
./dfl_test_bounds
```

Every program ends with a line that counts the checks and the failures (`dfl_test_fixed` also reads
`test/contour_ref.txt`, the reference values of the contour integral; `dfl_test_bounds` has its reference values in
the source: the points of the repairs of 0.5.0, where the line integral exceeds the bounds of the engine or, with a
fixed start, climbed, and those of the upper tail above the unit root of November 23, 2026). The workflow "C engine" (`.github/workflows/c.yml`)
runs the same tests on Linux, macOS and Windows, on a push that changes `c/`, on a tag `v*` and by hand.

The MATLAB programs: `ltratio_test` runs 337 checks in about seven minutes and returns the number of failures
(`ltratio_test(true)` is the quick version, `ltratio_test(false, true)` raises an error on a failure); `dfx_test` checks
the line integral against stored values of the Julia engine (`matlab/validation/python_exact_reference.json`); `kan_test` checks the devices of Bao and Kan against the values of
DFLimit.jl in `matlab/validation/kan_ref.txt` (a copy of `c/test/kan_ref.txt`); `ext_test` checks the extensions of the
limit against `matlab/validation/ext_ref.txt` (a copy of `c/test/ext_ref.txt`) and, where it finds the MEX file, against
the C engine. Start
MATLAB in a folder other than `matlab/`. With the folder of the MEX file ahead of `matlab/` on the path the same
`ltratio_test` runs on the C engine, in about 35 seconds, which is the check of a built MEX file.

`matlab/validation/` has the script behind the accuracy tables of earlier versions of the Supplementary Material
(`acc_table.m`), the timings behind the automatic choice of the representation (`method_sweep.m` with its logs) and a
record of older computations (`tratio_legacy/`).

The reference values in `c/test/*_ref.txt` and the tables `c/tables.c` are written from the Julia engine by
`c/test/*_ref.jl` and `c/gen_tables.jl`, in the Julia environment `c/test/Project.toml`, which loads DFLimit.jl from
its GitHub repository: once `julia --project=c/test -e 'using Pkg; Pkg.instantiate()'`, then for example
`julia --project=c/test -t auto c/test/single_ref.jl c/test/single_ref.txt` from the repository root. The reference
values of the general engine (`c/test/det_ref.jl`, which writes `c/test/det_ref.txt`) and its tables
(`c/gen_det_tables.jl`, which writes `c/dfl_det_tables.c`) need a DFLimit.jl with the general engine (branch
engine-merge or later) in the environment, and those of the devices of Bao and Kan (`c/test/kan_ref.jl`, which writes
`c/test/kan_ref.txt`, copied to `matlab/validation/kan_ref.txt`) and `c/test/prim_ref.txt` (the branch at the point) a
DFLimit.jl with the module DFKan (commit c9fecd6 of its main branch or later), and the sections of the rules of one
evaluation in `c/test/det_ref.jl` and `matlab/validation/det_ref.jl` a DFLimit.jl with those rules (commit d95305c or
later). `c/test/dbl_ref.txt` and `c/test/contour_ref.txt` hold checked values (`verified = true`), against which the
tests compare the checked computation and, at looser tolerances, one evaluation.

## The C engine upstream of R and Python

The C sources in `c/` are the upstream of the R package [`dflimit`](https://github.com/s-broda/dflimit) and of the
Python package dflimit-py: each copies the engine sources of `c/` and `c/faddeeva/` (not the test programs, `c/test/`
or `c/mex/`) into its own tree with its `tools/sync-engine.sh /path/to/dflimit-matlab` and records the commit it copied from. A change to the
engine is made here, tested here, and then synced into the two packages. The Julia engine, of which the C engine is a
port, lives in the Julia package [`DFLimit.jl`](https://github.com/s-broda/DFLimit.jl), which also has the same
functions with the same names and arguments.

The scripts that reproduce the tables and figures of the companion paper are in its replication package, which runs
them on these MATLAB programs and on DFLimit.jl.

## Licence and citation

MIT (see `LICENSE`); the Faddeeva package in `c/faddeeva/` is by Steven G. Johnson, MIT licence (see the header of
`c/faddeeva/Faddeeva.cc`). Please cite the lead paper when using the programs, and the companion paper as well where
you use the representations of the classical cases (`CITATION.cff` has both).

The code was developed in the companion paper's repository up to commit b484d4f (see `NEWS.md`).
