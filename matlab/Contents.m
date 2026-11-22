% The exact limiting distribution of the Dickey-Fuller t-ratio
% Version 0.5.0 21-Nov-2026
%
% MATLAB programs accompanying
%
%   Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of
%   the Dickey-Fuller t-Ratio and Some Related Statistics," Purdue University,
%   Lucerne University of Applied Sciences and University of Toronto.
%
% Model: y_t = alpha + delta*t + phi*y_{t-1} + sigma*e_t, t = 1..n, e_t i.i.d.
% N(0,1), no drift, phi = 1 + c/n, and t the least squares t-ratio of the
% hypothesis phi = 1 + c0/n (c0 = 0: the Dickey-Fuller tests of phi = 1) from
% the regression of y_t on y_{t-1} and the deterministic terms.  As n -> Inf,
% t => (U - c0 V)/sqrt(V), U = int J~ dJ, V = int J~^2, J the Ornstein-
% Uhlenbeck process dJ = cJ dr + dW started at g0 = y_0/(sigma sqrt(n)), J~
% its residual on the deterministic terms.  Seven cases, selected by a
% two-letter code (D: deterministic, i.e. fixed, start; S: stationary start):
%   'DK'  regression of y_t on y_{t-1}          (tau),      start g0 (opts.g0)
%   'DU'  regression on 1 and y_{t-1}            (tau_mu),   start g0
%   'DT'  regression on 1, t and y_{t-1}         (tau_tau),  start g0
%   'SK', 'SU', 'ST'  the same regressions with y_0 from the stationary
%         distribution, J(0) ~ N(0, -1/(2c)), c <= -1e-4
%   'DG'  the GLS-detrended trend case of Elliott, Rothenberg and Stock
%         (cbar = -13.5), zero start: the DF-GLS tests
% The default start is g0 = 0, the standard local-to-unity setting.
%
% Limiting distribution, quantiles and moments of t
%   ltratio_cdf   - distribution function P[t <= x]
%   ltratio_pdf   - density of t
%   ltratio_icdf  - quantiles (asymptotic critical values) of t
%   ltratio_mom   - moments E[t^j], j = 1..p
%   Options: c0 (the local null value, 0), g0 (fixed start, 0), method
%   ('auto', the default: the fastest exact representation that one
%   evaluation gives to well below 1e-5, point by point, see ltratio_method;
%   'single' (or 'line') and 'double' (or 'conditional') force the line
%   integral and the conditional integral representation; 'contour', the
%   contour integral, for a fixed start g0 ~= 0 in cases DU, DT only, where
%   'auto' takes it too; a fixed start in case DK uses 'double'), verified
%   (false: one evaluation by one method, with the warning ltratio:brittle in
%   regimes known to be brittle; true: error estimates, a status, fallbacks
%   and the recomputation of the tails).
%   The cases are DK, DU, DT (zero or fixed start), SK, SU, ST (stationary
%   start) and DG (the GLS-detrended trend case of Elliott, Rothenberg and
%   Stock, cbar = -13.5: zero start, line integral only, no moments).
%   ltratio_method - the automatic choice of the representation
%   ltratio_brittle - the warning of the one-evaluation default in a brittle regime
%
% The line integral representation (the inversion in the argument of U is
% explicit and one integral along the line Re t = tau remains)
%   dfx_core      - the engine: the transforms Phi_x(t), their partition into
%                   a smooth and an oscillating part, the trapezoidal and
%                   Filon rules on the line, the saddle factor for c > 0, the
%                   kernel along the parabola for c0 ~= 0, the moments
%   dfx_cdf, dfx_pdf, dfx_quantile, dfx_moments, dfx_tstar
%                 - front ends of dfx_core by model number (1, 2, 3 = K, U, T)
%                   and opts.start ('zero' or 'stationary'); a MATLAB port of
%                   the Julia implementation DFExact.jl
%   dfx_series    - Taylor coefficients of the determinants for the port
%   dfx_test      - the port against the stored reference values of DFExact.jl
%
% The conditional integral representation (the inner inversion in the
% argument of V along a path through a saddle point, the outer integral
% over V)
%   tr_limdist    - distribution function and density: the probability of a
%                   region of the plane in case DK (zero or fixed start; the
%                   inner probability is an error function), the form for
%                   the t-ratio with a zero or a stationary start in the
%                   other cases (closed-form Bessel integrands), the
%                   probability of a region with the noncentral probability
%                   for a fixed start in cases DU, DT
%   tr_limmom     - moments of the limit (any c, c0, start) as one-dimensional
%                   integrals of the moment generating function
%   tr_cflim      - the joint moment generating function of (U, V) with a
%                   zero start (the Nabeya-Tanaka determinants),
%                   its roots in the argument of U and its Girsanov forms
%   tr_nsstart    - the initial-value term of a fixed start in closed form
%                   (Nabeya and Sorensen's determinants in
%                   cases DU, DT; the Riccati solution in case DK)
%   tr_jbstart    - the tilted law of J(1) given the start, case DK
%   nt_series     - Taylor coefficients of the Nabeya-Tanaka determinants
%                   (generated by gen_series.py, Python with sympy)
%   faddeeva      - the Faddeeva function (complex error function)
%
% The contour integral representation (the inversion formula on paths
% around the singularities of the moment generating function, with the
% integral over v in closed form, the line in the argument of U replaced
% by two U-shaped paths around the roots of D(., t) and the line in t by a hyperbola: a double
% integral with an explicit integrand; a fixed start in cases DU, DT)
%   tr_contour    - distribution function and density (one evaluation; with
%                   opts.verified an error estimate from two evaluations on
%                   different paths with different rules, refined together; a port of tr_contour
%                   of julia/DFExact.jl; the roots and P(t) from dfx_core)
%
% The applications (local power, confidence sets and the weighted-average-power
% test in the paper; the same functions with the same arguments as in the
% Julia engine julia/DFExact.jl, with a zero or a stationary start as stated)
%   adfmax_cdf, adfmax_pdf, adfmax_icdf
%                 - the ADFmax test of Leybourne (1995), the larger of the
%                   forward and the reversed t-ratio: distribution function,
%                   density and quantiles of its limit (|R| - 1/2)/sqrt(V),
%                   cases DU, DT (zero start) and SU, ST (stationary start)
%   adfmax_dist   - argument checks of adfmax_cdf and adfmax_pdf
%   tr_adfmax     - the engine: the conditional integral over the region
%                   |R| <= 1/2 + x sqrt(V)
%   coef_cdf, coef_icdf
%                 - the coefficient statistic n(phi^ - phi_0), whose limit
%                   is (U - c0 V)/V: distribution function and quantiles
%   ltratio_region - the probability of a region {U <= g(V)} of the plane,
%                   case DK, zero start: a line, the parabola of the t-ratio,
%                   the boundary of the weighted-average-power test of
%                   the paper, or any function g
%   tr_region     - its engine (the conditional integral with the boundary
%                   a(V) = 1/2 + g(V))
%   ltratio_lincdf - the distribution of the linear form U - aV, or of V,
%                   case DK or the GLS trend case DG (inversion formula;
%                   for small |a| and c > 0 a conditional integral)
%   tr_gilpelaez  - its engine (Gil-Pelaez)
%   tr_linexact   - its route where the inversion is not trusted: the region
%                   {U <= aV + b} (case DK), tr_lincond (case DG)
%   tr_lincond    - the conditional integral given R = U + 1/2
%   power_envelope - the asymptotic power envelope of Elliott, Rothenberg and
%                   Stock (1996), case DK (= DU) or DG
%   ltratio_bisect - bisection (as ltratio_bisect of julia/DFExact.jl)
%   gls_parts, lam_gls - the GLS-detrended trend case at a given cbar
%
% The general engine (any deterministic terms: polynomial trends of any order, breaks in level and slope at
% several known dates, smooth functions of t/n; the joint moment generating function of U and V through a small
% determinant, its roots in the argument of U followed along an outer contour in t, and the inversion on paths
% around the roots, a double integral of functions in closed form; a port of the Julia engine DFLimit.DFDet)
%   dfdet_trend, dfdet_perron, dfdet_fdesign, dfdet_ltrend, dfdet_design_from_matrix, dfdet_detrend_first, dfdet_determ
%                 - the designs: polynomial trends, Perron's models A, B, C with one or several breaks, regressors
%                   given by a function (Chebyshev interpolants per piece), the orthonormal Legendre trend, the
%                   design of a regressor matrix (its breaks detected), the statistic that detrends first, and
%                   any piecewise polynomial regressors
%   dfdet_cdf, dfdet_pdf, dfdet_icdf, dfdet_mom
%                 - distribution function, density, quantiles and moments of the limiting t-ratio for a design
%                   (options c0, start 'zero', 'fixed' with g0, 'stationary', gls = cbar for GLS detrending,
%                   verified)
%   dfdet_coef_cdf, dfdet_coef_pdf, dfdet_coef_icdf - the coefficient statistic
%   dfdet_lin_cdf, dfdet_lin_pdf, dfdet_lin_icdf - the linear form U - aV
%   dfdet_v_cdf, dfdet_v_pdf - the denominator V
%   dfdet_adfmax_cdf, dfdet_adfmax_pdf, dfdet_adfmax_icdf - ADFmax
%   dfdet_power_envelope, dfdet_local_power, dfdet_region - the envelope (no deterministic terms, intercept,
%                   linear trend), local power of the left-tailed tests, regions bounded by lines and parabolas
%   pt_cdf, pt_pdf, pt_quantile, pt_power, pt_envelope
%                 - the point-optimal statistic of Elliott, Rothenberg and Stock and the Gaussian power envelope
%                   for any design (one inversion of Gil-Pelaez)
%   lm_cdf, lm_quantile, lm_power - the LM statistic (Schmidt and Phillips; with breaks, Lee and Strazicich)
%   mza_cdf, mza_quantile, msb_cdf, msb_quantile - MZ_alpha and MSB of Ng and Perron with GLS detrending
%   dfdet_prep    - the prepared design: the functionals ell, the matrices A, B of U + 1/2 and V - int J^2, the
%                   roots and their signs (least squares, GLS, the point-optimal and the LM statistics)
%   dfdet_model, dfdet_startmodel - the design with its start (the shift of the fixed start, the extended
%                   matrices of the stationary start), cached
%   dfdet_sigma   - the covariance Sigma(theta) of ell under the tilted Ornstein-Uhlenbeck law
%   dfdet_mgf, dfdet_dval, dfdet_logmgf, dfdet_fixe1 - the joint moment generating function: log D_V(t), the
%                   roots, the data of the fixed start, the determinant, its value at points (s, t)
%   dfdet_tstar   - t_*, the first zero of D_V
%   dfdet_cpath   - the outer contour (hyperbola or vertical line) with the roots and log P continued along it
%   dfdet_contour - the U-shapes around the root groups, the closed-form v-kernels, the outer sum, Filon's
%                   rule on the line
%   dfdet_verify  - the checked computation (opts.verified)
%   dfdet_dist, dfdet_stat, dfdet_args, dfdet_xspec, dfdet_qgrid, dfdet_lmpath, dfdet_newton, dfdet_cache,
%   dfdet_key, dfdet_peval, dfdet_msg - dispatch, argument checks, parabola terms, the grid of the point-
%                   optimal statistic, the safeguarded Newton's method, the caches, helpers
%   dfdet_make_paths - the path tables dfdet_paths.mat of cases DK, DU, DT at the unit root (from the tables of
%                   the Julia package)
%
% The two engines.  The user-level routines ltratio_cdf, ltratio_pdf, ltratio_icdf, ltratio_mom, coef_cdf,
% coef_icdf, adfmax_cdf, adfmax_pdf, adfmax_icdf, ltratio_lincdf, power_envelope and ltratio_region take
% opts.engine: 'classical' (the programs above), 'general' (the general engine) or 'auto' (the default: the
% engine of ltratio_engine, by the timings of the two MATLAB engines), and a design of the general engine in
% place of the case (a classical design is its case; any other goes to the general engine).
%   ltratio_engine - the rule of opts.engine = 'auto'
%   ltratio_spec, ltratio_route - the design, start and case of a computation; helpers of the router
%
% Validation
%   dfdet_test    - the general engine and the router against the Julia engine DFLimit.DFDet
%                   (validation/det_ref.txt, written by validation/det_ref.jl), with a failure count
%   ltratio_test  - validation suite with a failure count: the two
%                   representations against each other over the six cases,
%                   both sides of the unit root and the local null value;
%                   exact values (P[tau <= 0] in closed form, P[tau_mu <= 0]
%                   from the coefficient statistic, a one-dimensional
%                   inversion at c0 ~= 0, the asymptotic critical values of
%                   MacKinnon); the initial observation (continuity at a
%                   zero start, the stationary start as the Gauss-Hermite
%                   average over a fixed start, the closed-form noncentral
%                   probability in case DK, the rank-two structure of the
%                   initial-value term); moments (Nabeya's values, the two
%                   moment routines, the integral of the density); the
%                   density against differences of the distribution function;
%                   argument checks; the applications against the Julia
%                   engine; the contour integral against the Julia engine
%                   (a fixed start in cases DU, DT); the default against
%                   verified and the warnings.  About five minutes (30 s on
%                   the MEX drop-ins).
%   validation/   - the script behind the accuracy tables of earlier versions of the Supplementary Material,
%                   the timings behind the automatic choice of the representation and the reference
%                   values of dfx_test; the reference values of dfdet_test (det_ref.jl, det_ref.txt) and the
%                   timings behind the choice of engine (engine_table_matlab.m, engine_table_matlab.txt)
%
% Reproduction of the paper: not part of this package.  The paper's replication
% package has the scripts that recompute its tables, the numbers it quotes and
% the data of its figures on these programs (ltratio_paper.m, with
% mackinnon94_p.m) and on the Julia package DFLimit.jl.
%
% Example
%   ltratio_icdf([0.01 0.05 0.10],'DU')              % asymptotic critical values of tau_mu
%   ltratio_cdf(-2.86157,'DU',-5)                    % limiting power of the 5% test at c = -5
%   ltratio_cdf(-2.86157,'SU',-5)                    % the same with a stationary start
%   ltratio_cdf(-2.86157,'DU',-5,struct('g0',1.5))   % fixed start y_0 = 1.5 sigma sqrt(n) (contour integral)
%   ltratio_cdf(0,'DK',2,struct('c0',2))             % t-ratio centered at the true value
%   ltratio_mom(4,'DT')                              % first four moments of tau_tau
%   adfmax_icdf(0.05,'DU')                           % 5% critical value of ADFmax with a constant
%   coef_cdf(-14.0936,'DU',-10)                      % power of the 5% coefficient test at c = -10
%   dfdet_icdf(0.05,dfdet_perron('C',0.3))           % 5% critical value with a break in level and slope at 0.3
%   dfdet_cdf(-3.5,dfdet_perron('C',0.3),-10)        % its power at c = -10
%   ltratio_cdf(-2,'DU',-5,struct('engine','general')) % case DU by the general engine
%
% Requirements.  Base MATLAB (R2019b or later); runs on GNU Octave (tested
% with Octave 8), more slowly.
%
% License: MIT (see LICENSE).  Please cite the paper when using the programs.
