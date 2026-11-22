function F = dfdet_cdf(x,D,c,opts)
%DFDET_CDF  Limiting distribution function of the t-ratio for any deterministic terms (the general engine).
%
%   F = DFDET_CDF(x,D,c,opts) is the limit as n -> Inf of P[t <= x] for the
%   least squares t-ratio t of phi = 1 + c0/n (opts.c0, default 0) in the
%   regression of y_t on d(t/n) and y_{t-1} (or, for a design built by
%   DFDET_DETREND_FIRST, the statistic that detrends first), y_t = beta' d(t/n)
%   + u_t, u_t = phi u_{t-1} + e_t, under phi = 1 + c/n (c defaults to 0).
%   D is the design of the deterministic terms (DFDET_TREND, DFDET_PERRON,
%   DFDET_FDESIGN, DFDET_LTREND, DFDET_DESIGN_FROM_MATRIX, DFDET_DETERM): any
%   piecewise smooth regressors, polynomial trends of any order, breaks in
%   level and slope at several known dates.  x may be an array.
%
%   The limit is (U - c0 V)/sqrt(V) with U = int Jt dJ and V = int Jt^2, Jt the
%   residual of the Ornstein-Uhlenbeck process J on the regressors.  The joint
%   moment generating function of (U, V) is e^{-s/2} D(s,t)^{-1/2} with D a small
%   determinant (DFDET_MGF), a polynomial in s whose roots are followed along
%   an outer contour in t (DFDET_CPATH); the distribution function is the
%   inversion on paths around the roots, a double integral of functions in
%   closed form (DFDET_CONTOUR).
%
%   Options: c0 (the local null value, 0); start ('zero', the default, 'fixed'
%   with J(0) = g0, y_0 = g0 sigma sqrt(n), or 'stationary', J(0) ~ N(0,
%   -1/(2c)), c < 0); g0 (0; g0 ~= 0 implies the fixed start); gls (GLS
%   detrending at abar = 1 + gls/n, gls = cbar < 0, the data quasi-differenced
%   as by Elliott, Rothenberg and Stock; zero start only); verified (false: one
%   evaluation, with the warnings dfdet:rounding, dfdet:contour and dfdet:line
%   where the computation is known to be brittle; true: the checked
%   computation of DFDET_VERIFY, about twice the cost); vtol (1e-9); and the
%   rules phi, hfac, n, ntail, alpha, filon, tol (DFDET_CPATH, DFDET_CONTOUR).
%   method: 'auto' (the default) or 'paths', the inversion on paths; 'line',
%   the line integral with the explicit kernel, is in the MEX file only
%   (DFDET_METHOD).
%   Non-finite x: F = 0 at -Inf, 1 at +Inf, NaN at NaN.
%
%   Examples
%      dfdet_cdf(-3.5, dfdet_perron('C', 0.3))                  % break in level and slope at 0.3
%      dfdet_cdf(-3.5, dfdet_perron('C', 0.3), -10)             % its power at c = -10
%      dfdet_cdf(-2.9, dfdet_trend(1), -5, struct('gls', -13.5)) % the DF-GLS test
%      dfdet_cdf(-2, dfdet_trend(2), -5, struct('start', 'stationary'))
%
%   See also DFDET_PDF, DFDET_ICDF, DFDET_MOM, LTRATIO_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
F = dfdet_dist(x,D,c,opts,'cdf');
end
