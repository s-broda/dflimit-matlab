function [F,info] = ltratio_cdf(x,cas,c,opts)
%LTRATIO_CDF  Limiting distribution function of the Dickey-Fuller t-ratio.
%
%   F = LTRATIO_CDF(x,cas,c,opts) is the limit as n -> Inf of P[t <= x] for
%   the least squares t-ratio t of the hypothesis phi = phi0 in the Gaussian
%   first-order autoregression y_t = alpha + delta t + phi y_{t-1} + sigma e_t
%   under the local alternatives phi = 1 + c/n, phi0 = 1 + c0/n (opts.c0,
%   default 0: the Dickey-Fuller tests of phi = 1), no drift (alpha = 0,
%   delta = 0).  cas is the case: the initial observation, Deterministic
%   (fixed) or Stationary, and the deterministic terms of the regression,
%   intercept Known to be zero, Unknown intercept, or intercept and Trend:
%      'DK'  y_0 fixed, regression of y_t on y_{t-1}         (tau)
%      'DU'  y_0 fixed, regression on 1, y_{t-1}               (tau_mu)
%      'DT'  y_0 fixed, regression on 1, t, y_{t-1}            (tau_tau)
%      'SK', 'SU', 'ST'  the same regressions with y_0 drawn from the
%            stationary distribution (c < 0)
%   With a fixed start y_0 = g0 sigma sqrt(n) (opts.g0, default 0), the start
%   enters the limit through the Ornstein-Uhlenbeck process started at g0.
%   The case 'DG' is the GLS-detrended trend case of Elliott, Rothenberg and
%   Stock (1996), detrended at cbar = -13.5 (the DF-GLS tests of the
%   paper): zero start, the line integral representation only, no
%   moments, and no tail recomputation (the conditional integral has no
%   such case, so far in the tails the values carry the line integral's
%   absolute error, see below).  c defaults to 0 (the unit root).  x may be
%   an array.
%
%   The limit is tau = (U - c0 V)/sqrt(V), with U = int_0^1 J~ dJ and
%   V = int_0^1 J~^2 dr, J the Ornstein-Uhlenbeck process dJ = cJ dr + dW
%   and J~ its residual on the deterministic terms of the regression (the
%   paper's section on the t-ratios and their limits).  Three exact
%   representations are implemented:
%      'single'  the line integral representation (or 'line'): the inversion in
%                the argument of U is explicit (a scaled chi-square(1) density
%                in case DK, a Bessel K_0 density otherwise) and one integral
%                along Re t = tau remains (engine DFX_CORE);
%      'double'  the conditional integral representation (or 'conditional'): the
%                inner inversion in the argument of V along a path through a
%                saddle point, the outer integral over V = b (engine TR_LIMDIST);
%      'contour' the contour integral, for a fixed start g0 ~= 0 in cases DU,
%                DT only (an error elsewhere): the inversion formula on paths
%                around the singularities of the moment generating function,
%                with the integral over v in closed form and the paths in s
%                and t deformed, a double
%                integral with an explicit integrand (engine TR_CONTOUR).
%   opts.method selects it: 'auto' (the default) takes, point by point, the
%   fastest that one evaluation gives to well below 1e-5 (LTRATIO_METHOD:
%   the line integral for c0 = 0, with a zero start and c0 > 0, in case DK
%   with c0 < 0 for c > -5, in case SK for c > -0.5 and far from the unit
%   root; with c0 ~= 0 otherwise the conditional integral, in cases SU, ST
%   only for c0 < 0 at x >= 0; the contour integral with a fixed start in
%   cases DU, DT, the conditional integral there for c > 1 with c0 < 0).
%   With a fixed nonzero start (g0 ~= 0) these programs have no line
%   integral (the C and Julia implementations carry its noncentral kernel):
%   'single' is taken as 'auto' in cases DU, DT, and in case DK the
%   conditional integral is used whatever opts.method says.
%
%   opts.verified (false by default): false is one evaluation by one method
%   at each point (no second evaluation, no refinement driven by an error
%   estimate, no fallback to another representation, no recomputation),
%   with the warning ltratio:brittle where the parameters lie in a regime in
%   which one evaluation is known to be brittle (LTRATIO_BRITTLE); true is
%   the checked behaviour at 1 to 3.5 times the cost of the default (about
%   1.7 with the contour integral; more in the tails below 1e-3, which it
%   recomputes): the conditional
%   integral with the agreement of two paths and the refinement of their
%   rules, the contour integral as two evaluations on different paths
%   refined until they agree, error estimates (info.err) and a status
%   (info.status), and with 'auto' the rule before November 2026, the
%   flagged values given the conditional integral's value and the values
%   below 1e-3 recomputed by the conditional integral where that is checked
%   here (c >= -20, or c >= -40 with (c g0)^2/2 <= 100).  The
%   representations agree to about 1e-12 (LTRATIO_TEST).
%
%   Range of c.  The rules of both representations are set for c from about
%   -480 to 10; at c = -60, -100, -200 and -480 the line integral is checked
%   to 1e-11 by node doubling, a second contour and the moments of the
%   density, and agrees with the C and Julia engines to 1e-11.  For c < -20
%   the panels and steps of the line integral follow c (DFX_CORE: WMAXC,
%   HUFAC, FILONPANELS, K0EDGE1), and the conditional integral of cases DU,
%   DT starts at E V e^-3.  With 'auto' the line integral serves alone with a
%   stationary start for c < -50 and with c0 ~= 0 for c < -100, where the
%   conditional integral flags most points or is off by 1e-9 without a flag
%   (LTRATIO_METHOD).  Below c = -550 an error is raised (the moments fail
%   from about c = -600, the line integral from about c = -650).  With
%   a fixed start g0 ~= 0 an error is raised outside the region of
%   LTRATIO_ARGS (case DK c <= 2, |g0| <= 2, |g0| <= 1 for c > 0; case DU
%   -100 <= c <= 2, |g0| <= 1.5, (c g0)^2/2 <= 4050; case DT -100 <= c <= 2,
%   |g0| <= 3, (c g0)^2/2 <= 7200, the region of the C and Julia
%   implementations, which the contour integral covers), and where the
%   conditional integral is the method (case DK, and 'double' in cases DU,
%   DT) also outside c >= -20 (c >= -40 where (c g0)^2/2 <= 100); inside, a
%   value the method cannot certify comes with its status (and a warning
%   with 'double' or 'contour'), and the MEX file has the line integral in
%   case DK.
%
%   Further options are passed to the engine: for 'single' those of DFX_CDF
%   (tau, hu, N, ...), for 'double' those of TR_LIMDIST (M, Nv, otol, ...).
%   info: method, time, and the engine's diagnostics (for 'double' and
%   'contour' the status and error estimates of TR_LIMDIST and TR_CONTOUR,
%   where a nonzero status raises a warning; with 'auto' a point whose status
%   is nonzero gets the line integral's value instead, info.nfallback counting
%   them, and with a fixed start in cases DU, DT the conditional integral's
%   where that is clean, see LTRATIO_DIST).
%
%   Examples
%      ltratio_cdf(-2.86, 'DU')                       % P[tau_mu <= -2.86] at the unit root
%      ltratio_cdf(-1.94087, 'DK', -5)                % power of the 5% test of tau at c = -5
%      ltratio_cdf(-3.41, 'ST', -10)                  % trend case, stationary start
%      ltratio_cdf(0, 'DK', 2, struct('c0',2))        % t-ratio centered at the true value
%      ltratio_cdf(-2, 'DU', -5, struct('g0',1.5))    % fixed start of 1.5 sigma sqrt(n) (the contour integral)
%
%   See also LTRATIO_PDF, LTRATIO_ICDF, LTRATIO_MOM, LTRATIO_TEST, DFX_CDF, TR_LIMDIST, TR_CONTOUR.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[F,info] = ltratio_dist(x,cas,c,opts,'cdf');
F = ltratio_rangecheck('ltratio_cdf',cdf_args(cas,c,opts),'x',x,F);   % (the final check of the values)
end

function a = cdf_args(cas,c,opts)
% the arguments of the call, for the message of LTRATIO_RANGECHECK
if ischar(cas) || isstring(cas), a = sprintf('case %s, c = %g',upper(char(cas)),min(c)); else, a = sprintf('a design, c = %g',min(c)); end
if isfield(opts,'c0') && ~isempty(opts.c0) && opts.c0 ~= 0, a = [a sprintf(', c0 = %g',opts.c0)]; end
if isfield(opts,'g0') && ~isempty(opts.g0) && opts.g0 ~= 0, a = [a sprintf(', g0 = %g',opts.g0)]; end
end
