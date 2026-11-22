function y = dfdet_coef_icdf(p,D,c,opts)
%DFDET_COEF_ICDF  Quantiles of the limiting coefficient statistic for any deterministic terms.
%
%   y = DFDET_COEF_ICDF(p,D,c,opts): the p-quantiles (p an array in (0,1)) of
%   the limit of n(phi^ - phi0), by Newton's method on DFDET_COEF_CDF - p with
%   DFDET_COEF_PDF from the same evaluation, from opts.y0 (-10), each step
%   limited to 5, at most 40 steps, then the safeguarded Newton's method;
%   opts.tol (1e-9) the tolerance.  Other options as in DFDET_CDF.
%
%   Example: dfdet_coef_icdf(0.05, dfdet_trend(-1))   % -8.0391, the 5% critical value of rho^
%
%   See also DFDET_COEF_CDF, COEF_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','dfdet_coef_icdf: p must lie in (0,1).'); end
y = dfdet_stat('coef','icdf',p,D,c,opts);
end
