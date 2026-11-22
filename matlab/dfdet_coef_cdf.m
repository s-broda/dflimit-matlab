function F = dfdet_coef_cdf(y,D,c,opts)
%DFDET_COEF_CDF  Limiting distribution function of the coefficient statistic for any deterministic terms.
%
%   F = DFDET_COEF_CDF(y,D,c,opts) is the limit of P[n(phi^ - phi0) <= y] (y an
%   array) for the least squares estimate phi^ in the regression of DFDET_CDF,
%   phi0 = 1 + c0/n (opts.c0, 0), under phi = 1 + c/n.  The limit is
%   (U - c0 V)/V, and P[(U - c0 V)/V <= y] = P[U - (c0 + y) V <= 0] is the
%   distribution function of the t-ratio at 0 with the local null value c0 + y
%   (the contour integral; options as in DFDET_CDF).
%
%   Example: dfdet_coef_cdf(-20, dfdet_perron('C', 0.3))
%
%   See also DFDET_COEF_PDF, DFDET_COEF_ICDF, COEF_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
F = dfdet_stat('coef','cdf',y,D,c,opts);
end
