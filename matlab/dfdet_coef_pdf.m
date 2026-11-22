function f = dfdet_coef_pdf(y,D,c,opts)
%DFDET_COEF_PDF  Limiting density of the coefficient statistic for any deterministic terms.
%
%   f = DFDET_COEF_PDF(y,D,c,opts): the density of the limit (U - c0 V)/V of
%   n(phi^ - phi0) at y, the derivative of DFDET_COEF_CDF in y, by the kernel
%   2 v^3 of the contour integral (the derivative of the distribution function
%   of the t-ratio at 0 in its local null value).  Arguments as in DFDET_COEF_CDF.
%
%   Example: dfdet_coef_pdf(-10, dfdet_trend(2))
%
%   See also DFDET_COEF_CDF, DFDET_COEF_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
f = dfdet_stat('coef','pdf',y,D,c,opts);
end
