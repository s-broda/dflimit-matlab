function F = dfdet_adfmax_cdf(x,D,c,opts)
%DFDET_ADFMAX_CDF  Limiting distribution function of the ADFmax statistic of Leybourne for any deterministic terms.
%
%   F = DFDET_ADFMAX_CDF(x,D,c,opts) is the limit of P(ADFmax <= x) for the larger
%   of the Dickey-Fuller t-ratios of the forward regression (y_t on d(t/n) and
%   y_{t-1}) and of the reverse regression (y_t on d(t/n) and y_{t+1}).  In the
%   limit both regressions have the denominator V and the numerators U and
%   -U - 1, so that {ADFmax <= x} = {|U + 1/2| <= 1/2 + x sqrt(V)}: two
%   parabolas in sqrt(V), each a term of the contour integral (DFDET_XSPEC).
%   The limit of the numerators is the derivation for the constant and the
%   linear trend; for other designs and starts this is the same functional of
%   the design's (U, V).  Options start, g0 as in DFDET_CDF (no GLS detrending).
%
%   Example: dfdet_adfmax_cdf(-2.43542, dfdet_trend(0))   % 0.05
%
%   See also DFDET_ADFMAX_PDF, DFDET_ADFMAX_ICDF, ADFMAX_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
F = dfdet_stat('adfmax','cdf',x,D,c,opts);
end
