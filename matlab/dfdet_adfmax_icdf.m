function q = dfdet_adfmax_icdf(p,D,c,opts)
%DFDET_ADFMAX_ICDF  Quantiles of the limiting ADFmax statistic for any deterministic terms.
%
%   q = DFDET_ADFMAX_ICDF(p,D,c,opts): the safeguarded Newton's method on
%   DFDET_ADFMAX_CDF and DFDET_ADFMAX_PDF from opts.x0 (-2) to opts.tol (1e-10).
%
%   Example: dfdet_adfmax_icdf(0.05, dfdet_trend(0))   % -2.43542
%
%   See also DFDET_ADFMAX_CDF, ADFMAX_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','dfdet_adfmax_icdf: p must lie in (0,1).'); end
q = dfdet_stat('adfmax','icdf',p,D,c,opts);
end
