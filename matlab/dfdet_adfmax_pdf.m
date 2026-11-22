function f = dfdet_adfmax_pdf(x,D,c,opts)
%DFDET_ADFMAX_PDF  Limiting density of the ADFmax statistic for any deterministic terms; arguments as in DFDET_ADFMAX_CDF.
%
%   Example: dfdet_adfmax_pdf(-2.4, dfdet_perron('C', 0.3))
%
%   See also DFDET_ADFMAX_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
f = dfdet_stat('adfmax','pdf',x,D,c,opts);
end
