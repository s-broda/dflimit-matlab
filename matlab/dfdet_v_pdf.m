function f = dfdet_v_pdf(y,D,c,opts)
%DFDET_V_PDF  Density of the denominator V of the t-ratio for any deterministic terms.
%
%   f = DFDET_V_PDF(y,D,c,opts): the inversion (1/2 pi i) int_C M_V(t) e^{-t y} dt;
%   arguments as in DFDET_V_CDF.
%
%   Example: dfdet_v_pdf(0.05, dfdet_trend(1))
%
%   See also DFDET_V_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
f = dfdet_stat('v','pdf',y,D,c,opts);
end
