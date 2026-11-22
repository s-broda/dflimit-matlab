function f = dfdet_lin_pdf(b,a,D,c,opts)
%DFDET_LIN_PDF  Density of the linear form U - aV for any deterministic terms (a = NaN: of V).
%
%   f = DFDET_LIN_PDF(b,a,D,c,opts); arguments as in DFDET_LIN_CDF.
%
%   Example: dfdet_lin_pdf(0.3, -2, dfdet_trend(1), -5)
%
%   See also DFDET_LIN_CDF, DFDET_V_PDF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
if ~(isnumeric(a) && isscalar(a) && isreal(a) && ~isinf(a)), error('dfdet:args','dfdet_lin_pdf: a must be a finite real scalar (NaN: the distribution of V).'); end
f = dfdet_stat('lin','pdf',b,D,c,opts,a);
end
