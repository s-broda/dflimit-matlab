function b = dfdet_lin_icdf(p,a,D,c,opts)
%DFDET_LIN_ICDF  Quantiles of the linear form U - aV for any deterministic terms.
%
%   b = DFDET_LIN_ICDF(p,a,D,c,opts): b with P(U - aV <= b) = p (p an array in
%   (0,1)), by the safeguarded Newton's method on DFDET_LIN_CDF and DFDET_LIN_PDF
%   from opts.b0 (0) to opts.tol (1e-10).  Arguments as in DFDET_LIN_CDF.
%
%   Example: dfdet_lin_icdf(0.05, -3.5, dfdet_trend(-1))   % -0.2639
%
%   See also DFDET_LIN_CDF, DFDET_POWER_ENVELOPE.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
if ~(isnumeric(a) && isscalar(a) && isreal(a) && isfinite(a)), error('dfdet:args','dfdet_lin_icdf: a must be a finite real scalar.'); end
b = dfdet_stat('lin','icdf',p,D,c,opts,a);
end
