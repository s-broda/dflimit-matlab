function P = dfdet_lin_cdf(b,a,D,c,opts)
%DFDET_LIN_CDF  Distribution function of the linear form U - aV for any deterministic terms.
%
%   P = DFDET_LIN_CDF(b,a,D,c,opts) is the limit of P(U - aV <= b) (b an array)
%   under phi = 1 + c/n for the design D with its start and detrending (options
%   start, g0, gls as in DFDET_CDF); the point-optimal test of Elliott,
%   Rothenberg and Stock against cbar rejects for small U - aV.  With R = U + a0
%   (a0 the shift of Mt) the event is R <= X(sqrt V), X(v) = a v^2 + a0 + b, the
%   parabola of the t-ratio with c0 = a, x = 0 and the constant a0 + b
%   (DFDET_XSPEC; on the vertical line the part at the root of X by Filon's
%   rule).  a = NaN gives P(V <= b) (DFDET_V_CDF).
%
%   Example: dfdet_lin_cdf(-0.2639, -3.5, dfdet_trend(-1))   % 0.05 (the point-optimal test at cbar = -7)
%
%   See also DFDET_LIN_PDF, DFDET_LIN_ICDF, LTRATIO_LINCDF, DFDET_POWER_ENVELOPE.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
if ~(isnumeric(a) && isscalar(a) && isreal(a)), error('dfdet:args','dfdet_lin_cdf: a must be a real scalar (NaN: the distribution of V).'); end
P = dfdet_stat('lin','cdf',b,D,c,opts,a);
end
