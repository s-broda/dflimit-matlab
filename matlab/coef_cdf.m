function F = coef_cdf(y,cas,c,opts)
%COEF_CDF  Limiting distribution function of the Dickey-Fuller coefficient statistic.
%
%   F = COEF_CDF(y,cas,c,opts) is the limit as n -> Inf of P[n(phi^ - phi0) <= y] (y an array) for the least squares
%   estimate phi^ of phi in the regression of LTRATIO_CDF, under phi = 1 + c/n, with phi0 = 1 + c0/n (opts.c0,
%   default 0: the coefficient test of Dickey and Fuller, rho^ = n(phi^ - 1)).  The limit is (U - c0 V)/V, so that
%   P[(U - c0 V)/V <= y] = P[U - (c0 + y)V <= 0], the distribution function of the t-ratio at 0 with the local null
%   value c0 + y: LTRATIO_CDF(0,cas,c,opts) with c0 + y in place of c0, in every case and with every start and
%   option of LTRATIO_CDF (method, g0).
%
%   Example: coef_cdf(-14.094,'DU',0)    % 0.05: the 5% critical value of rho^_mu
%
%   See also COEF_ICDF, LTRATIO_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(y) && isreal(y)), error('ltratio:args','coef_cdf: y must be real.'); end
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
F = zeros(size(y));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (the regimes c0 + y passes, each warned once)
for k = 1:numel(y)
    if isnan(y(k)), F(k) = NaN; elseif y(k) == Inf, F(k) = 1; elseif y(k) == -Inf, F(k) = 0;
    else, o = opts; o.c0 = c0 + y(k); F(k) = ltratio_cdf(0,cas,c,o); end
end
end
