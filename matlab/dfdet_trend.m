function D = dfdet_trend(p)
%DFDET_TREND  The polynomial trend (1, r, ..., r^p) as a design of the general engine.
%
%   D = DFDET_TREND(p): p = -1 no deterministic terms (case DK), 0 an
%   intercept (DU), 1 a linear trend (DT), 2 a quadratic trend, and so on, in
%   the monomial basis.  For p >= 7 DFDET_LTREND (the same span in an
%   orthonormal basis) is the better conditioned form.
%
%   Example: the 5% critical value of the t-ratio with a cubic trend
%      dfdet_icdf(0.05, dfdet_trend(3))
%
%   See also DFDET_LTREND, DFDET_PERRON, DFDET_DETERM, DFDET_CDF.
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= -1), error('dfdet:args','dfdet_trend: p must be an integer >= -1.'); end
K = p + 1; deg = max(p,0);
C = zeros(deg + 1,K);
for i = 1:K, C(i,i) = 1; end
switch p
    case -1, name = 'no deterministic terms (DK)';
    case 0, name = 'intercept (DU)';
    case 1, name = 'linear trend (DT)';
    case 2, name = 'quadratic trend';
    case 3, name = 'cubic trend';
    otherwise, name = sprintf('polynomial trend of order %d',p);
end
D = dfdet_determ([0 1],{C},name);
end
