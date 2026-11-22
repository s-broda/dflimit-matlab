function q = lm_quantile(p,D,c)
%LM_QUANTILE  Critical values of the LM unit root statistic for the trend functions of D.
%
%   q = LM_QUANTILE(p,D,c): the p-quantile of the limit of the LM statistic
%   (LM_CDF; the critical value of the level-p test, which rejects for small
%   values), from the quantile of V by bisection in log V (70 steps).
%
%   Example: lm_quantile([0.01 0.05 0.10], dfdet_trend(1))   % -3.5617 -3.0230 -2.7501
%
%   See also LM_CDF, LM_POWER.
if nargin < 3 || isempty(c), c = 0; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','lm_quantile: p must lie in (0, 1).'); end
CP = dfdet_lmpath(D,c);
q = zeros(size(p));
for k = 1:numel(p)
    lo = 1e-6; hi = 2;
    while dfdet_lmpath(CP,hi) < p(k), hi = 2*hi; end
    for it = 1:70
        m = sqrt(lo*hi);
        if dfdet_lmpath(CP,m) > p(k), hi = m; else, lo = m; end
    end
    q(k) = -0.5/sqrt(sqrt(lo*hi));
end
end
