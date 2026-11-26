function [q,P] = lm_power(c,D,opts)
%LM_POWER  The critical value and the limiting power of the LM test for the trend functions of D.
%
%   [q,P] = LM_POWER(c,D,opts): q the critical value of the level-alpha test
%   (opts.alpha, 0.05; LM_QUANTILE) and P its limiting power at each c in the
%   array c.  An error for a design whose span does not contain the constant
%   (LM_CDF).
%
%   Example: [q,P] = lm_power([-5 -10 -15 -20], dfdet_trend(1))   % P = 0.1059 0.2656 0.5149 0.7599
%
%   See also LM_CDF, LM_QUANTILE.
if nargin < 3 || isempty(opts), opts = struct; end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
if ~(isnumeric(alpha) && isscalar(alpha) && alpha > 0 && alpha < 1), error('dfdet:args','lm_power: the level opts.alpha must lie in (0, 1).'); end
q = lm_quantile(alpha,D);
P = zeros(size(c));
for k = 1:numel(c), P(k) = dfdet_lmpath(dfdet_lmpath(D,c(k)),1/(4*q^2)); end
end
