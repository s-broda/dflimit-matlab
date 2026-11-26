function c0 = predictive_c0(c,rho,b)
%PREDICTIVE_C0  The local null value of the t-ratio of the slope in a predictive regression with a regressor local to unity.
%
%   c0 = PREDICTIVE_C0(c,rho,b) = c - b/rho for the regression y_t = mu +
%   beta x_{t-1} + u_t, x_t = phi x_{t-1} + v_t, phi = 1 + c/n, with rho the
%   correlation of the innovations (u_t, v_t) (in [-1, 1] and not 0) and the
%   local alternative beta = b sigma_u/(sigma_v n) (b = 0, the default: the
%   null hypothesis beta = 0).  With opts.rho = rho and opts.c0 =
%   PREDICTIVE_C0(c,rho,b), DFDET_CDF and LTRATIO_CDF of the design with an
%   intercept (dfdet_trend(0), case 'DU') give the limit of the statistic,
%   rho tau_c0 + sqrt(1 - rho^2) Z.
%
%   Example: ltratio_cdf(1.645, 'DU', 0, struct('rho', -0.95, 'c0', predictive_c0(0, -0.95)))   % 0.5787
%
%   See also CADF_C0, DFDET_CDF, LTRATIO_CDF.
if nargin < 3 || isempty(b), b = 0; end
if ~(isnumeric(rho) && isscalar(rho) && isreal(rho) && isfinite(rho) && rho >= -1 && rho <= 1 && rho ~= 0), error('dfdet:args','predictive_c0: rho must lie in [-1, 1] and not be 0'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c) && isnumeric(b) && isscalar(b) && isreal(b) && isfinite(b)), error('dfdet:args','predictive_c0: c and b must be finite'); end
c0 = double(c) - double(b)/double(rho);
end
