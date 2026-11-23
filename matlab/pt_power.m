function [b,P] = pt_power(c,D,cbar,opts)
%PT_POWER  The critical value and the limiting power of the point-optimal test against cbar for the design D.
%
%   [b,P] = PT_POWER(c,D,cbar,opts): b the critical value of the level-alpha
%   test (opts.alpha, 0.05; PT_QUANTILE) and P its limiting power at each c in
%   the array c (PT_CDF at b), by default both by the single inversion of
%   Bao and Kan (opts.method, opts.tol as in PT_CDF; above c = 4 the fixed
%   grid for both).
%
%   Example: [b,P] = pt_power(-10, dfdet_perron('C', 0.5), -10)
%
%   See also PT_QUANTILE, PT_ENVELOPE.
if nargin < 4 || isempty(opts), opts = struct; end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
if ~(isnumeric(alpha) && isscalar(alpha) && alpha > 0 && alpha < 1), error('dfdet:args','pt_power: the level opts.alpha must lie in (0, 1).'); end
% (the computation of the largest c for the critical value and the power alike, as pt_method of DFLimit.jl: above c = 4
% with opts.method 'auto' the fixed grid for both)
if dfdet_ptmethod(opts,max(c(:)),'pt_power'), opts.method = 'grid'; end
o = opts; if isfield(o,'q0'), o = rmfield(o,'q0'); end
if isfield(o,'alpha'), o = rmfield(o,'alpha'); end
b = pt_quantile(alpha,D,cbar,opts);
P = zeros(size(c));
for k = 1:numel(c), P(k) = pt_cdf(b,D,cbar,c(k),o); end
end
