function P = pt_envelope(c,D,opts)
%PT_ENVELOPE  The Gaussian power envelope for any design: the point-optimal test against cbar = c at c.
%
%   P = PT_ENVELOPE(c,D,opts): at each c < 0 (an array) the power at c of the
%   point-optimal level-alpha test (opts.alpha, 0.05) against cbar = c (PT_POWER);
%   alpha for c >= 0.  For no deterministic terms, an intercept and a linear
%   trend it is DFDET_POWER_ENVELOPE (to 1e-8).
%
%   Example: pt_envelope(-10, dfdet_perron('C', 0.5))   % 0.1789
%
%   See also PT_POWER, DFDET_POWER_ENVELOPE, POWER_ENVELOPE.
if nargin < 3 || isempty(opts), opts = struct; end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
if ~(isnumeric(alpha) && isscalar(alpha) && alpha > 0 && alpha < 1), error('dfdet:args','pt_envelope: the level opts.alpha must lie in (0, 1).'); end
P = zeros(size(c));
for k = 1:numel(c)
    if c(k) < 0, [~,P(k)] = pt_power(c(k),D,c(k),opts); elseif isnan(c(k)), P(k) = NaN; else, P(k) = alpha; end   % (NaN at NaN, as POWER_ENVELOPE)
end
end
