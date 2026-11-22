function q = msb_quantile(p,D,c,opts)
%MSB_QUANTILE  Quantiles of the limit of MSB with GLS detrending (by bisection on its distribution function, 60 steps).
%
%   Example: msb_quantile(0.05, dfdet_perron('B', 0.5), 0, struct('gls', -17.8))   % 0.14427
%
%   See also MSB_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'gls') || isempty(opts.gls), opts.gls = -13.5; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','msb_quantile: p must lie in (0, 1).'); end
q = zeros(size(p));
for k = 1:numel(p)
    lo = 1e-4; hi = 1;
    while dfdet_v_cdf(hi,D,c,opts) < p(k), hi = 2*hi; end
    for it = 1:60
        m = (lo + hi)/2;
        if dfdet_v_cdf(m,D,c,opts) > p(k), hi = m; else, lo = m; end
    end
    q(k) = sqrt((lo + hi)/2);
end
end
