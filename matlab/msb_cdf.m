function F = msb_cdf(x,D,c,opts)
%MSB_CDF  Limiting distribution function of MSB of Ng and Perron with GLS detrending, any design.
%
%   F = MSB_CDF(x,D,c,opts) is the limit of P(MSB <= x) for Stock's modified
%   Sargan-Bhargava statistic with the data GLS-detrended at cbar = opts.gls
%   (-13.5): MSB => sqrt(V_g), so that F(x) = P(V_g <= x^2) (DFDET_V_CDF) for
%   x > 0 and 0 for x <= 0.
%
%   Example: msb_cdf(0.144, dfdet_perron('B', 0.5), 0, struct('gls', -17.8))   % 0.0491
%
%   See also MSB_QUANTILE, MZA_CDF, DFDET_V_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'gls') || isempty(opts.gls), opts.gls = -13.5; end
F = zeros(size(x));
for k = 1:numel(x)
    if isnan(x(k)), F(k) = NaN; elseif x(k) <= 0, F(k) = 0; else, F(k) = dfdet_v_cdf(x(k)^2,D,c,opts); end
end
end
