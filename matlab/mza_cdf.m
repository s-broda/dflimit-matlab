function F = mza_cdf(y,D,c,opts)
%MZA_CDF  Limiting distribution function of MZ_alpha of Ng and Perron with GLS detrending, any design.
%
%   F = MZA_CDF(y,D,c,opts) is the limit of P(MZ_alpha <= y) for the modified
%   coefficient statistic of Ng and Perron with the data GLS-detrended at
%   cbar = opts.gls (-13.5): its limit is that of the coefficient statistic of
%   the GLS-detrended regression, (U - c0 V)/V (DFDET_COEF_CDF).  Other options
%   as in DFDET_COEF_CDF.
%
%   Example: mza_cdf(-23.659, dfdet_perron('B', 0.5), 0, struct('gls', -17.8))   % 0.05
%
%   See also MZA_QUANTILE, MSB_CDF, DFDET_COEF_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'gls') || isempty(opts.gls), opts.gls = -13.5; end
F = dfdet_coef_cdf(y,D,c,opts);
end
