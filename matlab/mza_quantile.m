function y = mza_quantile(p,D,c,opts)
%MZA_QUANTILE  Quantiles of the limit of MZ_alpha with GLS detrending (arguments as in MZA_CDF; opts.y0 = -20 the start).
%
%   Example: mza_quantile(0.05, dfdet_perron('B', 0.5), 0, struct('gls', -17.8, 'y0', -25))   % -23.659
%
%   See also MZA_CDF, DFDET_COEF_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'gls') || isempty(opts.gls), opts.gls = -13.5; end
if ~isfield(opts,'y0') || isempty(opts.y0), opts.y0 = -20; end
y = dfdet_coef_icdf(p,D,c,opts);
end
