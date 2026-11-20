function q = dfx_quantile(prob, model, c, opts)
%DFX_QUANTILE  Quantile of the limiting Dickey-Fuller t-ratio (port of DFExact.quantile):
%   Newton's method on DFX_CDF with DFX_PDF as derivative, safeguarded by bisection.
%   q = DFX_QUANTILE(prob, model, c, opts); opts as in DFX_CDF, plus q0 (start) and tol (1e-9).
if nargin < 3, c = 0; end
if nargin < 4, opts = struct(); end
q = zeros(size(prob));
for k = 1:numel(prob), q(k) = one(prob(k), model, c, opts); end
end

function q = one(prob, model, c, opts)
q0 = []; if isfield(opts, 'q0'), q0 = opts.q0; end
tol = 1e-9; if isfield(opts, 'tol'), tol = opts.tol; end
if isempty(q0), q0 = [-0.5, -1.6, -2.2]; q0 = q0(model); end
[q, ~, ~, conv] = ltratio_newton(@(y) dfx_cdf(y, model, c, opts), @(y) dfx_pdf(y, model, c, opts), q0, 1, tol, prob);
if ~conv, error('dfx_quantile: no convergence'); end
end
