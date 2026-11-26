function q = hegy_t3_quantile(p,D,c,opts)
%HEGY_T3_QUANTILE  Quantiles of the limit of the first seasonal t-ratio of Hylleberg, Engle, Granger and Yoo.
%
%   q = HEGY_T3_QUANTILE(p,D,c,opts): the p-quantiles (p in (0, 1), an array)
%   of the limit of t3 (HEGY_T3_CDF), by Newton's method with the density from
%   opts.x0 (-2), the step limited to 1, then the safeguarded search, to
%   opts.tol (1e-10).  opts.c0, opts.vol as in HEGY_T3_CDF.
%
%   Example: hegy_t3_quantile(0.05, dfdet_trend(0))      % -3.3442, the 5% critical value with seasonal dummies
%
%   See also HEGY_T3_CDF, HEGY_T3_PDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
fn = 'hegy_t3_quantile';
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','%s: p must lie in (0, 1)',fn); end
x0 = -2; if isfield(opts,'x0') && ~isempty(opts.x0), x0 = opts.x0; end
tol = 1e-10; if isfield(opts,'tol') && ~isempty(opts.tol), tol = opts.tol; end
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','%s: c0 must be a finite real scalar',fn); end
CP = dfdet_t3('path',fn,D,c,opts);
FF = @(y) ff(CP,y,c0);
q = zeros(size(p));
for k = 1:numel(p), q(k) = dfdet_newton(FF,p(k),x0,1,tol,100,30); end
end

function v = ff(CP,y,c0)
[F,f] = dfdet_t3('eval',CP,y,c0);
v = [F f];
end
