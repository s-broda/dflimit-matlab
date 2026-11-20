function [q,info] = adfmax_icdf(p,cas,c,opts)
%ADFMAX_ICDF  Quantiles of the limiting ADFmax statistic of Leybourne (1995).
%
%   q = ADFMAX_ICDF(p,cas,c,opts): the p-quantiles (p an array in (0, 1)) of the limit of ADFMAX_CDF, by Newton's
%   method with ADFMAX_PDF as the derivative, safeguarded by bisection on a bracket grown from [x0 - 1, x0 + 1]
%   (x0 = opts.x0, default -2) to opts.tol (1e-10), as adfmax_icdf of julia/DFExact.jl.  info: F at q.
%
%   Example: adfmax_icdf([0.01 0.05 0.10],'DT')   % -3.67949 -3.11891 -2.82890 (the quantiles of the ADFmax statistic in the paper)
%
%   See also ADFMAX_CDF, ADFMAX_PDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','adfmax_icdf: p must lie in (0,1).'); end
tol = 1e-10; x0 = -2; o = opts;
if isfield(o,'tol'), if ~isempty(o.tol), tol = o.tol; end, o = rmfield(o,'tol'); end
if isfield(o,'x0'), if ~isempty(o.x0), x0 = o.x0; end, o = rmfield(o,'x0'); end
q = zeros(size(p)); Fq = q;
for k = 1:numel(p)
    pk = p(k); x = x0; lo = x - 1; hi = x + 1;
    while adfmax_cdf(lo,cas,c,o) > pk, lo = lo - 1; end
    while adfmax_cdf(hi,cas,c,o) < pk, hi = hi + 1; end
    x = (lo + hi)/2; conv = false;
    for it = 1:100
        F = adfmax_cdf(x,cas,c,o); f = adfmax_pdf(x,cas,c,o);
        if F > pk, hi = x; else, lo = x; end
        xn = x - (F - pk)/f;
        if ~(lo < xn && xn < hi) || ~isfinite(xn), xn = (lo + hi)/2; end
        if abs(xn - x) < tol, x = xn; conv = true; break, end
        x = xn;
    end
    if ~conv, warning('ltratio:icdf','adfmax_icdf: no convergence at p = %g.',pk); end
    q(k) = x; Fq(k) = adfmax_cdf(x,cas,c,o);
end
info = struct('F',Fq);
end
