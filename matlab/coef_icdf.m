function y = coef_icdf(p,cas,c,opts)
%COEF_ICDF  Quantiles of the limiting Dickey-Fuller coefficient statistic.
%
%   y = COEF_ICDF(p,cas,c,opts): the p-quantiles (p an array in (0, 1)) of the limit of n(phi^ - phi0) (COEF_CDF),
%   by the Illinois variant of regula falsi on a bracket grown from [-40, 10] to the width opts.tol (1e-9), about
%   ten evaluations of COEF_CDF.  Options as in COEF_CDF.
%
%   Example: coef_icdf(0.05,'DK')         % -8.039, the 5% critical value of rho^
%
%   See also COEF_CDF, LTRATIO_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','coef_icdf: p must lie in (0,1).'); end
tol = 1e-9; o = opts; if isfield(o,'tol'), if ~isempty(o.tol), tol = o.tol; end, o = rmfield(o,'tol'); end
y = zeros(size(p));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (each brittle regime warned once)
for k = 1:numel(p)
    f = @(x) coef_cdf(x,cas,c,o) - p(k);
    lo = -40; hi = 10; flo = f(lo); fhi = f(hi);
    while flo > 0, hi = lo; fhi = flo; lo = 2*lo - 10; flo = f(lo); end
    while fhi < 0, lo = hi; flo = fhi; hi = 2*hi + 10; fhi = f(hi); end
    side = 0; x = (lo + hi)/2;
    for it = 1:200
        x = (lo*fhi - hi*flo)/(fhi - flo); fx = f(x);
        if fx == 0, break, end
        if fx < 0, lo = x; flo = fx; if side == -1, fhi = fhi/2; end, side = -1;
        else, hi = x; fhi = fx; if side == 1, flo = flo/2; end, side = 1; end
        if hi - lo < tol, x = (lo*fhi - hi*flo)/(fhi - flo); break, end
    end
    y(k) = x;
end
end
