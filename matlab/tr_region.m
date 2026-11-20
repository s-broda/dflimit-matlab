function [P,info] = tr_region(bnd,c,opts)
%TR_REGION  The engine of LTRATIO_REGION: P{U <= g(V)} for a general boundary, case DK, zero start.
%
%   [P,info] = TR_REGION(bnd,c,opts): the conditional integral of TR_LIMDIST (the first part of the conditional
%   integral representation, the probability of a region of the plane, in the form with the inner probability of R = U + 1/2 given V = b in closed form) with the boundary a(b) = 1/2 + g(b) in place
%   of the parabola of the t-ratio, the zeros of a(b) found numerically and kinks of g (bnd.brks) as breakpoints;
%   bnd is a structure as in LTRATIO_REGION (type 'fun', 'line', 'parabola', 'wap'), g = +Inf on (0, bnd.v0), whose
%   probability P(V < v0) is added (TR_GILPELAEZ).  The C engine has a drop-in of this function (c/mex) for the
%   types other than 'fun'.  A port of Pregion of julia/DFExact.jl (tr_region there).
%
%   See also LTRATIO_REGION, TR_LIMDIST.
if nargin < 3 || isempty(opts), opts = struct; end
g = region_g(bnd);
v0 = 0; if isfield(bnd,'v0') && ~isempty(bnd.v0), v0 = bnd.v0; end
brks = zeros(1,0); if isfield(bnd,'brks') && ~isempty(bnd.brks), brks = bnd.brks(:).'; end
o = opts; if isfield(o,'method'), o = rmfield(o,'method'); end
o.c0 = 0; o.g0 = 0; o.verified = true;                % (the applications keep the agreement rule)
o.bnd = struct('a',@(b) 1/2 + g(b),'v0',v0,'brks',brks);
[P,info] = tr_limdist(0,'DK',c,o,'cdf');
if v0 > 0, P = P + tr_gilpelaez(v0,NaN,'DK',c); end
end

function g = region_g(bnd)
% the boundary g(v) of the region {U <= g(V)}, elementwise in v
switch lower(bnd.type)
    case 'fun',      g = bnd.g;
    case 'line',     a = bnd.a; b0 = bnd.b; g = @(v) b0 + a*v;
    case 'parabola', x = bnd.x; c0 = 0; if isfield(bnd,'c0'), c0 = bnd.c0; end, g = @(v) x*sqrt(v) + c0*v;
    case 'wap',      k = bnd.k; w = bnd.omega; g = @(v) wapbnd(v,k,w);
        if ~(k > 0 && w > 0), error('ltratio:args','ltratio_region: the weighted-average-power boundary needs k > 0 and omega > 0.'); end
    otherwise, error('tr_region:args','tr_region: unknown boundary type %s.',bnd.type);
end
end

function g = wapbnd(v,k,w)
% the boundary of the weighted-average-power test (1 + w^2 V)^(-1/2) erfcx(U/sqrt(2 (V + 1/w^2))) >= k: g(v) =
% z (v + 1/w^2)^(1/2), z the root of log erfcx(z/sqrt 2) = log k + log(1 + w^2 v)/2 (bisection on [-60, 60] to 1e-11;
% z = 60 when the left side is above the right at 60), as g_wap in julia/DFExact.jl
sz = size(v); v = v(:);
lr = log(k) + log(1 + w^2*v)/2;
f = @(z) lerfcx(z/sqrt(2)) - lr;
lo = -60*ones(size(v)); hi = 60*ones(size(v));
flo = f(lo); top = f(hi) > 0;
for it = 1:200
    m = (lo + hi)/2; fm = f(m);
    s = fm.*flo < 0; hi(s) = m(s); lo(~s) = m(~s); flo(~s) = fm(~s);
    if all(hi - lo < 1e-11), break, end
end
z = (lo + hi)/2; z(top) = 60;
g = reshape(z.*sqrt(v + 1/w^2),sz);
end
function y = lerfcx(x)
% log erfcx(x) without overflow
y = zeros(size(x)); n = x < 0;
y(n) = x(n).^2 + log(erfc(x(n))); y(~n) = log(erfcx(x(~n)));
end
