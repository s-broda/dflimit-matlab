function [x,F,it,conv] = ltratio_newton(Ff,ff,x,sd,tol,p)
%LTRATIO_NEWTON  The p-quantile of a distribution function Ff with density ff by Newton's method, safeguarded by
%   bisection on a bracket grown from x -+ sd; shared by LTRATIO_ICDF and DFX_QUANTILE.  Returns the quantile, Ff at
%   it, the iteration count and a convergence flag (at most 60 iterations).
lo = x - sd; hi = x + sd;
while Ff(lo) > p, lo = lo - 2*(x - lo); end
while Ff(hi) < p, hi = hi + 2*(hi - x); end
conv = false;
for it = 1:60
    F = Ff(x); f = ff(x);
    if F > p, hi = x; else, lo = x; end
    xn = x - (F - p)/f;
    if ~(lo < xn && xn < hi) || ~isfinite(xn), xn = (lo + hi)/2; end
    if abs(xn - x) < tol, x = xn; conv = true; break; end
    x = xn;
end
F = Ff(x);
end
