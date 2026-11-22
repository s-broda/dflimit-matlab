function x = dfdet_newton(FF,p,x0,step,tol,maxit,clampit)
%DFDET_NEWTON  The solution of F(x) = p for an increasing F by Newton's method (the quantiles of the general engine).
%
%   x = DFDET_NEWTON(FF,p,x0,step,tol,maxit,clampit): FF(x) returns [F(x), f(x)] from one evaluation.  With clampit > 0
%   first up to clampit Newton steps from x0, each limited to step; where that does not converge to |dx| < tol (or
%   leaves the finite numbers), the safeguarded method of newton_safe of DFDet.jl: a bracket grown from
%   [x0 - step, x0 + step] by doubling, then Newton's method kept inside the bracket (bisection where a step leaves it),
%   to |dx| < tol, at most maxit steps (the warning dfdet:newton without convergence).
if nargin < 5 || isempty(tol), tol = 1e-10; end
if nargin < 6 || isempty(maxit), maxit = 100; end
if nargin < 7 || isempty(clampit), clampit = 0; end
if ~(isnumeric(x0) && isscalar(x0) && isreal(x0) && isfinite(x0)), error('dfdet:args','dfdet: the start value of the search (opts.x0, y0 or b0) must be a finite real scalar.'); end
if ~(isnumeric(tol) && isscalar(tol) && isreal(tol) && tol > 0 && isfinite(tol)), error('dfdet:args','dfdet: opts.tol must be a positive finite number.'); end
x = x0;
for it = 1:clampit
    Ff = FF(x);
    dx = min(max((Ff(1) - p)/Ff(2),-step),step);
    if ~isfinite(dx), break, end
    x = x - dx;
    if abs(dx) < tol, return, end
end
x = x0; lo = x - step; hi = x + step; w = step;
for it = 1:60
    Ff = FF(lo); if ~(Ff(1) > p), break, end
    lo = lo - w; w = 2*w;
end
w = step;
for it = 1:60
    Ff = FF(hi); if ~(Ff(1) < p), break, end
    hi = hi + w; w = 2*w;
end
for it = 1:maxit
    Ff = FF(x); F = Ff(1); f = Ff(2);
    if F > p, hi = x; else, lo = x; end
    xn = x - (F - p)/f;
    if ~(lo <= xn && xn <= hi && isfinite(xn)), xn = (lo + hi)/2; end   % (a step below half an ulp lands on the end of the bracket: converged)
    if abs(xn - x) < tol, x = xn; return, end
    x = xn;
end
warning('dfdet:newton','dfdet: newton_safe: no convergence (p = %g, x = %g).',p,x);
end
