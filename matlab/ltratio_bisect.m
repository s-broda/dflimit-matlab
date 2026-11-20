function x = ltratio_bisect(f,lo,hi,tol)
%LTRATIO_BISECT  The root of f on the bracket [lo, hi] (f(lo) f(hi) < 0) by bisection to the width tol (1e-9): the
%   same steps as ltratio_bisect of julia/DFExact.jl, used by POWER_ENVELOPE and by the replication scripts
%   (public/replication/ltratio_paper.m of the paper's repository).
if nargin < 4 || isempty(tol), tol = 1e-9; end
flo = f(lo); fhi = f(hi);
if ~(flo*fhi < 0), error('ltratio:bisect','ltratio_bisect: no sign change on [%g, %g] (%g, %g).',lo,hi,flo,fhi); end
for it = 1:200
    m = (lo + hi)/2; fm = f(m);
    if fm*flo < 0, hi = m; else, lo = m; flo = fm; end
    if hi - lo < tol, break, end
end
x = (lo + hi)/2;
end
