function [q,err] = pt_quantile(p,D,cbar,opts)
%PT_QUANTILE  Critical values of the point-optimal statistic against cbar for the design D.
%
%   [q,err] = PT_QUANTILE(p,D,cbar,opts): the p-quantile at c = 0 of the limit
%   of the point-optimal statistic (PT_CDF), the critical value of the level-p
%   test (it rejects for small values).  By default the Illinois variant of
%   regula falsi on the single inversion of Bao and Kan (PT_CDF), from the
%   bracket E Q -/+ its standard deviation (grown as needed) to the width 1e-12
%   times it, with the transform values along the ray kept across the search
%   (DFDET_KAN; err the error estimate of the distribution function at q);
%   with opts.method = 'grid' (or opts.N, opts.umax) Newton's method
%   safeguarded by bisection on [1e-3, 400] from opts.q0 (6) to |F - p| <
%   1e-12 on the fixed grid (err NaN).
%
%   Example: pt_quantile(0.05, dfdet_perron('B', 0.1), -16)   % 6.0936
%
%   See also PT_CDF, PT_POWER.
if nargin < 4 || isempty(opts), opts = struct; end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','pt_quantile: p must lie in (0, 1).'); end
[grid,tol] = dfdet_ptmethod(opts,0,'pt_quantile');
q = zeros(size(p)); err = NaN(size(p));
if grid
    G = dfdet_qgrid(D,double(cbar),0,getf(opts,'N',60000),getf(opts,'umax',4000));
    for k = 1:numel(p), q(k) = ptquantile(G,p(k),getf(opts,'q0',6)); end
    return
end
% (regula falsi on the single inversion of Bao and Kan, the transform values along the ray kept across the search)
K = dfdet_kan('cf',dfdet_model(dfdet_prep(D,'pt',double(cbar))),0);
for k = 1:numel(p), [q(k),err(k)] = dfdet_kan('pt_quantile1',K,double(cbar),p(k),tol); end
end

function q = ptquantile(G,p,q0)
lo = 1e-3; hi = 400; q = q0;
for it = 1:100
    [F,f] = dfdet_ptf(G,q);
    if F > p, hi = q; else, lo = q; end
    if abs(F - p) < 1e-12, break, end
    qn = q - (F - p)/f;
    if lo < qn && qn < hi, q = qn; else, q = (lo + hi)/2; end
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
