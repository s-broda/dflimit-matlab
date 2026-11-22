function q = pt_quantile(p,D,cbar,opts)
%PT_QUANTILE  Critical values of the point-optimal statistic against cbar for the design D.
%
%   q = PT_QUANTILE(p,D,cbar,opts): the p-quantile at c = 0 of the limit of the
%   point-optimal statistic (PT_CDF), the critical value of the level-p test
%   (it rejects for small values), by Newton's method safeguarded by
%   bisection on [1e-3, 400] from opts.q0 (6) to |F - p| < 1e-12.
%
%   Example: pt_quantile(0.05, dfdet_perron('B', 0.1), -16)   % 6.0936
%
%   See also PT_CDF, PT_POWER.
if nargin < 4 || isempty(opts), opts = struct; end
G = dfdet_qgrid(D,double(cbar),0,getf(opts,'N',60000),getf(opts,'umax',4000));
q = zeros(size(p));
for k = 1:numel(p), q(k) = ptquantile(G,p(k),getf(opts,'q0',6)); end
end

function q = ptquantile(G,p,q0)
lo = 1e-3; hi = 400; q = q0;
for it = 1:100
    z = exp(G.lM - 1i*G.u*q);
    F = 0.5 - sum(imag(z)./G.u.*G.w)/pi; f = sum(real(z).*G.w)/pi;
    if F > p, hi = q; else, lo = q; end
    if abs(F - p) < 1e-12, break, end
    qn = q - (F - p)/f;
    if lo < qn && qn < hi, q = qn; else, q = (lo + hi)/2; end
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
