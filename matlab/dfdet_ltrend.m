function D = dfdet_ltrend(p)
%DFDET_LTREND  The polynomial trend of order p in the orthonormal (shifted Legendre) basis.
%
%   D = DFDET_LTREND(p): the regressors sqrt(2k + 1) P_k(2r - 1), k = 0..p, the
%   same span as DFDET_TREND(p), hence the same statistic, with the Gram
%   matrix I instead of a Hilbert matrix (condition 5e11 at p = 8, where the
%   eigenvalue threshold of the preparation would drop directions of the
%   monomial basis).  A function design (DFDET_FDESIGN).
%
%   Example: dfdet_cdf(-4, dfdet_ltrend(8))
%
%   See also DFDET_TREND, DFDET_FDESIGN.
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 0), error('dfdet:args','dfdet_ltrend: p must be an integer >= 0.'); end
D = dfdet_fdesign(@(r) leg(r,p),[],struct('K',p + 1,'name',sprintf('polynomial trend of order %d (Legendre basis)',p)));
end

function v = leg(r,p)
x = 2*r - 1; P = zeros(p + 1,1); P(1) = 1;
if p >= 1, P(2) = x; end
for k = 1:p-1, P(k+2) = ((2*k + 1)*x*P(k+1) - k*P(k))/(k + 1); end
v = sqrt(2*(0:p)' + 1).*P;
end
