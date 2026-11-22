function V = dfdet_peval(C,r)
%DFDET_PEVAL  Values of polynomials given by their coefficients (Horner's rule).
%   V = DFDET_PEVAL(C,r): C the (deg+1) x K matrix of the coefficients of r^k (rows k = 0..deg) of K polynomials, r a
%   vector of points; V(i,a) is polynomial a at r(i) (numel(r) x K).  Shared by the general engine (DFDET_PREP,
%   DFDET_SIGMA).
r = r(:); n = size(C,1);
V = repmat(C(n,:),numel(r),1);
for k = n-1:-1:1
    V = V.*r + C(k,:);
end
end
