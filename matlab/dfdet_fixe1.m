function E = dfdet_fixe1(fx,s,g)
%DFDET_FIXE1  The bounded part E1(s) of the exponent of the fixed start, at many points s.
%
%   E = DFDET_FIXE1(fx,s,g): E1(s) = s [p1' (I - sH)^{-1} w1 + p2' (I - sH)^{-1}
%   w2] for the data fx of the fixed start (DFDET_MGF) at the points s, by
%   Gaussian elimination with partial pivoting (pivots by |.|^2) on
%   (I - sH) x = s w (|s| <= 1) or (I/s - H) x = w (|s| > 1), so that nothing
%   grows at large |s| (fixE1 of DFDet.jl).  The fields of fx may be stacked
%   over several outer nodes (H r x r x nn, the vectors r x nn); g gives the
%   node of each s (default 1).  All points are eliminated together.
%
%   See also DFDET_MGF.
s = s(:); ns = numel(s);
r = size(fx.H,1);
if r == 0 || ns == 0, E = complex(zeros(ns,1)); return, end
if nargin < 3 || isempty(g), g = ones(ns,1); end
H = reshape(permute(fx.H,[3 1 2]),[],r*r); H = reshape(H(g,:),ns,r,r);       % ns x r x r
w1 = fx.w1(:,g).'; w2 = fx.w2(:,g).'; p1 = fx.p1(:,g).'; p2 = fx.p2(:,g).';  % ns x r
big = real(s).^2 + imag(s).^2 > 1;
d = ones(ns,1); d(big) = 1./s(big);
f = s; f(big) = 1;
sc = s; sc(big) = 1;
B = complex(zeros(ns,r,r + 2));
B(:,:,1:r) = -sc.*H;
for j = 1:r, B(:,j,j) = B(:,j,j) + d; end
B(:,:,r+1) = f.*w1; B(:,:,r+2) = f.*w2;
rows = (1:ns)';
for k = 1:r
    % the pivot of each point among the rows k..r
    a2 = real(B(:,k:r,k)).^2 + imag(B(:,k:r,k)).^2;
    [~,ip] = max(a2,[],2); piv = ip + k - 1;
    sw = find(piv ~= k);
    if ~isempty(sw)
        for j = k:r+2
            i1 = rows(sw) + (k - 1)*ns + (j - 1)*ns*r; i2 = rows(sw) + (piv(sw) - 1)*ns + (j - 1)*ns*r;
            tmp = B(i1); B(i1) = B(i2); B(i2) = tmp;
        end
    end
    ipv = 1./B(:,k,k); B(:,k,k) = ipv;
    for i = k+1:r
        l = B(:,i,k).*ipv;
        B(:,i,k+1:r+2) = B(:,i,k+1:r+2) - l.*B(:,k,k+1:r+2);
    end
end
E = complex(zeros(ns,1));
for col = [r + 1, r + 2]
    if col == r + 1, pv = p1; else, pv = p2; end
    for k = r:-1:1
        x = B(:,k,col);
        for j = k+1:r, x = x - B(:,k,j).*B(:,j,col); end
        B(:,k,col) = x.*B(:,k,k);
        E = E + pv(:,k).*B(:,k,col);
    end
end
end
