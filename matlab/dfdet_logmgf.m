function [L0,LD] = dfdet_logmgf(Mo,c,s,t)
%DFDET_LOGMGF  The logarithm of the joint moment generating function of the general engine at points (s, t).
%
%   [L0,LD] = DFDET_LOGMGF(Mo,c,s,t): log M(s,t) = L0 - LD/2 at the points s
%   (a vector) and one t, with LD the principal logarithm of the determinant
%   det(I - 2 (s A + t B + gamma E11) Sigma(theta)) (no continuation: the caller
%   continues its branch) and L0 = -s/2 - (c + theta)/2, plus g0^2 (cf + E1(s) -
%   kappa s) with a fixed start (single-valued); the stationary start with the
%   extended matrices.  Sigma(theta) and the data of the start are computed
%   once for the t.  (logmgf_point of DFDet.jl.)
%
%   See also DFDET_MGF, DFDET_MOM.
P = Mo.P; s = s(:);
theta = sqrt(complex(c*c - 2*t)); gam = (c + theta)/2;
L0 = -s/2 - (c + theta)/2; LD = complex(zeros(size(s)));
if strcmp(Mo.start,'stationary')
    [S,h] = dfdet_sigma(P,theta);
    v = -1/(2*c); n = P.N; N = n + 1;
    Sx = complex(zeros(N)); Sx(1:n,1:n) = S + v*(h*h.'); Sx(N,1:n) = v*h.'; Sx(1:n,N) = v*h; Sx(N,N) = v;
    for k = 1:numel(s)
        Q = s(k)*Mo.A + t*Mo.B; Q(1,1) = Q(1,1) + gam; Q(N,N) = Q(N,N) - gam;
        LD(k) = ldet(eye(N) - 2*(Q*Sx));
    end
    return
end
if strcmp(Mo.start,'fixed')
    [~,~,~,fx] = dfdet_mgf(Mo,c,t);
    L0 = L0 + Mo.g0^2*(fx.cf + dfdet_fixe1(fx,s) - Mo.kappa*s);
end
S = dfdet_sigma(P,theta); N = P.N;
if abs(theta) < 0.5
    z = 2*theta;
    if abs(z) < 1e-3, e = 1 - z/2 + z*z/6 - z^3/24; else, e = (1 - exp(-z))/z; end
    m11 = 1 - (c + theta)*e;
else
    if c > 0, d = -2*t/(theta + c); else, d = theta - c; end
    m11 = (d + (c + theta)*exp(-2*theta))/(2*theta);
end
% (the (1,1) entry m11 - 2 (t B Sigma)_11 from these two terms alone, and the term of s added to it: above the unit
% root it is of the order e^{-2c}, and the form 1 - 2 (t B Sigma)_11 - 2 gamma Sigma_11 + (m11 - (1 - 2 gamma
% Sigma_11)) carried an absolute rounding error of 1e-16; dfdet_mgf's buildM)
AS = P.A*S; Q = -2*(t*P.B*S); M0 = eye(N) + Q;
M0(1,:) = M0(1,:) - 2*gam*S(1,:);
M0(1,1) = m11 + Q(1,1);
for k = 1:numel(s)
    LD(k) = ldet(M0 - 2*s(k)*AS);
end
end

function l = ldet(M)
% log det M from the LU factors (log|det| plus the principal logarithm of the phase)
[~,U,p] = lu(M,'vector');
u = diag(U); a = abs(u);
sg = prod(u./a);
n = numel(p); seen = false(1,n); odd = 0;
for i = 1:n
    if ~seen(i)
        j = i; len = 0;
        while ~seen(j), seen(j) = true; j = p(j); len = len + 1; end
        odd = odd + len - 1;
    end
end
if mod(odd,2) == 1, sg = -sg; end
l = sum(log(a)) + log(sg);
end
