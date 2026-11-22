function [S,h] = dfdet_sigma(P,theta,nq)
%DFDET_SIGMA  The covariance Sigma(theta) of the functionals ell under the Ornstein-Uhlenbeck law with parameter -theta.
%
%   [S,h] = DFDET_SIGMA(P,theta,nq): Sigma(theta) = int_0^1 F(u) F(u)' du
%   (bilinear, no conjugation) for the prepared design P (DFDET_PREP), F the
%   transfer functions of ell (ell_k = int_0^1 F_k(u) dW(u)): e^{-theta (1 - u)}
%   for J(1), 1{u < lambda} e^{-theta (lambda - u)} for J(lambda),
%   int_u^1 g(r) e^{-theta (r - u)} dr for z.  On piece j = [a, b] and with
%   y = b - u, F(u) = sum_k g^(k)(u)/k! y^(k+1) e_k(theta y) + e^{-theta y} F(b)
%   (exact), e_k(z) = int_0^1 x^k e^{-z x} dx (the power series for |z| < 2, the
%   upward recurrence for k <= |z| and the downward one from the series in
%   1 - x for k > |z|, no cancellation for small theta), and F(b) by a backward
%   recursion over the pieces.  The integral over y by Gauss-Legendre panels
%   (nq nodes, default 24) graded from y = 0: first length 1/|theta|, doubling,
%   at most 24/|theta| while |e^{-2 theta y}| > e^{-80}.  Var J(1) in closed
%   form.  h = F(0): the functionals applied to the path e^{-theta r} (the
%   fixed and the stationary start).  All nodes of one theta are evaluated
%   together.
%
%   See also DFDET_PREP, DFDET_MGF.
if nargin < 3 || isempty(nq), nq = 24; end
[zq,wq] = glrule(nq);
N = P.N; p = P.p; nj = P.nj; br = P.br; np = numel(br) - 1; Dg = P.deg;
Ib = complex(zeros(max(p,1),np));
if p > 0
    for m = 1:np
        am = br(m); Lm = br(m+1) - am;
        e = ek(theta*Lm,Dg);                              % 1 x (Dg+1)
        a = taylorc(P.g{m},am - P.off(m));                % (Dg+1) x p
        Ib(:,m) = ((Lm.^(1:Dg+1)).*e)*a;
    end
end
Rb = complex(zeros(N,np));
for j = np:-1:1
    b = br(j+1);
    Rb(1,j) = exp(-theta*(1 - b));
    for i = 1:nj
        lam = P.jl(i);
        if lam >= b - 1e-14, Rb(1+i,j) = exp(-theta*(lam - b)); end
    end
    if p > 0 && j < np
        el = exp(-theta*(br(j+2) - br(j+1)));
        Rb(2+nj:N,j) = Ib(:,j+1) + el*Rb(2+nj:N,j+1);
    end
end
at = abs(theta); rt = real(theta);
if at > 0, h0 = 1/at; else, h0 = Inf; end
Fc = cell(np,1); Wc = cell(np,1);
for j = 1:np
    b = br(j+1); Lj = b - br(j);
    % the panels in y
    pa = zeros(0,1); pb = pa; y0 = 0;
    while y0 < Lj*(1 - 1e-14)
        if at*Lj <= 2
            y1 = Lj;
        else
            wdt = max(y0,h0);
            if 2*rt*y0 < 80, wdt = min(wdt,24*h0); end
            y1 = min(y0 + wdt,Lj);
            if (Lj - y1) < 0.5*wdt, y1 = Lj; end
        end
        pa(end+1,1) = y0; pb(end+1,1) = y1; %#ok<AGROW>
        y0 = y1;
    end
    y = reshape(pa' + (pb - pa)'.*zq,[],1); w = reshape((pb - pa)'.*wq,[],1);
    ey = exp(-theta*y);
    F = complex(zeros(numel(y),N));
    F(:,1:1+nj) = ey.*Rb(1:1+nj,j).';
    if p > 0
        E = ek(theta*y,Dg);                               % numel(y) x (Dg+1)
        EY = E.*(y.^(1:Dg+1));
        a = taylorv(P.g{j},b - y - P.off(j));             % numel(y) x (Dg+1) x p
        for i = 1:p
            F(:,1+nj+i) = sum(a(:,:,i).*EY,2) + ey*Rb(1+nj+i,j);
        end
    end
    Fc{j} = F; Wc{j} = w;
end
F = vertcat(Fc{:}); w = vertcat(Wc{:});
S = F.'*(F.*w);
S = (S + S.')/2;
S(1,1) = e0(2*theta);
if nargout > 1
    h = complex(zeros(N,1));
    h(1) = exp(-theta);
    for i = 1:nj, h(1+i) = exp(-theta*P.jl(i)); end
    if p > 0
        el = exp(-theta*(br(2) - br(1)));
        h(2+nj:N) = Ib(:,1) + el*Rb(2+nj:N,1);
    end
end
end

function v = e0(z)
if abs(z) < 1e-3, v = 1 - z/2 + z*z/6 - z^3/24; else, v = (1 - exp(-z))/z; end
end

function a = taylorc(C,s)
% the Taylor coefficients at the point s of the polynomials in the columns of C (coefficients of r^k): a(k+1,i) =
% g_i^(k)(s)/k!, by repeated synthetic division
a = C; D = size(C,1) - 1;
for k = 0:D-1
    for j = D-1:-1:k
        a(j+1,:) = a(j+1,:) + s*a(j+2,:);
    end
end
end

function a = taylorv(C,s)
% the same at the points s (a column): a(q,k+1,i) = g_i^(k)(s_q)/k!
D = size(C,1) - 1; p = size(C,2);
a = repmat(reshape(C,1,D + 1,p),numel(s),1,1);
for k = 0:D-1
    for j = D-1:-1:k
        a(:,j+1,:) = a(:,j+1,:) + s.*a(:,j+2,:);
    end
end
end

function E = ek(z,D)
% e_k(z) = int_0^1 x^k e^{-zx} dx, k = 0..D, at the points z (rows): for |z| < 2 the power series (terms below 1e-18 at the
% end); for |z| >= max(D, 2) the upward recurrence e_k = (k e_{k-1} - e^{-z})/z (stable for |z| >= k); in between the
% upward recurrence for k <= |z| and, for k > |z|, the downward recurrence e_{k-1} = (z e_k + e^{-z})/k (stable for
% k >= |z|) from e_D = e^{-z} sum_m z^m D!/(D + m + 1)! (the series in 1 - x, whose terms decrease from the first for
% |z| < D + 2)
z = z(:); nz = numel(z);
E = complex(zeros(nz,D+1));
a2 = real(z).^2 + imag(z).^2;
sm = a2 < 4;
if any(sm)
    zs = z(sm);
    Es = repmat(1./(1:D+1),numel(zs),1);
    term = ones(numel(zs),1);
    for n = 1:300
        term = term.*(-zs*(1/n));
        Es = Es + term.*(1./(n + (1:D+1)));
        if all(real(term).^2 + imag(term).^2 < 1e-36), break, end
    end
    E(sm,:) = Es;
end
lg = ~sm;
if any(lg)
    zl = z(lg); em = exp(-zl); iz = 1./zl; al = a2(lg);
    El = complex(zeros(numel(zl),D+1));
    El(:,1) = (1 - em).*iz;
    for k = 1:D
        El(:,k+1) = (k*El(:,k) - em).*iz;
    end
    K1 = D*ones(numel(zl),1);
    md = al < D^2; K1(md) = floor(sqrt(al(md)));
    mid = K1 < D;
    if any(mid)
        zm = zl(mid); emm = em(mid); k1 = K1(mid);
        term = complex(1/(D + 1))*ones(numel(zm),1); acc = term;
        for m = 1:5000
            term = term.*zm/(D + m + 1);
            acc = acc + term;
            if all(real(term).^2 + imag(term).^2 < 1e-36*(real(acc).^2 + imag(acc).^2)), break, end
        end
        Ed = complex(zeros(numel(zm),D+1));
        Ed(:,D+1) = emm.*acc;
        for k = D:-1:min(k1)+2
            Ed(:,k) = (zm.*Ed(:,k+1) + emm)/k;
        end
        use = (1:D+1) >= k1 + 2;
        Em = El(mid,:); Em(use) = Ed(use);
        El(mid,:) = Em;
    end
    E(lg,:) = El;
end
end
