function [ts,mult] = dfdet_tstar(Mo,c)
%DFDET_TSTAR  t_* = sup{t: E e^{tV} < Inf}, the first zero of D_V(t) on the positive real axis.
%
%   [ts,mult] = DFDET_TSTAR(Mo,c) for a model Mo (DFDET_MODEL): t_* =
%   1/(2 lambda_1), lambda_1 the largest eigenvalue of the covariance
%   operator of the quadratic form V, estimated by the Nystrom method
%   (Gauss-Legendre with 16 nodes on each of 8 equal panels of every piece;
%   with GLS detrending and for the LM statistic, where V is not the norm of
%   a projection of J, the operator of the quadratic form with 48 nodes); D_V is then sampled on 81
%   points of +-2% around the estimate: at the first sign change bisection to
%   1e-13 relative, and where there is none (a double zero: lambda_1 =
%   lambda_2, as in model C at lambda = 1/2) the minimum of D_V by golden
%   section (mult = 2).  Far from the unit root, where the estimate is off by
%   more than 2% (the kernel of the covariance has the width 1/|c|, which the
%   rule does not resolve), D_V is scanned from the branch point c^2/2 in
%   omega (t = (c^2 + omega^2)/2) for its first sign change or a double zero.
%
%   See also DFDET_MGF, DFDET_CPATH.
if isfield(Mo,'vm'), [ts,mult] = vol_tstar(Mo); return, end
if strcmp(Mo.start,'stationary') && c >= 0, error('dfdet:args','%s',dfdet_msg('stat_c')); end
f = @(t) real(exp(ldv(Mo,c,t)));
P = Mo.P;
% (eigtop projects J orthogonally off the span of the design; the V of GLS detrending and of the LM statistic (detrend
% 'lm') is no such norm, and eigtop_q takes its quadratic form: with eigtop the LM statistic's t_* came out up to 2.3 times
% too large, above the unit root beyond twice the edge, where the path then failed)
if isnan(P.cbar) && ~strcmp(P.detrend,'lm'), [l1,l2] = eigtop(P,c,16,strcmp(Mo.start,'stationary')); else, [l1,l2] = eigtop_q(P,c,48); end
te = 1/(2*l1);
grid = linspace(0.98*te,1.02*te,81);
vals = zeros(size(grid));
for k = 1:numel(grid), vals(k) = f(grid(k)); end
k = find(vals <= 0,1);
if ~isempty(k) && k > 1
    ts = bisect(f,grid(k-1),grid(k)); mult = 1 + (abs(l1 - l2) < 1e-4*l1);
    return
end
% far from the unit root the estimate is off by more than 2%: t_* then lies beyond the branch point c^2/2 (theta = i
% omega imaginary, t = (c^2 + omega^2)/2), and D_V is scanned in omega from the branch point for its first sign change
% (at most 4000 steps of 0.05), or for a double zero (a local minimum below 1e-6 times the largest value before it)
found = false;
if c < 0 && f(c*c/2 + 5e-7) > 0
    tw = @(om) (c*c + om*om)/2;
    fprev = f(tw(1e-3)); fmax = abs(fprev); fpp = Inf;
    for j = 1:4000
        om = 0.05*j; fj = f(tw(om));
        if fj <= 0, ts = bisect(f,tw(om - 0.05),tw(om)); mult = 1; return, end
        if fprev < fpp && fprev <= fj && fprev < 1e-6*fmax
            a = tw(max(om - 0.1,1e-3)); b = tw(om); found = true; break
        end
        fmax = max(fmax,fj); fpp = fprev; fprev = fj;
    end
end
if ~found
    [~,k] = min(vals);
    a = grid(max(k - 1,1)); b = grid(min(k + 1,numel(grid)));
end
g = (sqrt(5) - 1)/2;
x1 = b - g*(b - a); x2 = a + g*(b - a); f1 = f(x1); f2 = f(x2);
for it = 1:80
    if f1 < f2
        b = x2; x2 = x1; f2 = f1; x1 = b - g*(b - a); f1 = f(x1);
    else
        a = x1; x1 = x2; f1 = f2; x2 = a + g*(b - a); f2 = f(x2);
    end
end
ts = (a + b)/2; mult = 2;
end

function l = ldv(Mo,c,t)
[~,l] = dfdet_mgf(Mo,c,complex(t),true);
end

function m = bisect(f,a,b)
for it = 1:100
    m = (a + b)/2;
    if f(m) > 0, a = m; else, b = m; end
    if b - a < 1e-13*b, break, end
end
m = (a + b)/2;
end

function K = oukernel(c,r,u,v0)
% the covariance Cov(J(r), J(u)) = e^{c(r+u)} (min(r,u) e_0(2c min(r,u)) + v0) of the Ornstein-Uhlenbeck process started at
% zero (v0 = 0) or with the stationary start's further term (v0 = -1/(2c)); far from the unit root, where e_0(2c min)
% overflows (-2c min(r,u) > 700), in the form (e^{c|r-u|} - e^{c(r+u)})/(-2c)
R = r + 0*u; U = u + 0*r; r = R; u = U;
mn = min(r,u); z = 2*c*mn;
e = (1 - exp(-z))./z; sm = abs(z) < 1e-3; zs = z(sm); e(sm) = 1 - zs/2 + zs.*zs/6 - zs.^3/24;
K = exp(c*(r + u)).*(mn.*e + v0);
far = -2*c*mn > 700;
if any(far(:))
    K(far) = (exp(c*abs(r(far) - u(far))) - exp(c*(r(far) + u(far))))/(-2*c) + exp(c*(r(far) + u(far)))*v0;
end
end

function [rs,ws,pj] = nystrom(br,m)
[z,w] = glrule(m);
rs = []; ws = []; pj = [];
for j = 1:numel(br)-1
    a = br(j); b = br(j+1);
    for kk = 0:7
        aa = a + (b - a)*kk/8;
        rs = [rs; aa + (b - a)/8*z]; ws = [ws; (b - a)/8*w]; pj = [pj; j*ones(m,1)]; %#ok<AGROW>
    end
end
end

function [l1,l2] = eigtop(P,c,m,stat,kern)
% the two largest eigenvalues of the covariance operator of the residual Jt = J - d' G^{-1} m on L2[0,1], by the
% Nystrom method on the kernel projected on the complement of span d in the discrete inner product (kern: another
% covariance kernel K(r, u) of a column r and a row u in place of the Ornstein-Uhlenbeck one, the volatility model's)
[rs,ws] = nystrom(P.br,m);
v0 = 0; if stat, v0 = -1/(2*c); end
sw = sqrt(ws);
if nargin < 5, Kw = sw.*oukernel(c,rs,rs',v0).*sw'; else, Kw = sw.*kern(rs,rs').*sw'; end
if P.K > 0
    D = P.D; X = zeros(numel(rs),size(D.C{1},2));
    for q = 1:numel(rs), X(q,:) = dval(D,rs(q)).'; end
    X = X.*sw;
    [Q,~] = qr(X,0);
    Kw = Kw - Q*(Q'*Kw); Kw = Kw - (Kw*Q)*Q';
end
ev = sort(eig((Kw + Kw')/2));
l1 = ev(end); l2 = ev(end-1);
end

function [l1,l2] = eigtop_q(P,c,m,kern)
% the two largest eigenvalues of the operator of the quadratic form V = int J^2 + ell' B ell for any prepared design
% (GLS detrending): the vector X = (J at the nodes, J(1), J(lambda_i)) with its covariance C (zero start), ell = Lambda X
% (z by the quadrature of int g J), V = X' (W + Lambda' B Lambda) X, and the eigenvalues of C^{1/2} (W + Lambda' B
% Lambda) C^{1/2} (kern: another covariance kernel K(r, u) of a column r and a row u in place of the Ornstein-Uhlenbeck
% one, the volatility model's)
[rs,ws,pj] = nystrom(P.br,m);
n = numel(rs);
pts = [rs; 1; P.jl(:)]; M = numel(pts);
if nargin < 4, C = oukernel(c,pts,pts',0); else, C = kern(pts,pts'); end
Lam = zeros(P.N,M); Lam(1,n+1) = 1;
for i = 1:P.nj, Lam(1+i,n+1+i) = 1; end
for q = 1:n
    j = pj(q);
    if P.p > 0, Lam(2+P.nj:end,q) = ws(q)*dfdet_peval(P.g{j},rs(q) - P.off(j)).'; end
end
Q = Lam'*P.B*Lam;
id = (0:n-1)*(M + 1) + 1; Q(id) = Q(id) + ws';
[V,E] = eig((C + C')/2); sq = V*diag(sqrt(max(diag(E),0)))*V';
Q = (Q + Q')/2; T = sq*Q*sq;
ev = sort(eig((T + T')/2));
l1 = ev(end); l2 = ev(end-1);
end

function v = dval(D,r)
np = numel(D.br) - 1; j = 1;
while j < np && r > D.br(j+1), j = j + 1; end
off = 0; if D.loc, off = D.br(j); end
v = dfdet_peval(D.C{j},r - off).';
end

function [ts,mult] = vol_tstar(VM)
% t_* of the volatility model (vol_tstar of DFDet_vol.jl): the Nystrom estimate te = 1/(2 lambda_1) with the covariance of
% the regimes (with the stationary start's term; with GLS detrending, where V is not the norm of a projection of the path,
% the operator of the quadratic form itself, eigtop_q, on the nodes of the least squares designs, since the scan from zero
% below covers an estimate that is off); D_V on 81 points of [0.98 te, 1.02 te], bisection at the first change of sign,
% the golden section for a double zero; where the estimate is off, D_V scanned from zero in steps te/1000 (at most 8000)
% for its first change of sign or a double zero
f = @(t) real(exp(ldv(VM,0,t)));
kern = @(r,u) dfdet_vol('kernel',VM,r,u);
if isnan(VM.P.cbar), [l1,l2] = eigtop(VM.P,min(VM.c),16,false,kern); else, [l1,l2] = eigtop_q(VM.P,min(VM.c),16,kern); end
te = 1/(2*l1);
grid = linspace(0.98*te,1.02*te,81);
vals = zeros(size(grid));
for k = 1:numel(grid), vals(k) = f(grid(k)); end
k = find(vals <= 0,1);
if ~isempty(k) && k > 1, ts = bisect(f,grid(k-1),grid(k)); mult = 1 + (abs(l1 - l2) < 1e-4*l1); return, end
if isempty(k)
    [~,j] = min(vals);
    if j > 1 && j < numel(grid) && vals(j) < 1e-6*max(vals), ts = golden(f,grid(j-1),grid(j+1)); mult = 2; return, end
end
h = te/1000; fprev = f(h); fmax = abs(fprev); fpp = Inf; tprev = h;
for j = 2:8000
    tj = j*h; fj = f(tj);
    if fj <= 0, ts = bisect(f,tprev,tj); mult = 1; return, end
    if fprev < fpp && fprev <= fj && fprev < 1e-6*fmax, ts = golden(f,tprev - h,tj); mult = 2; return, end
    fmax = max(fmax,fj); fpp = fprev; fprev = fj; tprev = tj;
end
error('dfdet:tstar','tstar: no zero of D_V of the volatility model below %g (the Nystrom estimate of t_* is %g)',8*te,te);
end

function m = golden(f,a,b)
g = (sqrt(5) - 1)/2;
x1 = b - g*(b - a); x2 = a + g*(b - a); f1 = f(x1); f2 = f(x2);
for it = 1:80
    if f1 < f2
        b = x2; x2 = x1; f2 = f1; x1 = b - g*(b - a); f1 = f(x1);
    else
        a = x1; x1 = x2; f1 = f2; x2 = a + g*(b - a); f2 = f(x2);
    end
end
m = (a + b)/2;
end
