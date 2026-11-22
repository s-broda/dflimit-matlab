function [P,a1] = dfdet_prep(D,kind,cbar)
%DFDET_PREP  The prepared design of the general engine: the functionals and the quadratic forms of U and V.
%
%   P = DFDET_PREP(D) for a design D (DFDET_DETERM and its constructors): the
%   orthonormal basis g of span{d, d'} (coefficients P.g{j}, (deg+1) x p, on
%   piece j), the jump points P.jl, the matrices A, B of U + 1/2 = ell' A ell
%   and V - int J^2 = ell' B ell for the Gaussian vector
%   ell = (J(1), J(lambda_i) at the jumps, z = int g J), and the nonzero
%   eigenvalues P.Dv (with eigenvectors P.L) of A: P.r = numel(P.Dv) roots in
%   the argument of U, P.rp positive and P.rm negative at real t (Sylvester's
%   law of inertia).  P.b is the coefficient of J(0) in U + 1/2 =
%   ell' A ell - J(0)^2/2 + J(0) b' ell (a fixed or a stationary start).  The
%   basis removes the linear dependences between int d J and int d' J, so that
%   the covariance of ell is nonsingular.  For the statistic that detrends
%   first (D.detrend = 'first') the block of z of A is sym((Pm Phi_n)(Pm Phi_d))
%   in place of sym(Pm Phi_n).
%
%   P = DFDET_PREP(D,'gls',cbar): the data GLS-detrended at abar = 1 + cbar/n
%   (cbar < 0, Elliott, Rothenberg and Stock): the first observation as it is,
%   the quasi-differences y_t - abar y_{t-1} for t >= 2, the same for the
%   regressors, and the t-ratio of the regression of Delta ytilde_t on
%   ytilde_{t-1} without deterministic terms.  Only the directions b with
%   d(0)' b = 0 and Delta_i' b = 0 are estimated at the rate that matters; with
%   delta = Z' d for a basis Z of that space, psi = delta' - cbar delta and
%   eta = int psi (dJ - cbar J dr), J_g = J - delta' Psi^{-1} eta,
%   U_g + 1/2 = J_g(1)^2/2 (rank one) and V_g = int J_g^2.  Zero start only.
%
%   P = DFDET_PREP(D,'pt',cbar): the limit Q = cbar^2 V_g + ell' A ell of the
%   point-optimal statistic against cbar (PT_CDF), A in the field A.
%   [P,a1] = DFDET_PREP(D,'lm'): the denominator V = int Sbar^2 of the LM
%   statistic (LM_CDF; the numerator a placeholder), and the row a1 with
%   J_0(1) = a1' ell.
%
%   See also DFDET_MODEL, DFDET_SIGMA, DFDET_CDF.
if nargin < 2 || isempty(kind), kind = 'ls'; end
if ~strcmp(D.detrend,'onestep') && ~strcmp(D.detrend,'first'), error('dfdet:args','the statistic of the design must be ''onestep'' or ''first''.'); end
a1 = [];
switch kind
    case 'ls', P = prep_ls(D);
    case 'gls'
        if strcmp(D.detrend,'first'), error('dfdet:args','%s',dfdet_msg('first_gls')); end
        if ~(cbar < 0 && isfinite(cbar)), error('dfdet:args','GLS detrending needs cbar < 0 (a finite number).'); end
        P = prep_gls(D,double(cbar));
    case 'pt'
        if ~(cbar < 0 && isfinite(cbar)), error('dfdet:args','the point-optimal statistic needs cbar < 0 (a finite number).'); end
        P = prep_pt(D,double(cbar));
    case 'lm', [P,a1] = prep_lm(D);
    otherwise, error('dfdet:args','dfdet_prep: kind must be ''ls'', ''gls'', ''pt'' or ''lm''.');
end
end

% ---------------------------------------------------------------------------------------------------------------------
function P = prep_ls(D)
K = size(D.C{1},2); np = numel(D.br) - 1; deg = size(D.C{1},1) - 1; br = D.br; off = offs(D);
[jl,Dl] = jumps(D); nj = numel(jl);
if K == 0
    g = repmat({zeros(deg + 1,0)},1,np); p = 0; G = zeros(0);
else
    F = cell(1,np); for j = 1:np, F{j} = [D.C{j}, dcoef(D.C{j})]; end
    Gf = gram(F,br,D.loc);
    [U,lam] = symeig(Gf);
    keep = lam > 1e-11*max(lam);
    U = U(:,keep); lam = lam(keep); p = numel(lam);
    g = cell(1,np); for j = 1:np, g{j} = F{j}*U*diag(1./sqrt(lam)); end
    Phi = U*diag(sqrt(lam));
    Phid = Phi(1:K,:); Phin = Phi(K+1:2*K,:);
    G = Gf(1:K,1:K);
    Pm = Phid'/G;                                       % p x K: m' G^{-1} = z' Pm
    if ~(max(max(abs(Phid*Phid' - G))) < 1e-10*max(abs(G(:)))), error('dfdet:prep','dfdet_prep: basis does not reproduce G.'); end
end
N = 1 + nj + p;
A = zeros(N); B = zeros(N); bv = zeros(N,1);
A(1,1) = 0.5;
if K > 0
    iz = (2 + nj):N;
    bv(iz) = Pm*dval(D,0);
    v = -0.5*(Pm*dval(D,1));
    A(iz,1) = v; A(1,iz) = v';
    % the one-step regression: z' Pm Phi_n z from m' G^{-1} int d' J; detrending first subtracts m' G^{-1} int d' Jtilde,
    % Jtilde = J - d' G^{-1} m, which leaves z' (Pm Phi_n)(Pm Phi_d) z
    if strcmp(D.detrend,'first'), M = (Pm*Phin)*(Pm*Phid); else, M = Pm*Phin; end
    A(iz,iz) = 0.5*(M + M');
    for i = 1:nj
        v = 0.5*(Pm*Dl(i,:)');
        A(iz,1+i) = v; A(1+i,iz) = v';
    end
    B(iz,iz) = -(Pm*Phid);
    B = (B + B')/2;
end
[L,Dv] = rangeeig(A);
P = mkprep(D,K,p,nj,N,deg,br,g,jl,A,B,L,Dv,G,bv,off,NaN,D.detrend);
end

function P = prep_gls(D,cbar)
K = size(D.C{1},2); np = numel(D.br) - 1; deg = size(D.C{1},1) - 1; br = D.br; off = offs(D);
[Z,~,~] = deltabasis(D);
k = size(Z,2);
if k == 0
    P = mkprep(D,0,0,0,1,deg,br,repmat({zeros(deg + 1,0)},1,np),zeros(1,0),0.5,0,1,0.5,zeros(0),0,off,cbar,'onestep');
    return
end
Cd = cell(1,np); Cd1 = Cd; Cd2 = Cd; Cps = Cd;
for j = 1:np, Cd{j} = D.C{j}*Z; Cd1{j} = dcoef(Cd{j}); Cd2{j} = dcoef(Cd1{j}); Cps{j} = Cd1{j} - cbar*Cd{j}; end
ev = @(Cs,j,r) dfdet_peval(Cs{j},r - off(j)).';
Psi = gram(Cps,br,D.loc);
Gd = gram(Cd,br,D.loc);
psi1 = ev(Cd1,np,1) - cbar*ev(Cd,np,1); delta1 = ev(Cd,np,1);
% the jumps of delta' (slope breaks); delta itself is continuous
jl = zeros(1,0); vs = zeros(k,0);
for j = 1:np-1
    lam = br(j+1);
    dd = ev(Cd,j + 1,lam) - ev(Cd,j,lam);
    if ~(max(abs(dd)) < 1e-9*max(1,max(abs(ev(Cd,j,lam))))), error('dfdet:prep','dfdet_prep: delta not continuous at %g.',lam); end
    s1 = ev(Cd1,j + 1,lam) - ev(Cd1,j,lam);
    if max(abs(s1)) > 1e-12, jl(end+1) = lam; vs(:,end+1) = s1; end %#ok<AGROW>
end
nj = numel(jl);
% orthonormal basis of span{delta, delta''}
F = cell(1,np); for j = 1:np, F{j} = [Cd{j}, Cd2{j}]; end
Gf = gram(F,br,D.loc);
[U,lv] = symeig(Gf);
keep = lv > 1e-11*max(lv);
U = U(:,keep); lv = lv(keep); p = numel(lv);
g = cell(1,np); for j = 1:np, g{j} = F{j}*U*diag(1./sqrt(lv)); end
Phi = U*diag(sqrt(lv)); Phid = Phi(1:k,:); Phi2 = Phi(k+1:2*k,:);
if ~(max(max(abs(Phid*Phid' - Gd))) < 1e-10*max(abs(Gd(:)))), error('dfdet:prep','dfdet_prep: basis does not reproduce G_d.'); end
N = 1 + nj + p;
H = zeros(k,N); H(:,1) = psi1;
for i = 1:nj, H(:,1+i) = -vs(:,i); end
iz = (2 + nj):N;
H(:,iz) = -(Phi2 - cbar^2*Phid);
PiH = Psi\H;
av = zeros(N,1); av(1) = 1; av = av - H'*(Psi\delta1);
A = av*av'/2;
S = zeros(p,N); S(:,iz) = eye(p);
B = -(H'*(Psi\(Phid*S)) + (S'*Phid')*PiH) + PiH'*Gd*PiH;
B = (B + B')/2;
[L,Dv] = rangeeig(A);
P = mkprep(D,k,p,nj,N,size(F{1},1) - 1,br,g,jl,A,B,L,Dv,Gd,zeros(N,1),off,cbar,'onestep');
end

function P = prep_pt(D,cbar)
% the limit of the point-optimal statistic against cbar, {S(abar) - abar S(1)}/s^2: Q = cbar^2 int J^2 - cbar J(1)^2
% - eta' Xi^{-1} eta + eta_0' Xi_0^{-1} eta_0 = cbar^2 V_g + ell' A ell, V_g the denominator of the GLS-detrended statistic
% and A = C - cbar^2 B_g (the construction of the GLS preparation with another numerator)
np = numel(D.br) - 1; deg = size(D.C{1},1) - 1; br = D.br; off = offs(D);
[Z,~,~] = deltabasis(D); k = size(Z,2);
if k == 0                                               % only the path at one: Q = cbar^2 int J^2 - cbar J(1)^2
    P = mkprep(D,0,0,0,1,deg,br,repmat({zeros(deg + 1,0)},1,np),zeros(1,0),-cbar,0,1,-cbar,zeros(0),0,off,cbar,'onestep');
    return
end
Cd = cell(1,np); Cd1 = Cd; Cd2 = Cd; Cps = Cd;
for j = 1:np, Cd{j} = D.C{j}*Z; Cd1{j} = dcoef(Cd{j}); Cd2{j} = dcoef(Cd1{j}); Cps{j} = Cd1{j} - cbar*Cd{j}; end
ev = @(Cs,j,r) dfdet_peval(Cs{j},r - off(j)).';
Xi = gram(Cps,br,D.loc); Xi0 = gram(Cd1,br,D.loc); Gd = gram(Cd,br,D.loc);
psi1 = ev(Cd1,np,1) - cbar*ev(Cd,np,1); dot1 = ev(Cd1,np,1);
jl = zeros(1,0); vs = zeros(k,0);
for j = 1:np-1
    lam = br(j+1); s1 = ev(Cd1,j + 1,lam) - ev(Cd1,j,lam);
    if max(abs(s1)) > 1e-12, jl(end+1) = lam; vs(:,end+1) = s1; end %#ok<AGROW>
end
nj = numel(jl);
F = cell(1,np); for j = 1:np, F{j} = [Cd{j}, Cd2{j}]; end
[U,lv] = symeig(gram(F,br,D.loc));
keep = lv > 1e-11*max(lv); U = U(:,keep); lv = lv(keep); p = numel(lv);
g = cell(1,np); for j = 1:np, g{j} = F{j}*U*diag(1./sqrt(lv)); end
Phi = U*diag(sqrt(lv)); Phid = Phi(1:k,:); Phi2 = Phi(k+1:2*k,:);
N = 1 + nj + p; iz = (2 + nj):N;
H = zeros(k,N); H0 = zeros(k,N); H(:,1) = psi1; H0(:,1) = dot1;
for i = 1:nj, H(:,1+i) = -vs(:,i); H0(:,1+i) = -vs(:,i); end
H(:,iz) = -(Phi2 - cbar^2*Phid); H0(:,iz) = -Phi2;
C = -(H'*(Xi\H)) + H0'*(Xi0\H0); C(1,1) = C(1,1) - cbar;
Pg = prep_gls(D,cbar); C = C - cbar^2*Pg.B; C = (C + C')/2;
[L,Dv] = rangeeig(C);
P = mkprep(D,k,p,nj,N,size(F{1},1) - 1,br,g,jl,C,Pg.B,L,Dv,Gd,zeros(N,1),off,cbar,'onestep');
end

function [P,a1] = prep_lm(D)
% the denominator V = int Sbar^2 of the LM statistic (Schmidt and Phillips; with breaks, Amsler and Lee, Lee and
% Strazicich): the trend coefficients estimated from the first differences, S_t = y_t - psi - Z_t delta; with delta the
% trend functions that vanish at the origin and do not jump, Xi_0 = int delta' delta'' and eta_0 = int delta' dJ, the
% detrended path is J_0 = J - delta' Xi_0^{-1} eta_0, the second regression projects it off the derivatives,
% Sbar = J_0 - delta'' Xi_0^{-1} int delta' J_0, and the limit is {J_0(1)^2 - 1}/2/sqrt(V) for every c
np = numel(D.br) - 1; br = D.br; off = offs(D);
[Z,~,~] = deltabasis(D); k = size(Z,2);
if k == 0, error('dfdet:args','prep_lm: no trend function that vanishes at the origin and does not jump.'); end
Cd = cell(1,np); Cd1 = Cd; Cd2 = Cd;
for j = 1:np, Cd{j} = D.C{j}*Z; Cd1{j} = dcoef(Cd{j}); Cd2{j} = dcoef(Cd1{j}); end
ev = @(Cs,j,r) dfdet_peval(Cs{j},r - off(j)).';
jl = zeros(1,0); vs = zeros(k,0);
for j = 1:np-1
    lam = br(j+1); s1 = ev(Cd1,j + 1,lam) - ev(Cd1,j,lam);
    if max(abs(s1)) > 1e-12, jl(end+1) = lam; vs(:,end+1) = s1; end %#ok<AGROW>
end
nj = numel(jl);
F = cell(1,np); for j = 1:np, F{j} = [Cd{j}, Cd1{j}, Cd2{j}]; end
[U,lv] = symeig(gram(F,br,D.loc));
keep = lv > 1e-11*max(lv); U = U(:,keep); lv = lv(keep); p = numel(lv);
g = cell(1,np); for j = 1:np, g{j} = F{j}*U*diag(1./sqrt(lv)); end
Phi = U*diag(sqrt(lv)); Pd = Phi(1:k,:); P1 = Phi(k+1:2*k,:); P2 = Phi(2*k+1:3*k,:);
N = 1 + nj + p; iz = (2 + nj):N;
H0 = zeros(k,N); H0(:,1) = ev(Cd1,np,1);
for i = 1:nj, H0(:,1+i) = -vs(:,i); end
H0(:,iz) = -P2;
Xi0 = P1*P1'; S = zeros(p,N); S(:,iz) = eye(p);
Bb = Xi0\H0;                                            % the coefficients of delta: b = Bb ell
Ba = Xi0\(P1*S - (P1*Pd')*Bb);                          % those of delta': a = Ba ell
W = [Bb; Ba]; T = [Pd; P1];                             % the fitted function is (delta, delta')' W ell = g' T' W ell
B = -(S'*T'*W + W'*T*S) + W'*(T*T')*W; B = (B + B')/2;
a1 = zeros(N,1); a1(1) = 1; a1 = a1 - Bb'*ev(Cd,np,1);   % J_0(1) = a1' ell
A = zeros(N); A(1,1) = 0.5;                             % placeholder
% (detrend 'lm': V is not the norm of an orthogonal projection of J, and DFDET_TSTAR takes the estimate of the edge t_* from
% the quadratic form itself)
P = mkprep(D,k,p,nj,N,size(F{1},1) - 1,br,g,jl,A,B,[1; zeros(N - 1,1)],0.5,Pd*Pd',zeros(N,1),off,NaN,'lm');
end

% ---------------------------------------------------------------------------------------------------------------------
function P = mkprep(D,K,p,nj,N,deg,br,g,jl,A,B,L,Dv,G,b,off,cbar,detrend)
Dv = Dv(:);
P = struct('D',D,'K',K,'p',p,'nj',nj,'N',N,'deg',deg,'br',br,'g',{g},'jl',jl(:).','A',A,'B',B,'L',L,'Dv',Dv, ...
    'r',numel(Dv),'rp',nnz(Dv > 0),'rm',nnz(Dv < 0),'G',G,'b',b(:),'off',off,'cbar',cbar,'detrend',detrend);
end

function off = offs(D)
% the origins of the variable of the coefficients of the pieces: 0 (global r) or br(j) (local)
np = numel(D.br) - 1;
if D.loc, off = D.br(1:np); else, off = zeros(1,np); end
end

function v = dval(D,r)
% the regressors at r (a column)
np = numel(D.br) - 1; j = 1;
while j < np && r > D.br(j+1), j = j + 1; end
off = 0; if D.loc, off = D.br(j); end
v = dfdet_peval(D.C{j},r - off).';
end

function [jl,Dl] = jumps(D)
% the break points at which the regressors jump, and the jumps (rows)
np = numel(D.br) - 1; K = size(D.C{1},2); off = offs(D);
jl = zeros(1,0); Dl = zeros(0,K);
for j = 1:np-1
    lam = D.br(j+1);
    dl = dfdet_peval(D.C{j+1},lam - off(j+1)) - dfdet_peval(D.C{j},lam - off(j));
    if ~isempty(dl) && max(abs(dl)) > 1e-12, jl(end+1) = lam; Dl(end+1,:) = dl; end %#ok<AGROW>
end
end

function [Z,Mr,k] = deltabasis(D)
% a basis Z of the directions b with d(0)' b = 0 and Delta_i' b = 0 (the trend functions that vanish at the origin and do
% not jump), by the singular value decomposition with the absolute tolerance 1e-10 max(1, max |Mr|)
K = size(D.C{1},2);
if K == 0, Z = zeros(0,0); Mr = zeros(1,0); k = 0; return, end
[~,Dl] = jumps(D);
Mr = [dval(D,0).'; Dl];
tol = 1e-10*max(1,max(abs(Mr(:))));
[~,~,V] = svd(Mr); sv = svd(Mr);
rk = nnz(sv > tol);
Z = V(:,rk+1:end); k = size(Z,2);
end

function Dm = dcoef(C)
% the derivative of a coefficient matrix (rows: powers), same size, last row zero
Dm = zeros(size(C));
n = size(C,1);
if n > 1, Dm(1:n-1,:) = (1:n-1)'.*C(2:n,:); end
end

function G = gram(F,br,loc)
% the Gram matrix of the columns of the piecewise polynomials F{j} on [0, 1]
n = size(F{1},2); deg = size(F{1},1) - 1;
G = zeros(n);
if loc
    % local coefficients: Gauss-Legendre with deg + 2 nodes per piece (exact for the products), on the values (the sum of
    % the products of the monomial coefficients would cancel quadratically at high degree)
    [z,w] = glrule(deg + 2);
    for j = 1:numel(F)
        L = br(j+1) - br(j);
        V = dfdet_peval(F{j},L*z);
        G = G + V'*(V.*(L*w));
    end
    G = (G + G')/2;
    return
end
kl = (0:deg)' + (0:deg) + 1;
for j = 1:numel(F)
    M = (br(j+1).^kl - br(j).^kl)./kl;
    G = G + F{j}'*M*F{j};
end
G = (G + G')/2;
end

function [U,lam] = symeig(G)
% the eigen-decomposition of a symmetric matrix, ascending
[U,E] = eig((G + G')/2);
[lam,ix] = sort(diag(E)); U = U(:,ix);
end

function [L,Dv] = rangeeig(A)
% the nonzero eigenvalues of the symmetric A (|.| > 1e-10 times the largest) and their eigenvectors
[U,ev] = symeig(A);
keep = abs(ev) > 1e-10*max(abs(ev));
Dv = ev(keep); L = U(:,keep);
end
