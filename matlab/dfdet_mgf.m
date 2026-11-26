function [theta,ld,mu,fx] = dfdet_mgf(Mo,c,t,light)
%DFDET_MGF  The joint moment generating function of the general engine at one t: log D_V(t), the roots, the fixed start.
%
%   [theta,ld,mu,fx] = DFDET_MGF(Mo,c,t) for a model Mo (DFDET_MODEL) at the
%   local alternative c and one complex t: theta = sqrt(c^2 - 2t) (principal),
%   ld = log D_V(t) = c + theta + log det M0 with M0 = I - 2 (t B + gamma E11)
%   Sigma(theta), gamma = (c + theta)/2 (D_V = M_V(t)^-2, an entire function of
%   t; the imaginary part is the principal one, not continued), and the r
%   eigenvalues mu_k of H = 2 Dv L' Sigma M0^{-1} L, so that
%   det(I - 2 (s A + t B + gamma E11) Sigma) = det M0 prod_k (1 - s mu_k): the
%   roots in the argument of U are s_k = 1/mu_k.  The (1,1) entry of M0 is
%   computed without the cancellation of theta - c for c > 0.  With a
%   stationary start the extended matrices (DFDET_MODEL: the covariance
%   [Sigma + v h h', v h; v h', v], v = -1/(2c), and gamma (E11 - ENN)).  With a
%   fixed start fx holds the data of the start at t (the joint moment
%   generating function with a fixed start: log Mt = log Mt_0 + g0^2 (cf +
%   E1(s)), E1(s) = s [p1' (I - sH)^{-1} w1 + p2' (I - sH)^{-1} w2], by
%   Woodbury's identity from the Gaussian integral over zeta ~ N(x mt, Sigma)):
%   the fields cf, H, w1, p1, w2, p2 and mt (the mean of L'(ell + g0 yb) under
%   the tilt per unit g0); fx = [] otherwise.  DFDET_MGF(Mo,c,t,true) returns
%   theta and ld only (and mu, fx empty), the zero-start parts of a fixed start.
%
%   See also DFDET_MODEL, DFDET_SIGMA, DFDET_DVAL.
if nargin < 4, light = false; end
if isfield(Mo,'vm'), [theta,ld,mu] = vol_parts(Mo,t,light); fx = []; return, end
P = Mo.P;
theta = sqrt(complex(c*c - 2*t)); gam = (c + theta)/2;
fx = []; mu = [];
if strcmp(Mo.start,'stationary')
    [S,h] = dfdet_sigma(P,theta);
    v = -1/(2*c); n = P.N; N = n + 1;
    Sx = complex(zeros(N)); Sx(1:n,1:n) = S + v*(h*h.'); Sx(N,1:n) = v*h.'; Sx(1:n,N) = v*h; Sx(N,N) = v;
    Q = t*Mo.B; Q(1,1) = Q(1,1) + gam; Q(N,N) = Q(N,N) - gam;
    M0 = eye(N) - 2*(Q*Sx);
    ld = c + theta + ldet(M0);
    if light, return, end
    X = M0\Mo.L;
    H = 2*(Mo.Dv.*(Mo.L.'*(Sx*X)));
    mu = eig(H);
    return
end
if strcmp(Mo.start,'fixed') && ~light
    [S,h] = dfdet_sigma(P,theta);
else
    S = dfdet_sigma(P,theta);
end
M0 = buildM(P,S,c,0,t,theta);
ld = c + theta + ldet(M0);
if light, return, end
L = P.L; Dv = P.Dv;
X = M0\L;
H = 2*(Dv.*(L.'*(S*X)));
mu = eig(H);
if strcmp(Mo.start,'fixed')
    yb = Mo.yb;
    Q0 = t*P.B; Q0(1,1) = Q0(1,1) + gam;
    mt = h + yb;
    v0 = Q0*(h - yb); qt = -2*(Q0*yb);
    a0 = M0\v0; aq = M0\qt;
    Sq = S*qt;
    cf = -gam + yb.'*(Q0*yb) + mt.'*a0 + qt.'*(S*aq)/2;
    Mt = M0.';
    p1 = L.'*(Mt\mt);
    w1 = 2*Dv.*(L.'*(S*a0)) + Dv.*(L.'*mt);
    p2 = L.'*(Mt\Sq);
    w2 = Dv.*(L.'*(S*aq));
    mm = L.'*((Mt\h) + yb);
    fx = struct('cf',cf,'H',H,'w1',w1,'p1',p1,'w2',w2,'p2',p2,'mt',mm);
end
end

function M = buildM(P,S,c,s,t,theta)
% I - 2 (s A + t B + gamma E11) Sigma, its (1,1) entry without cancellation: m011 - 2 ((s A + t B) Sigma)_11, formed from
% these two terms alone (above the unit root it is of the order e^{-2c} where the transform matters, and the earlier
% form 1 - 2 ((s A + t B) Sigma)_11 - 2 gamma Sigma_11 + (m011 - (1 - 2 gamma Sigma_11)) carried an absolute rounding
% error of 1e-16, which put the moments off by 6e-4 relative at c = 15 and by 0.95 at c = 20)
N = P.N; gam = (c + theta)/2;
Q = -2*((s*P.A + t*P.B)*S);
M = eye(N) + Q;
M(1,:) = M(1,:) - 2*gam*S(1,:);
M(1,1) = m011(c,t,theta) + Q(1,1);
end

function v = m011(c,t,theta)
% 1 - (c + theta) e_0(2 theta) = 1 - 2 gamma Var J(1), without the cancellation of theta - c for c > 0
if abs(theta) < 0.5
    z = 2*theta;
    if abs(z) < 1e-3, e = 1 - z/2 + z*z/6 - z^3/24; else, e = (1 - exp(-z))/z; end
    v = 1 - (c + theta)*e; return
end
if c > 0, d = -2*t/(theta + c); else, d = theta - c; end
v = (d + (c + theta)*exp(-2*theta))/(2*theta);
end

function l = ldet(M)
% log det M from the LU factors: log|det| plus the principal logarithm of the phase (LinearAlgebra's logdet of an LU)
[~,U,p] = lu(M,'vector');
u = diag(U); a = abs(u);
sg = prod(u./a);
% the parity of the permutation p: n minus the number of its cycles
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

function [thl,ld,mu] = vol_parts(VM,t,light)
% the volatility transform at t (mgfparts of a VolModel): theta_j = sqrt(c_j^2 - 2 t omega_j^2), the tilt Gamma (g_last at
% J(1), g_j - g_{j+1} at the boundary between pieces j and j + 1, g_j = (c_j + theta_j)/(2 omega_j^2)), log D_V(t) =
% sum_j (c_j + theta_j) L_j + log det{I - 2 (t B + Gamma) Sigma} and the reciprocal roots mu, the eigenvalues of
% 2 Dv L' Sigma M0^{-1} L
[th,gam,ex] = dfdet_vol('tilt',VM,t);
S = dfdet_sigma(VM,th);
M0 = eye(VM.Ne) - 2*((t*VM.B + diag(gam))*S);
ld = ex + ldet(M0); thl = th(end); mu = [];
if light, return, end
X = M0\VM.Lx;
H = 2*(VM.Dv.*(VM.Lx.'*(S*X)));
mu = eig(H);
end
