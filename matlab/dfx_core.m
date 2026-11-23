function varargout = dfx_core(action, varargin)
%DFX_CORE  Engine of the DFX routines: a MATLAB port of DFExact.jl (Bao, Broda and Kan,
%   "The exact limiting distribution of the Dickey-Fuller t-ratio and some related
%   statistics", repository df_limiting).
%
%   Not called directly; use DFX_CDF, DFX_PDF, DFX_QUANTILE, DFX_MOMENTS and DFX_TSTAR.
%
%   The limit tau = U/sqrt(V) of the Dickey-Fuller t-ratio, models 1, 2, 3 (no deterministic
%   terms, intercept, intercept and trend; cases DK, DU, DT, and SK, SU, ST with a stationary
%   start), under phi = 1 + c/n, zero drift.  The joint mgf of (U,V) is e^(-s/2) D(s,t)^(-1/2)
%   with D a polynomial of degree one (model 1) or two in s; the inversion in s is explicit (a
%   scaled chi-square(1) density, or K_0 for a difference of two), and the density and the
%   distribution function are single integrals along Re t = tau (the line integral representation), with
%   inner Gaussian transforms of the s-inverse evaluated by Gauss-Legendre rules on complex
%   paths.  The port follows the Julia code function by function; the inner integrands are
%   vectorized over their nodes, and the loop over the nodes of the outer integral is serial
%   (the Julia code runs it on several threads).
%
%   Actions:
%     dfx_core('cdf', q, model, c, opts), dfx_core('pdf', q, model, c, opts)
%     dfx_core('tstar', model, c, stat), dfx_core('moments', model, c, opts)
%     dfx_core('nodefac', k)      multiplies every node count by k (k = 2: convergence check)
%     dfx_core('k0x', z)          e^z K_0(z), Re z >= 0 (for testing)
%     dfx_core('D', model, c, s, t, stat)   the determinant D(s,t)
%     dfx_core('count')           [kernel evaluations, inner calls] since the last call (testing)
%     dfx_core('tune', [NGROW NLEV NCON])   node counts of the segment pieces (testing)
%     [P, s1, s2, M0] = dfx_core('tparts', model, c, t, stat)   P(t), the roots s1(t), s2(t) of D(., t)
%                                 (s1 the one with Re < 0 on the line) and M0(t) at the points t, the branches
%                                 fixed at each point (after Bao and Kan); dfx_core('gstat', model, z), the entire function g of the
%                                 initial-value terms; dfx_core('tau', model, c, stat), the default line
%                                 (the pieces of the engine used by TR_CONTOUR)
%   opts: start (zero or stationary), c0 (0; the local null value: the t-ratio for
%         phi0 = 1 + c0/n, with the limit (U - c0 V)/sqrt(V); the kernel then runs along the
%         parabola u = q v + c0 v^2, see BREAKPTS, INNER_GEN, SEGINT, OUTER_GEN), gen (false;
%         true: the general kernel also at c0 = 0, for checking), tau (line Re t = tau),
%         hu (0.075), N (20), pfac (1; subdivides the Filon panels), glfac (1; multiplies the
%         node counts of the Gauss-Legendre route of the smooth part), wmax (4900), edge_hw
%         (trapezoidal rule for the edge part, for checking), uniform (0; hw > 0: the
%         trapezoidal rule in w with step hw for the whole integrand, for checking).
switch lower(action)
    case {'cdf','pdf'}
        q = varargin{1}; model = varargin{2};
        if ~(isnumeric(model) && isscalar(model) && any(model == [1 2 3 4]))
            error('dfx: model must be 1, 2, 3 or 4 (no deterministic terms, intercept, intercept and trend, GLS-detrended trend)');
        end
        c = 0; if numel(varargin) >= 3 && ~isempty(varargin{3}), c = varargin{3}; end
        o = struct(); if numel(varargin) >= 4 && ~isempty(varargin{4}), o = varargin{4}; end
        out = zeros(size(q));
        for k = 1:numel(q)
            out(k) = onevalue(lower(action), q(k), model, c, o);
        end
        varargout{1} = out;
    case 'tstar'
        stat = false; if numel(varargin) >= 3, stat = varargin{3}; end
        varargout{1} = tstar(varargin{1}, varargin{2}, stat);
    case 'moments'
        model = varargin{1}; c = 0; if numel(varargin) >= 2 && ~isempty(varargin{2}), c = varargin{2}; end
        o = struct(); if numel(varargin) >= 3 && ~isempty(varargin{3}), o = varargin{3}; end
        varargout{1} = moments(model, c, o);
    case 'nodefac'
        nodefac(varargin{1});
    case 'k0x'
        varargout{1} = k0x(varargin{1});
    case 'count'
        varargout{1} = k0x([]); varargout{2} = inner_gen();
    case 'tune'
        varargout{1} = tune(varargin{:});
    case 'd'
        varargout{1} = Dfun(varargin{:});
    case 'tparts'
        stat = false; if numel(varargin) >= 4, stat = varargin{4}; end
        [varargout{1:max(nargout,1)}] = tparts(varargin{1}, varargin{2}, varargin{3}, stat);
    case 'gstat'
        varargout{1} = gstat(varargin{1}, varargin{2});
    case 'tau'
        stat = false; if numel(varargin) >= 3, stat = varargin{3}; end
        varargout{1} = default_tau(varargin{1}, varargin{2}, stat);
    case 'dbg'
        [varargout{1:nargout}] = feval(varargin{:});
    otherwise
        error('dfx_core: unknown action %s', action);
end
end

% =====================================================================================
function v = onevalue(kind, q, model, c, o)
stat = startflag(getf(o, 'start', 'zero'), c);
tau = getf(o, 'tau', []);
if isempty(tau), tau = default_tau(model, c, stat); end
kern = kernel(model, stat);
c0 = getf(o, 'c0', 0);
if q == 0 && c0 == 0 && getf(o, 'gen', false), q = 1e-12; end   % the general kernel forced at x = 0 with c0 = 0; the standard route and the c0 ~= 0 kernel take x = 0 exactly
if strcmp(kind, 'pdf'), weight = @(t) 1; kd = 1; else, weight = @(t) 1 / t; kd = 2; end
hw = getf(o, 'uniform', 0);
if c0 ~= 0 || getf(o, 'gen', false)
    % the local null value c0: the kernel runs along the parabola u = q v + c0 v^2 (BREAKPTS,
    % INNER_GEN), with one Filon part of the outer integral per breakpoint
    R = breakpts(q, c0, tau);
    if R.double && kern == 1 && kd == 1, v = Inf; return; end   % case DK: the density is infinite at x = -sqrt(2 c0)
    ed = [0 R.b]; mid = (ed(1:end-1) + ed(2:end)) / 2;
    R.negl = 3 * tau * mid.^2 / 4 > 40;                % the edge part is below e^-40 there
    R.vsad = saddle_loc(model, c, q, c0, tau, stat, kern, R);
    inner = @(t, s1, s2, M0, part) inner_gen(kern, t, R, s1, s2, M0, kd, part);
    if hw > 0
        v = uniform(model, c, tau, stat, inner, weight, hw);
    else
        v = outer_gen(model, c, q, c0, R, tau, inner, weight, o, stat);
    end
    return
end
e1 = k0edge1(model, c, stat, tau);
if kd == 1
    inner = @(t, s1, s2, M0, part) inner_density(kern, t, q, s1, s2, M0, part, e1);
else
    inner = @(t, s1, s2, M0, part) inner_cdf(kern, t, q, s1, s2, M0, part, e1);
end
if hw > 0
    v = uniform(model, c, tau, stat, inner, weight, hw);
else
    v = outer(model, c, q, tau, inner, weight, o, stat);
end
end

function x = getf(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), x = o.(f); else, x = d; end
end

function k = nodefac(set)
persistent NF
if isempty(NF), NF = 1; end
if nargin, NF = set; end
k = NF;
end

function s = startflag(start, c)
switch lower(start)
    case 'zero', s = false;
    case 'stationary'
        if ~(c <= -1e-4), error('ltratio:args','ltratio: a stationary start (cases SK, SU, ST) requires c <= -1e-4: closer to zero the variance -1/(2c) of the start is too large for the quadrature.'); end   % (LTRATIO_ARGS)
        s = true;
    otherwise, error('dfx: start must be ''zero'' or ''stationary''');
end
end

function k = kernel(model, stat)
% rank of the quadratic form R under the tilt (model 4, the GLS-detrended trend case: one)
if (model == 1 || model == 4) && ~stat, k = 1; else, k = 2; end
end

function S = series()
persistent SS
if isempty(SS), SS = dfx_series; end
S = SS;
end

% ---------- entire functions of z = mu^2 ----------
function [C, Q, S] = pieces(z)
% cosh mu, (1 - cosh mu)/z, sinh mu / mu, each entire in z, with series near z = 0
mu = sqrt(z);
C = cosh(mu);
S = sinh(mu) ./ mu;
Q = (1 - C) ./ z;
m = abs(z) < 1e-3;
if any(m(:))
    zm = z(m);
    S(m) = 1 + zm / 6 + zm .* zm / 120 + zm.^3 / 5040;
    Q(m) = -1/2 - zm / 24 - zm .* zm / 720 - zm.^3 / 40320;
end
end

% ---------- the stationary start ----------
function g = gstat(model, z)
S = series();
if model == 2, ser = S.G2; else, ser = S.G3; end
g = complex(zeros(size(z)));
m = abs(z) < 8;
if any(m(:))
    zm = z(m); acc = zeros(size(zm));
    for j = numel(ser):-1:1, acc = acc .* zm + ser(j); end
    g(m) = acc;
end
if any(~m(:))
    zz = z(~m);
    [C, Q, Sv] = pieces(zz);
    if model == 2
        g(~m) = (2 * Q + Sv) ./ zz;
    else
        g(~m) = (-12 * (1 - zz .* Q) + Q .* (24 - 4 * zz) + Sv .* (24 + zz)) ./ zz.^3;
    end
end
end

function a = Dstat_add(model, c, s, t)
if model == 2, ck = c; else, ck = c^3; end
a = ck * exp(c) * (t + s .* s / 2 + c * s) .* gstat(model, c * c - 2 * t);
end

function D = Dfun(model, c, s, t, stat)
% D(s,t); log M = -s/2 + extra - log(D)/2, extra = -c/2 for model 1 with a zero start
if nargin < 5, stat = false; end
if isscalar(s) && ~isscalar(t), s = s * ones(size(t)); end
if isscalar(t) && ~isscalar(s), t = t * ones(size(s)); end
z = c * c - 2 * t;
e = exp(c);
if stat
    if model == 1
        [C, ~, S] = pieces(z);
        D = e * (C - c * S + S .* (2 * t + s .* s) / (2 * c));
    else
        D = Dfun(model, c, s, t, false) + Dstat_add(model, c, s, t);
    end
    return
end
if model == 1
    [C, ~, S] = pieces(z);
    D = C - (c + s) .* S;
    return
end
lam = -s; lx = t;
D = complex(zeros(size(z)));
m = abs(z) < 8;                     % Taylor series in z (cancellation in num/z^k)
if any(m(:))
    SS = series();
    if model == 2, T = SS.D1; else, T = SS.D2; end
    nj = size(T, 1); ni = size(T, 2); nk = size(T, 3);
    cp = c .^ (0:nk-1);
    A = zeros(nj, ni);              % A(j,i): coefficient of z^(j-1) lam^(i-1)
    for i = 1:ni, A(:, i) = reshape(T(:, i, :), nj, nk) * cp(:); end
    lm = lam(m); zm = z(m); acc = zeros(size(zm));
    for j = nj:-1:1
        cj = A(j, ni) * ones(size(lm));
        for i = ni-1:-1:1, cj = cj .* lm + A(j, i); end
        acc = acc .* zm + cj;
    end
    D(m) = e * acc;
end
if any(~m(:))
    zz = z(~m); l = lam(~m); x = lx(~m);
    [C, Q, S] = pieces(zz);
    if model == 2
        A = 2 * c * c * l - 4 * c * x - 2 * l .* l;
        B = -c^3 + c * c * l - 2 * x - l .* l;
        D(~m) = e * (c * c * (1 - zz .* Q) + A .* Q + B .* S) ./ zz;
    else
        c2 = c^2; c3 = c^3; c4 = c^4; c5 = c^5; c6 = c^6; c7 = c^7;
        P1 = c6 - 12 * c4 * l + 2 * x * (-c4 + 12 * c3 - 12 * c2) + 48 * l .* l * (c2 - 3 * c + 3);
        P2 = (4 * l * (-c6 + 6 * c4) + 16 * l .* l * (c4 - 3 * c3 - 3 * c2 + 18 * c - 18) ...
            + 8 * x .* (c5 - 3 * c4 + c4 * l - 4 * l .* l * (c2 - 3 * c + 3)) + 16 * x .* x * (-c3 + 3 * c2 - 6 * c + 6));
        P3 = (-c7 + l * (c6 + 24 * c4) + 4 * l .* l * (-c4 + 3 * c3 - 27 * c2 + 72 * c - 72) ...
            + 2 * x .* (c5 - 4 * c4 - 12 * c3 + 12 * c2 - c4 * l + 4 * l .* l * (c2 - 3 * c + 3)) + 16 * x .* x * (c2 - 3 * c + 3));
        D(~m) = e * ((1 - zz .* Q) .* P1 + Q .* P2 + S .* P3) ./ zz.^3;
    end
end
end

function ph = unwrapc(ph)
% the unwrapping of DFExact.unwrap!: ph(j) -= 2 pi round((ph(j) - ph(j-1))/(2 pi)), sequentially
if numel(ph) < 2, return; end
k = [0; cumsum(round(diff(ph(:)) / (2 * pi)))];
ph = ph(:) - 2 * pi * k;
end

% ---------- P(t), the roots and M0 at each point: the branch after Bao and Kan ----------
% The branch of P(t) = (-a_2(t))^(-1/2), of M_0(t) = (-s_1 s_2)^(-1/2) and the labels of the roots s_1(t), s_2(t) of D(., t)
% were fixed by continuity along a grid from the real point until November 2026 (an unwrapping of the phases).  They are
% determined at each point, after Yong Bao and Raymond Kan ("Exact Finite-Sample Distribution and Moments of the Least
% Squares Estimator in First-Order Autoregressions: O(1) Algorithms"), who write each argument under a root as a product
% of factors whose arguments are known (tparts_point of DFExact.jl, dfl_tparts of the C engine).  For Im t > 0 (the
% conjugate for Im t < 0), with mu = (c^2 - 2t)^(1/2) the principal root: the reciprocal roots are the eigenvalues of
% K Sigma_t, y*Ky / y*Sigma_t^(-1)y, a real number over a number with a positive real part for Re t < t_* (the lemma on
% accretive matrices of Bao and Kan) and a negative imaginary part for Im t > 0, so that s_1 (negative type) lies in the
% open left half-plane on the line and in the open upper half-plane, s_2 in the right and the lower one: the labels are
% read off the signs; M_0 = (-s_1)^(-1/2) s_2^(-1/2) with principal roots on the line; D(0, t) = a_0 = (-a_2)(-s_1) s_2
% has the continuous argument of a canonical product over the positive zeros of the moment generating function of V:
% without deterministic terms e^c cosh(mu) G, G = 1 - c tanh(mu)/mu, with the argument Im mu + Arg(1 + e^(-2 mu)) of
% cosh(mu) and G in a half-plane; a regressor projected out, or a start integrated out, changes V by a quadratic form of
% rank one of one sign, and the principal argument of the ratio of the two determinants, in (0, pi) for a projection
% and in (-pi, 0) for the stationary start, is the continuous one (two projections: (0, 2 pi); GLS detrending, one update
% of each sign: the principal argument).  Then arg(-a_2) = arg D(0, t) - Arg(-s_1) - Arg(s_2), all at the point.
function a = argpos(a, two)
% an argument known to lie in (0, pi) (one projection) or (0, 2 pi) (two) from its principal value
if two, m = a < -1e-6; else, m = a < -pi / 2; end
a(m) = a(m) + 2 * pi;
end
function a = argneg(a), m = a > pi / 2; a(m) = a(m) - 2 * pi; end
function [aK, mu, e2, G] = argK(c, t)
% the continuous argument (zero at real t below t_*) of cosh(mu) - c sinh(mu)/mu = cosh(mu) G at Im t > 0
z = c * c - 2 * t; mu = sqrt(z); e2 = exp(-2 * mu);
if c > 0, del = -2 * t ./ (mu + c); else, del = mu - c; end
G = (del + (mu + c) .* e2) ./ (mu .* (1 + e2));
sm = abs(mu) < 1e-3; G(sm) = 1 - c * (1 - z(sm) / 3 + 2 * z(sm) .* z(sm) / 15);
aK = imag(mu) + angle(1 + e2) + angle(G);
end
function [s1, s2] = labelroots(ra, rb)
% the two roots labelled by their signs: s_1 in the left half-plane and s_2 in the right one; where both real parts have the
% same sign (beyond t_*), s_1 the root in the upper half-plane
d = (real(ra) < 0) ~= (real(rb) < 0);
first = (d & real(ra) < real(rb)) | (~d & imag(ra) > imag(rb));
s1 = rb; s2 = ra; s1(first) = ra(first); s2(first) = rb(first);
end

function [P, s1, s2, M0] = tparts(model, c, t, stat)
% P(t), s1(t), s2(t) and M0(t) = (-s1 s2)^(-1/2) at the points t, in any order: each point by itself, the branches of the
% square roots and the labels of the roots fixed at the point (after Bao and Kan, above) with no reference to its
% neighbours.  For c > 2 with a zero start the roots and P(t) come from the Girsanov forms (IIPAR), which are free of the
% cancellation of D.  (Against the continuation along the grid that this replaced, within 5.7e-14 relative on the
% engine's grids, November 2026, as in the Julia and the C engine.)
t = complex(t(:)); nt = numel(t);
neg = imag(t) < 0; t(neg) = conj(t(neg));
re = imag(t) == 0;
s2 = complex(zeros(nt, 1)); M0 = s2;
if model == 1 && ~stat
    z = c * c - 2 * t;
    mu = sqrt(z); e2 = exp(-2 * mu); om = 1 - e2;
    if c > 0
        % s1 = mu coth mu - c without the cancellation of C - c S at large c: mu - c = -2t/(mu + c), mu coth mu - mu =
        % 2 mu e^(-2mu)/(1 - e^(-2mu))
        ep = 2 * mu .* e2 ./ om;
        sm = abs(mu) < 1e-4; ep(sm) = 1 - mu(sm) + mu(sm).^2 / 3;
        s1 = -2 * t ./ (mu + c) + ep;
    else
        s1 = mu .* (1 + e2) ./ om - c;
    end
    % sinh(mu)/mu = e^mu (1 - e^(-2 mu))/(2 mu): its modulus and its continuous argument Im mu + Arg(1 - e^(-2 mu)) - Arg mu
    lS = real(mu) + log(abs(om)) - log(2 * abs(mu));
    ph = imag(mu) + angle(om) - angle(mu); ph(re) = 0;
    P = exp(-c / 2 - lS / 2) .* exp(-1i * ph / 2);
    sz = abs(z) < 1e-3;                                % (the series of pieces; S near 1, its principal argument)
    if any(sz)
        [C, ~, S] = pieces(z(sz));
        if c <= 0, s1(sz) = (C - c * S) ./ S; end
        P(sz) = exp(-c / 2) * abs(S).^(-1/2) .* exp(-1i * angle(S) / 2);
    end
elseif model == 4 && ~stat
    % the GLS-detrended trend case: M(s,t) = e^(-s/2) P(t) (s1 - s)^(-1/2) as in case DK, with M_Vg(t) = e^(-c/2)
    % (cosh th det(I - Q Sig))^(-1/2) and s1 from GLS_PARTS; the argument of D from that of 1 - c tanh(th)/th and the
    % principal argument of the ratio (one update of each sign)
    [lch, D, s1, GK] = gls_parts(t, c);
    ph = angle(GK) + angle(D ./ GK); ph(re) = 0;
    P = exp(-c / 2 - (lch + log(abs(D)) + 1i * ph) / 2 + (log(abs(s1)) + 1i * angle(s1)) / 2);
elseif ~stat && c > 2
    % cases DU, DT at large c: D(s,t)/e^c cancels catastrophically (terms of size e^(2c) against a value near 1): the
    % roots and the prefactor from the Girsanov forms (IIPAR): A0 s^2 - A1 s + A2 = 0 has the roots s_k = 1/(2 kap_k),
    % M(0,t) = e^(-c/2) (cosh th det(I - Q Sig))^(-1/2), and P = M(0,t) (-s1 s2)^(1/2)
    cas = 'DU'; if model == 3, cas = 'DT'; end
    Q = iipar(t, c, cas);
    disc = sqrt(Q.A1 .* Q.A1 - 4 * Q.A0 .* Q.A2);
    sw = sign(real(conj(Q.A1) .* disc)); sw(sw == 0) = 1;
    rb = (Q.A1 + sw .* disc) ./ (2 * Q.A0); rs = Q.A2 ./ (Q.A0 .* rb);
    [s1, s2] = labelroots(rb, rs);
    ph = angle(Q.N11) + argpos(angle(Q.D ./ Q.N11), model == 3); ph(re) = 0;
    lm = log(abs(s1 .* s2)) + 1i * (angle(-s1) + angle(s2));
    M0 = exp(-lm / 2);
    P = exp(-c / 2 - (Q.lch + log(abs(Q.D)) + 1i * ph) / 2 + lm / 2);
else
    % the quadratic D(s, t) = a2 lam^2 + a1 lam + a0, lam = -s, from its values at s = 0, -1, 1
    if model == 1
        d0 = Dfun(1, c, complex(0), t, true); d0z = d0;
    else
        d0z = Dfun(model, c, complex(0), t, false);
        d0 = d0z; if stat, d0 = d0z + Dstat_add(model, c, complex(0), t); end
    end
    dp = Dfun(model, c, complex(-1), t, stat);
    dm = Dfun(model, c, complex(1), t, stat);
    a2 = (dp + dm) / 2 - d0; a1 = (dp - dm) / 2; a0 = d0;
    disc = sqrt(a1 .* a1 - 4 * a2 .* a0);
    % the root of larger magnitude from the formula, the other from the product a0/a2 (the formula loses the small root
    % to cancellation when the roots differ by orders of magnitude, as at large c)
    sw = sign(real(conj(a1) .* disc)); sw(sw == 0) = 1;
    rb = (-a1 - sw .* disc) ./ (2 * a2); rs = a0 ./ (a2 .* rb);
    [s1, s2] = labelroots(-rb, -rs);
    M0 = 1 ./ (sqrt(-s1) .* sqrt(s2));
    ph0 = zeros(nt, 1);
    if any(~re)
        [aK, mu, e2, G] = argK(c, t(~re));
        DK = exp(c + mu) .* (1 + e2) / 2 .* G;            % e^c (cosh mu - c sinh(mu)/mu)
        ph0(~re) = aK;
        if model ~= 1, ph0(~re) = ph0(~re) + argpos(angle(d0z(~re) ./ DK), model == 3); end
        if stat
            if model == 1, ref = DK; else, ref = d0z(~re); end
            ph0(~re) = ph0(~re) + argneg(angle(d0(~re) ./ ref));
        end
    end
    P = abs(a2).^(-1/2) .* exp(-1i * (ph0 - angle(-s1) - angle(s2)) / 2);
end
P(neg) = conj(P(neg)); s1(neg) = conj(s1(neg)); s2(neg) = conj(s2(neg)); M0(neg) = conj(M0(neg));
end

% ---------- the convergence abscissa in t ----------
function ts = tstar(model, c, stat)
persistent CACHE
if isempty(CACHE), CACHE = containers.Map('KeyType', 'char', 'ValueType', 'double'); end
key = sprintf('%d|%.17g|%d', model, c, stat);
if isKey(CACHE, key), ts = CACHE(key); return; end
fac = 1.05;
if model == 4 && ~stat
    % the GLS-detrended trend case: the sign of cosh(th) det(I - Q Sig), scanned in steps of 0.5%
    g = @(t) gls_sign(t, c); fac = 1.005;
elseif model == 1 && ~stat
    % D(0,t) e^-c = S eta with S = sinh(mu)/mu: the sign is that of eta = mu coth mu - c, formed
    % without the cancellation of C - c S at large c (as in TPARTS), times the sign of S, which is
    % sin(th)/th for 2t > c^2.  eta alone changes sign at the poles of coth as well, and for c < 0
    % the first zero of D lies just below the pole (c^2 + pi^2)/2, inside a stretch that the 5%
    % steps of the bracketing can skip (c = -9.5: 61.6 for 49.2).
    g = @(t) eta1(c, t) * sgnS(c, t);
elseif ~stat && c > 2
    % cases DU, DT at large c: the Girsanov form det(I - Q Sig) of IIPAR_DFX, whose first
    % real zero is t_*, in place of D(0,t), whose terms of size e^(2c) cancel to a value of
    % order e^(-c) (at c = 20 D(0,t) is lost altogether; t_* is then about 1e-13)
    cas = 'DU'; if model == 3, cas = 'DT'; end
    g = @(t) real(getfield(iipar(complex(t), c, cas), 'D'));
else
    g = @(t) real(Dfun(model, c, complex(0), complex(t), stat));
end
if c < -20 && model ~= 4
    % for c < -20, t_* lies above c^2/2, where th = sqrt(c^2 - 2t) = i phi and the zeros of D(0,t) are
    % spaced by about pi in phi (the first at phi = 3.1, 6.3, 9.0 in cases K, U, T); the 5% steps in t
    % below (5% of c^2/2) skip pairs of them from c ~ -45 on, so the bracket is searched in phi from 0
    % in steps of 0.05 (as TSTAR in DFExact.jl; D(0,c^2/2), of order e^-|c| in cases DU, DT, underflows
    % beyond c ~ -745: an error there).  Not in case DG, whose first zero lies below c^2/2 from c ~ -100
    % on (and whose t_* serves c > 0 only)
    if ~(g(c * c / 2) > 0), error('dfx: D(0,t) is not positive at t = c^2/2 (|c| too large)'); end
    ph = 0;
    while true
        tn = (c * c + (ph + 0.05)^2) / 2;
        if ~(g(tn) > 0), a = (c * c + ph * ph) / 2; b = tn; break; end
        ph = ph + 0.05;
        if ph > 100, error('dfx: no zero of D(0,t) below t = (c^2 + 100^2)/2'); end
    end
else
    a = 1e-6;
    while ~(g(a) > 0) && a > 1e-300, a = a / 100; end      % (t_* is about 4e-7 at c = 10)
    if ~(g(a) > 0), error('dfx: D(0,t) not positive near t = 0'); end
    b = a;
    while g(b) > 0
        a = b; b = b * fac;
        if b > max(1e4, 2 * c * c), error('dfx: no zero of D(0,t) below t = max(1e4, 2 c^2)'); end
    end
end
for it = 1:100
    m = (a + b) / 2;
    if g(m) > 0, a = m; else, b = m; end
    if b - a < 1e-15 * b, break; end
end
ts = (a + b) / 2;
CACHE(key) = ts;
end

function s = sgnS(c, t)
% the sign of S = sinh(mu)/mu at real t (sin(th)/th for 2t > c^2)
q = c * c - 2 * t;
if q >= 0, s = 1; else, s = sign(sin(sqrt(-q))); end
end

function g = gls_sign(t, c)
[lch, D] = gls_parts(complex(t), c);
g = real(exp(lch) * D);
end

function e = eta1(c, t)
% eta = mu coth mu - c, mu = (c^2 - 2t)^(1/2), real t, without cancellation at large c
mu = sqrt(complex(c * c - 2 * t));
if abs(mu) < 1e-4
    e = real(1 + mu^2 / 3 - c);
    return
end
e2 = exp(-2 * mu);
e = real(-2 * t / (mu + c) + 2 * mu * e2 / (1 - e2));
end

function tau = default_tau(model, c, stat)
if c <= 0 && ~stat, tau = 1; else, tau = min(1, tstar(model, c, stat) / 2); end
end

% ---------- e^z K_0(z) ----------
function y = k0x(z)
% e^z K_0(z) for Re z >= 0: power series for |z| < 2, Taylor expansions about the centres of
% cells of side 1/2 for |z| < 25, the asymptotic series beyond (DFExact.k0x)
persistent K0A T CNT
if isempty(CNT), CNT = [0 0]; end
if isempty(z), y = CNT; CNT = [0 0]; return; end
CNT = CNT + [1 numel(z)];
if isempty(K0A)
    K0A = zeros(21, 1); K0A(1) = 1;
    for k = 1:20, K0A(k+1) = -K0A(k) * (2*k - 1)^2 / (8*k); end
    [I, J] = ndgrid(0:49, -50:49);
    z0 = complex((I(:).' + 0.5) * 0.5, (J(:).' + 0.5) * 0.5);
    T = complex(zeros(25, numel(z0)));
    T(1, :) = besselk(0, z0); T(2, :) = -besselk(1, z0);
    for n = 0:22
        if n >= 1, cm1 = T(n, :); else, cm1 = 0; end
        if n >= 2, cm2 = T(n-1, :); else, cm2 = 0; end
        T(n+3, :) = (z0 .* z0 .* T(n+1, :) + 2 * z0 .* cm1 + cm2 - 2 * z0 * (n+1) * n .* T(n+2, :) ...
            - z0 * (n+1) .* T(n+2, :) - n * (n-1) * T(n+1, :) - n * T(n+1, :)) ./ (z0 .* z0 * (n+2) * (n+1));
    end
end
y = complex(zeros(size(z)));
a = abs(z);
m0 = real(z) < 0;
if any(m0(:)), y(m0) = besselk(0, z(m0), 1); end
m1 = ~m0 & a >= 25;
if any(m1(:))
    zz = z(m1); w = 1 ./ zz; acc = K0A(21) * ones(size(zz));
    for k = 20:-1:1, acc = acc .* w + K0A(k); end
    y(m1) = sqrt(pi ./ (2 * zz)) .* acc;
end
m2 = ~m0 & a < 2;
if any(m2(:))
    zz = z(m2); q4 = zz .* zz / 4; term = ones(size(zz)); I0 = term; S = zeros(size(zz)); H = 0;
    for k = 1:14
        term = term .* q4 / (k * k); H = H + 1 / k;
        I0 = I0 + term; S = S + term * H;
    end
    y(m2) = exp(zz) .* (-(log(zz / 2) + 0.5772156649015329) .* I0 + S);
end
m3 = ~m0 & ~m1 & ~m2;
if any(m3(:))
    zz = z(m3);
    i = floor(real(zz) / 0.5); j = floor(imag(zz) / 0.5);
    idx = i + 1 + 50 * (j + 50);
    u = zz - complex((i + 0.5) * 0.5, (j + 0.5) * 0.5);
    acc = T(25, idx).'; acc = reshape(acc, size(zz));
    for n = 24:-1:1, acc = acc .* u + reshape(T(n, idx).', size(zz)); end
    y(m3) = exp(zz) .* acc;
end
end

% ---------- the inverse s-transform phi_t(X) ----------
function p = phi_pos(kern, X, s1, s2)
if kern == 1
    p = exp(-s1 * X) ./ sqrt(pi * X);
else
    p = k0x((s2 - s1) / 2 * X) .* exp(-s2 * X) / pi;
end
end
function p = phi_neg(kern, X, s1, s2)
if kern == 1
    p = zeros(size(X));
else
    p = k0x(-(s2 - s1) / 2 * X) .* exp(-s1 * X) / pi;
end
end
function p = phi_alg(kern, X, s1, s2)
if kern == 1
    p = 1 ./ sqrt(pi * X);
else
    p = k0x((s2 - s1) / 2 * X) / pi;
end
end

% ---------- quadrature ----------
function [y, wt] = glnodes(a, b, n, p)
% nodes and weights of the Gauss-Legendre rule on [a, b] with y = a + (b - a) u^p
n = n * nodefac();
[u, gw] = glrule(n);
L = b - a;
switch p
    case 1, y = a + L * u; dy = L * ones(n, 1);
    case 2, y = a + L * u.^2; dy = 2 * L * u;
    otherwise, y = a + L * u.^4; dy = 4 * L * u.^3;
end
wt = dy .* gw;
end

function acc = glsum(f, a, b, n, p)
[y, wt] = glnodes(a, b, n, p);
acc = sum(f(y) .* wt);
end

function n = nodes(n, osc, L)
% the node count of an outer piece of RAYINT, at most 0.85 rad per node; a rule beyond the bound of the engine (4096 nodes,
% GLRULE) raises the error dfx:bound:rule
m = 1.2 * osc * L;
if m > 4096, error('dfx:bound:rule', 'a Gauss-Legendre rule of the inner integrals needs %.3g nodes, at most 4096', m); end
n = max(n, 16 * ceil(ceil(m) / 16));
end

function acc = rayint(f, R, l, n1, n2, p, osc)
% the Gaussian ray from r = 0: a first piece on the scale l, then geometric panels of ratio 4
if nargin < 7, osc = 0; end
a = min(R, 8 * l);
[y, wt] = glnodes(0, a, n1, p);
while a < R && a > 0
    b = min(R, 4 * a);
    [y2, w2] = glnodes(a, b, nodes(n2, osc, b - a), 1);
    y = [y; y2]; wt = [wt; w2]; %#ok<AGROW>
    a = b;
end
acc = sum(f(y) .* wt);
end

function y = cap(x)
y = x ./ sqrt(1 + x .* x);
end

function acc = gaussexp(A, t, q, sd, lA, force, p)
% int_0^inf A(x) exp(-t x^2/q^2 - sd x) dx along a path on which neither factor grows
% (DFExact.gaussexp: the saddle path for Re z < 0, the level curve otherwise); p = 4 removes
% an algebraic or logarithmic singularity of A at x = 0
if nargin < 6, force = 0; end
if nargin < 7, p = 1; end
EXPO = 60;
at = abs(t); e = exp(-1i * angle(t) / 2);
x0 = -sd * q * q / (2 * t); u0 = -sd * sd * q * q / (4 * t);
if force == 1 || (force == 0 && real(sd * e) < 0)
    ru = real(u0);
    if ru > EXPO, lam_max = 1 - sqrt(1 - EXPO / ru); else, lam_max = 1; end
    f1 = @(lam) A(x0 * lam) .* exp(-u0 * lam .* (2 - lam)) * x0;
    l = min(lA / abs(x0), 1 / (2 * abs(u0) + 1));
    acc = rayint(f1, lam_max, l, 48, 32, p, 2 * abs(imag(u0)));
    if ru > EXPO, return; end
    f2 = @(r) A(x0 + r * e) .* exp(-u0 - at * r .* r / (q * q)) * e;
    acc = acc + rayint(f2, q * sqrt(EXPO / at), lA, 48, 32, 1, 0);
    return
end
v0 = abs(sd) * q / (2 * sqrt(at)); vA = sqrt(lA * abs(sd));
vmin = max(min([v0, vA, 1]) / 2, 1e-4);
S = asinh(sqrt(EXPO) / vmin);
if p == 1
    acc = glsum(@(s) levelf(s, A, vmin, sd, t, q), 0, S, 16 * ceil(ceil(13 * S) / 16), 1);
else
    S1 = min(1, S);
    acc = glsum(@(s) levelf(s, A, vmin, sd, t, q), 0, S1, 32, p);
    if S > 1, acc = acc + glsum(@(s) levelf(s, A, vmin, sd, t, q), 1, S, 16 * ceil(ceil(13 * (S - 1)) / 16), 1); end
end
end

function y = levelf(s, A, vmin, sd, t, q)
v = vmin * sinh(s); dv = vmin * cosh(s);
rt = sd * sqrt(1 + 4 * t * v .* v / (q * q * sd * sd));
x = 2 * v .* v ./ (sd + rt);
y = A(x) .* exp(-v .* v) .* 2 .* v ./ rt .* dv;
end

function ok = ray_ok(sd, e, q, at, cap)
% the Gaussian ray from v = 0 serves when the exponential grows by at most e^2 and turns by at most cap radians
% along it (60 by default, as in the Julia and C engines; RAYEXP uses 150 for the general kernel, see there)
EXPO = 60;
if nargin < 5, cap = 60; end
ok = max(0, -real(sd * e))^2 * q * q / (4 * at) <= 2 && abs(imag(sd * e)) * q * sqrt(EXPO / at) <= cap;
end

% ---------- the inner integrals ----------
function [Is, Ie] = inner_density(kern, tj, q, s1, s2, M0, part, e1) %#ok<INUSL>
% I_f(t) = int_0^inf 2 v^2 e^(-t v^2) phi_t(q v + 1/2) dv, split into the smooth part and, for
% q < 0, the edge part (the neighborhood of the support edge X = 0)
EXPO = 60; KAPPA = 6;
th = angle(tj) / 2; at = abs(tj); e = exp(-1i * th); w = abs(imag(tj));
vs = 1 / (2 * abs(q));
Is = complex(0); Ie = complex(0);
if q == 0                                                   % x = 0: I_f(t) = phi_t(1/2) int 2 v^2 e^(-t v^2) dv = (sqrt(pi)/2) t^(-3/2) phi_t(1/2), no edge part
    if part ~= 2, Is = sqrt(pi) / 2 * tj^(-1.5) * phi_pos(kern, 0.5, s1, s2); end
    return
end
if q > 0 || 3 * real(tj) * vs * vs / 16 > 40
    if part == 2, return; end
    if kern == 1, sd = s1; lA = 1 / (abs(s1) + 2); else, sd = s2; lA = 1 / (abs((s2 - s1) / 2) + 2); end
    if q < 0 || ray_ok(sd, e, q, at)
        if q > 0, R = sqrt(EXPO / at); else, R = min(vs, sqrt(EXPO / at)); end
        f = @(r) 2 * (r * e).^2 .* exp(-at * r .* r) .* phi_pos(kern, q * (r * e) + 0.5, s1, s2) * e;
        Is = rayint(f, R, lA / abs(q), 48, 32, 1, abs(q) * abs(imag(sd * e)));
        return
    end
    A = @(x) 2 * x .* x / q^3 .* phi_alg(kern, 0.5 + x, s1, s2) * exp(-sd / 2);
    Is = gaussexp(A, tj, q, sd, lA, 0);
    return
end
if imag(tj) >= 0, sg = 1; else, sg = -1; end
kap = KAPPA / (vs / 2); H = 1.5 / kap;
pathd = @(d) pathd_f(d, vs, H, sg, q, tj, kern, s1, s2);
chi = @(d) (1 + erfz(kap * (vs / 2 - d - 1i * H * cap((vs - d) / H) .* cap(d / H) * sg))) / 2;
if part ~= 2
    fs = @(r) pathd(vs - r) .* (1 - chi(vs - r));
    if w > 0, r1 = min(vs, sqrt(EXPO / (1.9 * w))); else, r1 = vs; end
    Is = glsum(fs, 0, r1, 64, 1);
    if r1 < vs, Is = Is + glsum(fs, r1, vs, 48, 1); end
end
if part == 1, return; end
fe = @(d) pathd(d) .* chi(d);
if w > 0, delta = min(vs, EXPO / (1.9 * w * vs)); else, delta = vs; end
if nargin < 8, e1 = Inf; end
Ie = edgehead(fe, delta / 8, e1 / abs(q)) + glsum(fe, delta / 8, delta, 64, 1);
if delta < vs, Ie = Ie + glsum(fe, delta, vs, 48, 1); end
if kern ~= 1
    rho = real((s2 - 3 * s1) / 2 * e) * abs(q);
    R = min(sqrt(EXPO / at), EXPO / (2 * at * vs * cos(th)));
    if rho > 0, R = min(R, EXPO / rho); end
    f = @(r) 2 * (vs + r * e).^2 .* exp(-tj * (vs + r * e).^2) .* phi_neg(kern, -abs(q) * r * e, s1, s2) * e;
    Ie = Ie + rayint(f, R, 1 / (rho + 2 * at * vs * cos(th)), 64, 32, 4, ...
        abs(q) * abs(imag(s1 * e)) + 2 * at * vs * abs(sin(th)));
end
end

function y = pathd_f(d, vs, H, sg, q, tj, kern, s1, s2)
% the bulged segment of inner_density, parametrized by d = v* - r (X = |q| d exact near the edge)
r = vs - d;
aa = cap(r / H); bb = cap(d / H);
g = H * aa .* bb;
dg = bb ./ sqrt(1 + (r / H).^2).^3 - aa ./ sqrt(1 + (d / H).^2).^3;
v = r - 1i * g * sg;
dv = 1 - 1i * dg * sg;
X = abs(q) * (d + 1i * g * sg);
y = 2 * v .* v .* exp(-tj * v .* v) .* phi_pos(kern, X, s1, s2) .* dv;
end

function y = pth_f(Xr, H, sg, q, tj, kern, s1, s2)
% the bulged segment (0, 1/2) of inner_cdf
aa = cap(Xr / H); bb = cap((0.5 - Xr) / H);
g = H * aa .* bb;
dg = bb ./ sqrt(1 + (Xr / H).^2).^3 - aa ./ sqrt(1 + ((0.5 - Xr) / H).^2).^3;
X = Xr + 1i * g * sg;
dX = 1 + 1i * dg * sg;
x = X - 0.5;
y = -phi_pos(kern, X, s1, s2) .* exp(-tj * x .* x / (q * q)) .* dX;
end

function [Is, Ie] = inner_cdf(kern, tj, q, s1, s2, M0, part, e1)
% I_F(t) = int phi_t(X) W_q(X - 1/2, t) dX without the factor 1/t
EXPO = 60; KAPPA = 6;
th = angle(tj) / 2; at = abs(tj); e = exp(-1i * th); w = abs(imag(tj));
Is = complex(0); Ie = complex(0);
vs2 = 1 / (4 * q * q);
if part ~= 2, Is = psi_half(kern, s1, s2, M0); end          % (i) the part of W_q that is constant in x: Psi_t(1/2)
if q == 0, return; end                                      % x = 0: I_F(t) = Psi_t(1/2), the Gaussian part vanishes and there is no edge part
if q > 0 || 3 * real(tj) * vs2 / 16 > 40
    if part == 2, return; end
    if kern == 1, sd = s1; lA = 1 / (abs(s1) + 2); else, sd = s2; lA = 1 / (abs((s2 - s1) / 2) + 2); end
    if q < 0 || ray_ok(sd, e, q, at)
        if q > 0, R = abs(q) * sqrt(EXPO / at); else, R = min(0.5, abs(q) * sqrt(EXPO / at)); end
        sq = sign(q);
        f = @(r) phi_pos(kern, 0.5 + sq * r * e, s1, s2) .* exp(-at * r .* r / (q * q)) * sq * e;
        Is = Is + rayint(f, R, lA, 48, 32, 1, abs(imag(sd * e)));
        return
    end
    A = @(x) phi_alg(kern, 0.5 + x, s1, s2) * exp(-sd / 2);
    Is = Is + gaussexp(A, tj, q, sd, lA, 0);
    return
end
if imag(tj) >= 0, sg = 1; else, sg = -1; end
kap = KAPPA / 0.25; H = 1.5 / kap;
pth = @(Xr) pth_f(Xr, H, sg, q, tj, kern, s1, s2);
chi = @(Xr) (1 + erfz(kap * (0.25 - Xr - 1i * H * cap(Xr / H) .* cap((0.5 - Xr) / H) * sg))) / 2;
if part ~= 2
    fs = @(y) pth(0.5 - y) .* (1 - chi(0.5 - y));
    if w > 0, ep = min(0.5, abs(q) * sqrt(EXPO / (1.9 * w))); else, ep = 0.5; end
    Is = Is + glsum(fs, 0, ep, 64, 1);
    if ep < 0.5, Is = Is + glsum(fs, ep, 0.5, 48, 1); end
end
if part == 1, return; end
fe = @(Xr) pth(Xr) .* chi(Xr);
if w > 0, delta = min(0.5, EXPO * q * q / (0.95 * w)); else, delta = 0.5; end
if nargin < 8, e1 = Inf; end
Ie = edgehead(fe, delta / 8, e1) + glsum(fe, delta / 8, delta, 64, 1);
if delta < 0.5, Ie = Ie + glsum(fe, delta, 0.5, 48, 1); end
if kern ~= 1
    rho = real((s2 - 3 * s1) / 2 * e);
    R = min(abs(q) * sqrt(EXPO / at), EXPO * q * q / (at * cos(th)));
    if rho > 0, R = min(R, EXPO / rho); end
    f = @(r) -phi_neg(kern, -r * e, s1, s2) .* exp(-tj * (-r * e - 0.5).^2 / (q * q)) * e;
    Ie = Ie + rayint(f, R, 1 / (rho + at * cos(th) / (q * q)), 64, 32, 4, ...
        abs(imag(s1 * e)) + at * abs(sin(th)) / (q * q));
end
end

% ---------- the outer integrals ----------
function v = filon_panel(G, a, b, y, N)
% int_a^b G(w) e^(-i y w) dw, G interpolated by Legendre polynomials at the N Gauss-Legendre
% nodes, the oscillation integrated exactly: int P_n(x) e^(-i Om x) dx = 2 (-i)^n j_n(Om)
persistent PM
if isempty(PM), PM = {}; end
if N > numel(PM) || isempty(PM{N})
    [u, gw] = glrule(N); x = 2 * u - 1; ww = 2 * gw;
    Pm = zeros(N, N); Pm(:, 1) = 1; if N > 1, Pm(:, 2) = x; end
    for n = 1:N-2, Pm(:, n+2) = ((2*n + 1) * x .* Pm(:, n+1) - n * Pm(:, n)) / (n + 1); end
    nc = (0:N-1).';
    PM{N} = (Pm .* ww).' .* ((2 * nc + 1) / 2) .* (2 * (-1i).^nc);    % maps G to the coefficients times 2(-i)^n
end
L = b - a; Om = y * L / 2;
nn = (0:N-1)';
if Om == 0
    jn = [1; zeros(N-1, 1)];
else
    jn = sqrt(pi / (2 * Om)) * besselj(nn + 0.5, Om);
end
v = L / 2 * exp(-1i * y * (a + b) / 2) * (jn.' * (PM{N} * G(:)));
end

function psi = saddle_phase(kern, q, ts, s1, s2)
if kern == 1, sd = s1; else, sd = s2; end
z = sd * q ./ (2 * sqrt(ts));
psi = nan(size(ts));
m = real(z) < 0 & real(z .* z) > -30;
psi(m) = imag(z(m) .^ 2);
end

function st = phase_steps(psi)
st = abs(diff(psi(:)));
st(isnan(st)) = 0;
end

function W = wmaxc(c)
% the end of the outer integrals: for c < -20 the kernel decays along the line like
% exp(-(Re sqrt(c^2 - 2 i w) - |c|)/2), more slowly, and the end moves to W with
% Re sqrt(d^2 - 2 i W) - d = sqrt(4900), d = -(c + 20) (unchanged on [-20, inf); as wmaxc in DFExact.jl)
W = 4900;
if c < -20
    d = -(c + 20); a = d + 70; W = max(W, a * sqrt(a * a - d * d));
end
end

function f = hufac(c)
% the factor of the step of the trapezoidal rule in u: 20/|c| for c < -20 (the phase of the kernel,
% about w/(2|c|), outruns its decay; as hufac in DFExact.jl)
f = 1;
if c < -20, f = 20 / (-c); end
end

function m = filonpanels(a, b, y, c)
% the number of Filon panels on [a, b], a panel of the geometric partition of an edge part with the
% frequency y: panels of at most 25/y up to w = 300/y; for c < -20 at most half the period 4 pi |c| of
% the factor e^(i w E V), E V ~ 1/(2|c|), of the edge integrand, which panels of 25/y ~ 50|c| span
% several times (unchanged on [-20, inf); as filonpanels in DFExact.jl)
if b <= 300 / y, m = ceil((b - a) * y / 25); else, m = 1; end
if c < -20, m = max(m, ceil((b - a) / (2 * pi * (-c)))); end
end

function e1 = k0edge1(model, c, stat, tau)
% the K_0 kernel for c < -20: the bound 1/b, b = (s_2 - s_1)/2 at t = tau, on the length of the first
% piece of the edge segment, the scale on which K_0(b X) turns from its logarithm to its exponential
% decay (with a zero start b ~ c^2/2 in case DU; the 48 nodes of a first piece of up to 1/16 lose 5e-10
% in F at c = -480); Inf on [-20, inf) and for the kernel of case DK (as k0first in DFExact.jl)
e1 = Inf;
if kernel(model, stat) == 2 && c < -20
    [~, s1, s2] = tparts(model, c, complex(tau), stat);
    e1 = 1 / abs(real(s2 - s1) / 2);
end
end

function acc = edgehead(f, a, e1)
% the first piece [0, a] of the edge segment, from the support edge, with the substitution u^4; with a
% bound e1 < a on its length (K0EDGE1) the rest of it in geometric panels of ratio 4 (EDGEHEAD in DFExact.jl)
if ~(e1 < a), acc = glsum(f, 0, a, 48, 4); return; end
acc = glsum(f, 0, e1, 48, 4); lo = e1;
while lo < a
    hi = min(a, 4 * lo); acc = acc + glsum(f, lo, hi, 48, 1); lo = hi;
end
end

function val = outer(model, c, q, tau, inner, weight, o, stat)
% f(q) = (1/pi) int_0^inf Re[P(t) I_f(t)] dw, P(tau <= q) = (1/pi) int_0^inf Re[P(t) I_F(t)/t] dw,
% t = tau + i w: the smooth part by the trapezoidal rule in u, w = tau sinh(u) (or Gauss-Legendre
% panels when the saddle factor oscillates), the edge part (q < 0) by Filon's method
WMAX = wmaxc(c);
hu = getf(o, 'hu', 0.075) * hufac(c); N = getf(o, 'N', 20); cde = tau;
NF = nodefac();
kern = kernel(model, stat);
umax = asinh(WMAX / cde);
us = (0:hu:umax)';
ts = tau + 1i * cde * sinh(us);
[P, s1, s2, M0] = tparts(model, c, ts, stat);
psi = saddle_phase(kern, q, ts, s1, s2);
if q > 0 && max(phase_steps(psi)) > 0.25
    w0 = min(0.5, tau / 2);
    edges0 = [0, w0 * 2.^(0:ceil(log2(WMAX / w0)))]; edges0(end) = WMAX;
    K = 8; wk = [];
    for p = 1:numel(edges0)-1
        tmp = linspace(edges0(p), edges0(p+1), K + 1); wk = [wk, tmp(1:end-1)]; %#ok<AGROW>
    end
    wk(end+1) = WMAX;
    te0 = tau + 1i * wk(:);
    [~, s1e, s2e] = tparts(model, c, te0, stat);
    steps = phase_steps(saddle_phase(kern, q, te0, s1e, s2e));
    wn = []; wt = [];
    for p = 1:numel(edges0)-1
        a = edges0(p); b = edges0(p+1);
        dp = sum(steps((p-1)*K + (1:K)));
        m = max(1, ceil(dp / 50));
        n = max(24, 16 * ceil(ceil(3 * dp / m) / 16)) * NF * getf(o, 'glfac', 1);
        [u, gw] = glrule(n);
        for k = 1:m
            aa = a + (b - a) * (k - 1) / m; bb = a + (b - a) * k / m;
            wn = [wn; aa + (bb - aa) * u]; wt = [wt; (bb - aa) * gw]; %#ok<AGROW>
        end
    end
    ts = tau + 1i * [0; wn];
    [P, s1, s2, M0] = tparts(model, c, ts, stat);
    gauss_check(kern, q, ts, s1, s2);
    gs = zeros(numel(wn), 1);
    for j = 1:numel(wn)
        Is = inner(ts(j+1), s1(j+1), s2(j+1), M0(j+1), 1);
        gs(j) = real(P(j+1) * Is * weight(ts(j+1))) * wt(j);
    end
    val = sum(gs) / pi;
else
    gauss_check(kern, q, ts, s1, s2);
    gs = zeros(numel(ts), 1);
    for j = 1:numel(ts)
        Is = inner(ts(j), s1(j), s2(j), M0(j), 1);
        gs(j) = real(P(j) * Is * weight(ts(j))) * cde * cosh(us(j));
    end
    gs(1) = gs(1) * 0.5;
    val = sum(gs) * hu / pi;
end
% the edge part
y = 1 / (4 * q * q);
if q > 0 || 3 * tau * y / 16 > 40, return; end
ehw = getf(o, 'edge_hw', []);
if ~isempty(ehw)
    ws = (0:ehw:WMAX)'; te = tau + 1i * ws;
    [P, s1, s2, M0] = tparts(model, c, te, stat);
    ge = zeros(numel(te), 1);
    for j = 1:numel(te)
        [~, Ie] = inner(te(j), s1(j), s2(j), M0(j), 2);
        ge(j) = real(P(j) * Ie * weight(te(j)));
    end
    ge(1) = ge(1) * 0.5;
    val = val + sum(ge) * ehw / pi;
    return
end
w0 = min([0.5, 1 / y, tau / 2]);
edges0 = [0, w0 * 2.^(0:ceil(log2(WMAX / w0)))]; edges0(end) = WMAX;
edges = 0;
for p = 1:numel(edges0)-1
    a = edges0(p); b = edges0(p+1);
    m = filonpanels(a, b, y, c);
    tmp = linspace(a, b, m + 1); edges = [edges, tmp(2:end)]; %#ok<AGROW>
end
pfac = getf(o, 'pfac', 1);
if pfac > 1
    e2 = [];
    for p = 1:numel(edges)-1
        tmp = linspace(edges(p), edges(p+1), pfac + 1); e2 = [e2, tmp(1:end-1)]; %#ok<AGROW>
    end
    edges = [e2, edges(end)];
end
u = glrule(N);
np = numel(edges) - 1;
wn = reshape(edges(1:np) + (edges(2:end) - edges(1:np)) .* u, [], 1);   % N x np, column-major
te = tau + 1i * [0; wn];
[P, s1, s2, M0] = tparts(model, c, te, stat);
G = complex(zeros(numel(wn), 1));
for j = 1:numel(wn)
    t = te(j+1);
    [~, Ie] = inner(t, s1(j+1), s2(j+1), M0(j+1), 2);
    G(j) = P(j+1) * Ie * weight(t) * exp(t * y);
end
acc = 0;
for p = 1:np
    acc = acc + filon_panel(G((p-1)*N+1:p*N), edges(p), edges(p+1), y, N);
end
val = val + real(exp(-tau * y) * acc) / pi;
end

function m = gauss_nmax(kern, t, q, s1, s2)
% The largest node count of a piece of the segment to the saddle that GAUSSEXP takes for the inner integrals of the
% central kernels at t and q > 0 (INNER_DENSITY, INNER_CDF), before rounding; 0 where they take another path.  Its phase
% 2 Im u0, u0 = -sd^2 q^2/(4t), grows like q^2, and RAYINT gives every piece 1.2 nodes per radian: at c = 7 in case DT
% the largest piece has more than 4096 nodes from q ~ 660 (gauss_nmax of DFExact.jl and the C engine)
EXPO = 60;
if kern == 1, sd = s1; lA = 1 / (abs(s1) + 2); else, sd = s2; lA = 1 / (abs((s2 - s1) / 2) + 2); end
at = abs(t); e = exp(-1i * angle(t) / 2);
m = 0;
if q <= 0 || ray_ok(sd, e, q, at) || real(sd * e) >= 0, return; end
x0 = -sd * q * q / (2 * t); u0 = -sd * sd * q * q / (4 * t);
ru = real(u0); if ru > EXPO, R = 1 - sqrt(1 - EXPO / ru); else, R = 1; end
l = min(lA / abs(x0), 1 / (2 * abs(u0) + 1)); osc = 2 * abs(imag(u0));
a = min(R, 8 * l);
while a < R
    b = min(R, 4 * a); m = max(m, 1.2 * osc * (b - a)); a = b;
end
end

function gauss_check(kern, q, ts, s1, s2)
% the error dfx:bound:rule before the smooth part of OUTER is evaluated when the inner integrals at one of its nodes ts
% need a rule beyond 4096 nodes (GAUSS_NMAX; only for q > 0, where GAUSSEXP serves)
if ~(q > 0), return; end
m = 0;
for j = 1:numel(ts), m = max(m, gauss_nmax(kern, ts(j), q, s1(j), s2(j))); end
if m > 4096
    error('dfx:bound:rule', ['the segment to the saddle of the Gaussian factor of the inner integrals, whose phase grows ' ...
          'like x^2, needs %.3g nodes, at most 4096'], m);
end
end

function v = uniform(model, c, tau, stat, inner, weight, hw)
% the trapezoidal rule in w for the whole integrand (DFExact.density_uniform, cdf_uniform)
WMAX = 4900;
w = (0:hw:WMAX)'; t = tau + 1i * w;
[P, s1, s2, M0] = tparts(model, c, t, stat);
g = zeros(numel(t), 1);
for j = 1:numel(t)
    [Is, Ie] = inner(t(j), s1(j), s2(j), M0(j), 0);
    g(j) = real(P(j) * (Is + sum(Ie)) * weight(t(j)));
end
g(1) = g(1) * 0.5;
v = sum(g) * hw / pi;
end

% ---------- moments ----------
function m = moments(model, c, o)
% E tau^p = Gamma(p/2)^(-1) int_0^inf s^(p/2-1) m_p(s) ds, m_p(s) = d^p/du^p M(u,-s) at u = 0, by
% Cauchy integrals on a circle of radius 0.4 min|s_k(-s)| (N = 64) and the trapezoidal rule in log s
stat = startflag(getf(o, 'start', 'zero'), c);
if model == 4, error('dfx: the moments are not available for the GLS-detrended case (model 4)'); end
pmax = getf(o, 'pmax', 4); h = getf(o, 'h', 0.05); N = getf(o, 'Nc', 64); rfrac = getf(o, 'rfrac', 0.4);
if model == 1 && ~stat, extra = -c / 2; else, extra = 0; end
kern = kernel(model, stat);
ls = (-60:h:8)';
acc = zeros(pmax, 1);
th = 2 * pi * (0:N-1)' / N;
for il = 1:numel(ls)
    s = exp(ls(il)); t = complex(-s);
    [~, s1, s2] = tparts(model, c, t, stat);
    if kern == 1, rad = rfrac * abs(s1); else, rad = rfrac * min(abs(s1), abs(s2)); end
    u = rad * exp(1i * th);
    d = Dfun(model, c, u, t * ones(N, 1), stat);
    ph = unwrapc(angle(d));
    M = exp(-u / 2 + extra) .* abs(d).^(-0.5) .* exp(-1i * ph / 2);
    for p = 1:pmax
        mp = real(sum(M .* exp(-1i * p * th))) / N / rad^p * factorial(p);
        acc(p) = acc(p) + s^(p / 2) * mp;
    end
end
m = acc * h ./ gamma((1:pmax)' / 2);
end

% =====================================================================================
% The local null value c0: the general kernel
%
% The t-ratio for phi0 = 1 + c0/n has the limit (U - c0 V)/sqrt(V).  In the derivation of
% the line integral representation only the region changes: {U <= x sqrt(V)} becomes {U <= x sqrt(V) + c0 V}, so with
% V = v^2 the kernel runs along the parabola u = x v + c0 v^2 instead of the line u = x v,
% with the same Jacobian 2 v^2.  With X(v) = c0 v^2 + x v + 1/2 the transforms are
%   density:  I_f(t) = int_0^inf 2 v^2 e^(-t v^2) phi_t(X(v)) dv,
%   cdf:      I_F(t) = int_0^inf 2 t v e^(-t v^2) Psi_t(X(v)) dv
%                    = Psi_t(1/2) + int_0^inf e^(-t v^2) (x + 2 c0 v) phi_t(X(v)) dv,
% Psi_t the distribution function of phi_t, and f = (1/2 pi i) int P I_f dt, F = (1/2 pi i)
% int P I_F dt/t as before.  The support edge X = 0 of phi_t (case DK), or its logarithmic
% singularity, is met at the positive roots of X(v): none, one or two of them.  Each root b
% contributes a part carrying exp(-t b^2), which oscillates in Im t with the period 2 pi/b^2
% and gets its own Filon part in the outer integral; for c0 < 0 there is always one root,
% also for x > 0.  Between the roots X < 0 (phi_neg; nothing in case DK).  For c0 > 0 and
% x < 0 without real roots, the vertex v0 = -x/(2 c0) of the parabola is a breakpoint too: the
% complex roots v0 +- i l of X are at the distance l = sqrt(Xmin/c0) from the real axis, and
% the paths of the quadrature, which leave the real axis, must not cross the cuts of
% X^(-1/2) or log X that run from them away from the real axis; the segment [0, v0] is
% bulged and the ray starts at v0, as at an edge, and the part near v0 carries exp(-t v0^2).
% Along each ray from a breakpoint b the exponent of the kernel and of phi_t is
%   -t v^2 - s_d X(v) = -t b^2 - s_d X(b) - beta r - alpha r^2,   v = b + r,
% alpha = t + c0 s_d, beta = 2 alpha b + s_d x, s_d = s_1 (X < 0, or case DK) or s_2, so the
% ray and saddle paths of GAUSSEXP apply with (alpha, beta) in place of (t/q^2, s_d).

function R = breakpts(q, c0, tau)
% The breakpoints of the v-integrals of INNER_GEN for the kernel along u = q v + c0 v^2:
% X(v) = c0 v^2 + q v + 1/2 = c0 (v - r1)(v - r2) (or q (v - r1) for c0 = 0, cancellation-free
% near the roots), its positive roots b (edges, kindb = 1) and, for c0 > 0, q < 0 without real
% roots, the vertex (kindb = 0); y = b^2 the frequencies of the parts; ell the scale of the
% singularity near each breakpoint (the distance of the nearest other root, real or complex),
% ell0 the same at v = 0;
% sgn the sign of X on the segments [0, b_1], ..., [b_ne, inf).
R.q = q; R.c0 = c0;
if c0 == 0
    rr = -1 / (2 * q); lead = q; r = rr;
else
    d = q * q - 2 * c0;
    if d >= 0, sq = sqrt(d); else, sq = 1i * sqrt(-d); end
    rr = [(-q - sq) / (2 * c0), (-q + sq) / (2 * c0)]; lead = c0;
    if d >= 0, r = sort(real(rr)); else, r = []; end
end
% X at v = b + dv, cancellation-free near a root b: lead * prod_k (dv + (b - r_k))
R.rr = rr; R.Xat = @(dv, b) xat(dv, b, rr, lead);
b = r(r > 0); b = b(:).'; noroot = isempty(b);
% a root so far out that e^(-t v^2) is below e^-40 there on the line Re t = tau (3 tau b^2/4 > 40)
% is dropped for c0 = 0 only: the part near it is negligible and the segment up to it would be
% integrated at a cost that grows with b (for q -> 0- the root -1/(2 q) runs off to infinity).
% With c0 ~= 0 a root is kept however far out it is: X changes sign there and the segment beyond
% it has the other form (phi_neg instead of phi_pos, or the reverse) and the other exponential
% factor (s1 instead of s2), whose Gaussian decay alpha = t + c0 s_d the form of the near side
% does not have.  Dropping it ran the near-side form out to infinity: Re alpha <= 0 on the ray
% (an error), or K_0 outside its domain.  (November 2026; ltratio_pdf(9,'SU',-1.2,
% struct('c0',-1.2)), and every x at c0 = -0.05, where the root q/|c0| is far out.)  The edge
% parts of far roots are skipped as negligible by ONEVALUE (R.negl).
if c0 == 0, b = b(3 * tau * b.^2 / 4 <= 40); end
kindb = ones(size(b)); ell = inf(size(b)); R.double = false;
if numel(b) == 2 && b(2) - b(1) <= 1e-9 * b(2)              % a double root (x = -sqrt(2 c0), c0 > 0)
    b = mean(b); kindb = 1; ell = Inf; R.double = true;
elseif numel(b) == 2
    ell(:) = b(2) - b(1);
end
if c0 > 0 && q < 0 && noroot
    b = -q / (2 * c0); kindb = 0; ell = sqrt(max(1/2 - q * q / (4 * c0), 0) / c0);
end
R.dX = @(v) 2 * c0 * v + q;
for k = 1:numel(b)                                        % the distance of the nearest other root
    dr = abs(b(k) - rr); dr = dr(dr > 1e-9 * max(b(k), 1));
    if ~isempty(dr), ell(k) = min(ell(k), min(dr)); end
end
ell = max(ell, 1e-8 * max(b, 1));                         % (a double root: the density is infinite there)
R.ell0 = min(abs(rr));                                    % the distance of the nearest root from v = 0
R.b = b; R.y = b.^2; R.kindb = kindb; R.ell = ell; R.ne = numel(b);
ed = [0 b Inf]; R.sgn = ones(1, numel(b) + 1);
for j = 1:numel(b) + 1
    if isfinite(ed(j+1)), vm = (ed(j) + ed(j+1)) / 2; else, vm = ed(j) + 1; end
    sv = sign(real(R.Xat(vm, 0))); if sv ~= 0, R.sgn(j) = sv; end
end
R.negl = false(1, numel(b));
end

function X = xat(dv, b, rr, lead)
X = lead * ones(size(dv));
for k = 1:numel(rr), X = X .* (dv + (b - rr(k))); end
end

function p = phialg(kern, X, sgn, s1, s2)
% phi_t without its exponential factor e^(-s_d X): algebraic, bounded away from X = 0
if kern == 1
    p = 1 ./ sqrt(pi * X);
elseif sgn > 0
    p = k0x((s2 - s1) / 2 * X) / pi;
else
    p = k0x(-(s2 - s1) / 2 * X) / pi;
end
end

function p = phisg(kern, X, sgn, s1, s2)
if sgn > 0, p = phi_pos(kern, X, s1, s2); else, p = phi_neg(kern, X, s1, s2); end
end

function Is = psi_half(kern, s1, s2, M0)
% Psi_t(1/2) = int_{-inf}^{1/2} phi_t(X) dX: closed form in case DK, otherwise the total mass
% minus a short integral over X > 1/2 (or the integral over X < 1/2 when |s2| is small)
if kern == 1
    Is = erfz(sqrt(s1 / 2)) / sqrt(s1); return
end
l = 1 / (abs(s2 - s1) + 2);
if abs(s2) >= abs(s1) / 4
    e2 = exp(-1i * angle(s2));
    Is = M0 - rayint(@(r) phi_pos(kern, 0.5 + r * e2, s1, s2) * e2, 36 / abs(s2), l, 40, 32, 1, 0);
else
    e1 = exp(-1i * angle(-s1));
    Is = rayint(@(X) phi_pos(kern, X, s1, s2), 0.5, l, 40, 32, 4, 0) + ...
        rayint(@(r) phi_neg(kern, -r * e1, s1, s2) * e1, 36 / abs(s1), l, 40, 32, 4, 0);
end
end

function val = rayexp(A, alpha, beta, lA, p)
% int_0^inf A(r) exp(-alpha r^2 - beta r) dr, Re alpha > 0: along the ray r e^(-i arg(alpha)/2),
% on which the Gaussian decays, while the exponential grows by at most e^2 and turns by at
% most 150 rad along it; otherwise by the saddle or level-curve path of GAUSSEXP.  The cap is
% 150 rad, not the 60 of RAY_OK: the level-curve path, taken when the saddle is absent, passes
% the branch point of its square root at |z| with the weight e^(-|z|^2); at 60 rad it took over
% at |Im z| ~ 4, where that is 1e-7, and the inner transform lost seven digits on a short
% stretch of the line (1e-10 in the distribution function in case DK with c0 > 0 at x >= 0);
% at 150 rad the level curve begins at |Im z| ~ 10, where the branch point is below 1e-40
% (the same change as in the Julia and C engines, November 2026).
EXPO = 60;
if ~(real(alpha) > 0), error('dfx: Re alpha <= 0 on a ray'); end
e = exp(-1i * angle(alpha) / 2); aa = abs(alpha);
if max(0, -real(beta * e))^2 / (4 * aa) <= 2 && abs(imag(beta * e)) * sqrt(EXPO / aa) <= 150
    Rr = sqrt(EXPO / aa);
    if real(beta * e) > 0, Rr = min(Rr, EXPO / real(beta * e)); end
    f = @(r) A(r * e) .* exp(-aa * r .* r - beta * e * r) * e;
    val = rayint(f, Rr, lA, 48, 32, p, abs(imag(beta * e)));
else
    val = gaussexp(A, alpha, 1, beta, lA, 0, p);
end
end

function [Is, Ie] = inner_gen(kern, tj, R, s1, s2, M0, kind, part)
% The transforms of the general kernel (see above), kind 1 (density) or 2 (cdf).  Is: the
% smooth part (Psi_t(1/2) for the cdf, the bump of e^(-t v^2) at v = 0, and the ray from 0
% when there is no breakpoint); Ie(k): the part near the breakpoint R.b(k), carrying
% exp(-t b_k^2), from the segment(s) adjacent to it and, for the last one, the ray beyond.
% part: 0 all, 1 the smooth part only, 2 all edge parts, 10 + k the edge part k only.
EXPO = 60; KAPPA = 6;
persistent NCALL
if isempty(NCALL), NCALL = zeros(1, 3); end
if nargin == 0, Is = NCALL; NCALL = zeros(1, 3); return; end
NCALL(min(part, 2) + 1) = NCALL(min(part, 2) + 1) + 1;
q = R.q; c0 = R.c0; ne = R.ne; b = R.b;
w = abs(imag(tj));
if imag(tj) >= 0, sg = 1; else, sg = -1; end
Is = complex(0); Ie = complex(zeros(1, ne));
if kind == 1, W = @(v) 2 * v .* v; else, W = @(v) q + 2 * c0 * v; end
wants = part == 0 || part == 1;
want = @(k) (part == 0 || part == 2 || part == 10 + k) && ~R.negl(k);
if kind == 2 && wants, Is = psi_half(kern, s1, s2, M0); end
if kern == 1, lX = 1 / (abs(s1) + 2); else, lX = 1 / (abs((s2 - s1) / 2) + 2); end   % scale of the algebraic factor in X
ed = [0 b Inf];
for j = 1:ne + 1
    a = ed(j); bb = ed(j+1); sgn = R.sgn(j);
    if kern == 1 && sgn < 0, continue; end
    if kern == 1 || sgn < 0, sd = s1; else, sd = s2; end
    if isinf(bb)
        % the ray from a (a = 0: smooth; a = b_ne: edge ne)
        if j == 1, if ~wants, continue; end, else, if ~want(j - 1), continue; end, end
        alpha = tj + c0 * sd; beta = 2 * alpha * a + sd * q;
        cst = exp(-tj * a * a - sd * R.Xat(0, a));
        A = @(rho) W(a + rho) .* phialg(kern, R.Xat(rho, a), sgn, s1, s2);
        lv = lX / max(abs(R.dX(a)), 1e-300);
        if c0 ~= 0, lv = min(lv, sqrt(lX / abs(c0))); end
        lv = min(lv, R.ell0 / 16);                        % (RAYINT's first piece is 8 lv)
        p = 1;
        if j > 1 && R.kindb(j - 1) == 0
            % from a pseudo-edge: the straight ray in the direction e of the Gaussian's decay (within
            % 45 degrees of the real axis, so at least 0.7 ell from the branch points v0 -+ i ell,
            % whose cuts the level-curve path of GAUSSEXP would cross), the oscillation of
            % e^(-beta rho e) along it resolved by PIECES_NODES
            e = exp(-1i * angle(alpha) / 2); aa = abs(alpha); be = beta * e;
            Rr = (-real(be) + sqrt(real(be)^2 + 4 * aa * EXPO)) / (2 * aa);       % fallen by EXPO
            [s, w] = pieces_nodes(@(r) -aa * r .* r - be * r, Rr, min(Rr / 8, R.ell(j - 1)), 1, 48);
            val = cst * sum(A(s * e) .* exp(-aa * s .* s - be * s) .* w) * e;
        else
            if j > 1, lv = min(lv, R.ell(j - 1) / 16); p = 4; end
            val = cst * rayexp(A, alpha, beta, lv, p);
        end
        if j == 1, Is = Is + val; else, Ie(j - 1) = Ie(j - 1) + val; end
    else
        % a finite segment [a, bb]: SEGINT
        L = bb - a; h = L / 2; kap = KAPPA / h; m = a + h;
        dep = L / 8;
        if j == 1 && R.vsad > 0 && R.vsad < bb                % (R.vsad is NaN when the saddle is irrelevant)
            % the saddle of the exponent inside the first segment (c0 < 0, large q): the partition
            % is placed midway between the saddle and bb, and sharpened, so that the saddle's
            % oscillation, e^(i psi) with psi ~ q^2 Im t, which the smooth part's rule follows,
            % leaks into the edge part by at most erf(-6), and the edge into the smooth part
            % likewise; the depth of the paths is reduced with it, keeping |chi| <= e^2.25
            m = max(m, (R.vsad + bb) / 2); kap = max(kap, KAPPA / (bb - m));
            kap = min(kap, 4 * KAPPA / L); m = min(m, bb - KAPPA / kap);     % at most twice as sharp
            dep = min(dep, KAPPA / 4 / kap);
        end
        chi = @(v) (1 + erfz(kap * (v - m))) / 2;                       % 1 near bb, 0 near a
        pa = 1; if j > 1, pa = 4; end
        ws_ = j == 1 && wants; wl_ = j > 1 && want(j - 1); wr_ = want(j);
        if ~(ws_ || wl_ || wr_), continue; end
        ella = R.ell0; if j > 1, ella = R.ell(j - 1); end
        ka = 1; if j > 1, ka = R.kindb(j - 1); end           % the kinds of the two ends (0: a pseudo-edge)
        [vs_, vl_, vr_] = segint(kern, tj, sd, sgn, s1, s2, W, R, a, bb, chi, pa, 1 + 3 * R.kindb(j), ella, R.ell(j), ws_, wl_, wr_, sg, ka, R.kindb(j), dep);
        if ws_, Is = Is + vs_; end
        if wl_, Ie(j - 1) = Ie(j - 1) + vl_; end
        if wr_, Ie(j) = Ie(j) + vr_; end
    end
end
end

function [vs, vl, vr] = segint(kern, tj, sd, sgn, s1, s2, W, R, a, b, chi, pa, pb, ella, ellb, ws, wl, wr, sg, ka, kb, dep)
% The integral over the segment [a, b] of W(v) e^(-t v^2) phi_t(X(v)), split by the partition
% chi into the part near b (vr) and the part near a (vs: the bump at a = 0, or vl: the edge at
% a).  The path has the fixed endpoints a and b (so that each part is independent of it and
% analytic in t) and three straight pieces: from a in the direction of steepest descent of the
% exponent E = -t v^2 - s_d X(v) at a, from b likewise, both kept within 80 degrees of the
% segment and within the depth L/8 (|chi| <= e^(2.25) there) and stopped where the quadratic
% term turns the exponent up, and the straight connection between their ends.  Near an end
% where e^(-s_d X) is still large and the Gaussian oscillates, the descent direction may point
% into the other half plane than the Gaussian's; the connection then crosses the real axis
% where both factors are small.  The nodes come from PIECES_NODES (geometric pieces, node
% counts from the phase of the exponent, stopping where the exponent is below -60) and, for
% the connection, PHASE_NODES; the integrand and the partition are evaluated once at all of
% them.
EXPO = 60;
L = b - a;
Xat = R.Xat; c0 = R.c0;
% the integrand at the offset z from the end o (a or b), X formed from z without cancellation
Ez = @(z, o) -tj * (o + z) .* (o + z) - sd * Xat(z, o);
Fz = @(z, o) W(o + z) .* exp(Ez(z, o)) .* phialg(kern, Xat(z, o), sgn, s1, s2);
Eb = Ez(0, b); Ea = Ez(0, a);
alpha = tj + c0 * sd;
% the descent directions, kept within 80 degrees of the segment, or 45 degrees at a pseudo-edge,
% whose branch points v0 -+ i ell have their cuts running vertically away from the real axis;
% the first pieces on the scale of the nearest singularity (ell/2 when it lies on the real axis,
% ell at a pseudo-edge, where it is at the perpendicular distance >= 0.7 ell from the path)
dega = 80; fa = ella / 2; if ka == 0, dega = 45; fa = ella; end
degb = 80; fb = ellb / 2; if kb == 0, degb = 45; fb = ellb; end
[ua, la] = descent(-2 * tj * a - sd * R.dX(a), alpha, 1, L, dep, sg, EXPO, dega);
[ub, lb] = descent(-2 * tj * b - sd * R.dX(b), alpha, -1, L, dep, sg, EXPO, degb);
% a descent piece that climbs above the larger end of the segment is cut back (NOCLIMB): its terms would cancel against
% the connection's
if isnan(real(Ea)) || isnan(real(Eb)), Etop = NaN; else, Etop = max(real(Ea), real(Eb)) + 4; end
la = noclimb(@(z) Ez(z * ua, a), la, Etop, L);
lb = noclimb(@(z) Ez(z * ub, b), lb, Etop, L);
pa_ = a + la * ua; pb_ = b + lb * ub;
% the three pieces: the integrand at the nodes (from the offsets along the two descent pieces,
% so that X is exact near the ends) and the partition
[s, w] = pieces_nodes(@(s) Ez(s * ua, a) - Ea, la, min(la / 8, fa), pa, 48 + 16 * (pa == 1));
fv = Fz(s * ua, a) .* (w * ua); v = a + s * ua;
[s, w] = pieces_nodes(@(s) Ez(s * ub, b) - Eb, lb, min(lb / 8, fb), pb, 48);
fv = [fv; -Fz(s * ub, b) .* (w * ub)]; v = [v; b + s * ub];     % from b toward the segment: reversed
% the connection: the straight line, unless the exponent rises along it above the ends (the
% Gaussian grows in the half plane a descent piece has entered under the exponential factor,
% as with Re alpha ~ 0 and a steep -beta v); then the polyline through the real axis, pa_ ->
% Re pa_ -> Re pb_ -> pb_, on which the exponent is that of the integrand itself
if abs(pb_ - pa_) > 0
    Emax = max(real(Ea), real(Eb)) + 2;
    legs = [pa_, pb_];
    if max(real(Ez(pa_ + (0:32)' / 32 * (pb_ - pa_) - a, a))) > Emax
        legs = [pa_, real(pa_), real(pb_), pb_];
    end
    for k = 1:numel(legs) - 1
        dv = legs(k+1) - legs(k);
        if abs(dv) == 0, continue; end
        [s, w] = phase_nodes(@(s) Ez(legs(k) + s * dv - a, a), Ea, Eb, EXPO);
        vc = legs(k) + s * dv;
        fv = [fv; Fnear(vc, a, b, Fz) .* (w * dv)]; v = [v; vc]; %#ok<AGROW>
    end
end
cw = chi(v);
vs = complex(0); vl = complex(0); vr = complex(0);
if ws || wl
    val = sum(fv .* (1 - cw));
    if ws, vs = val; else, vl = val; end
end
if wr, vr = sum(fv .* cw); end
end

function l = noclimb(Ef, l, Etop, L)
% The length l of a descent piece of SEGINT, halved until the real part of its exponent Ef (at 8 points along it) is at
% most Etop, the larger exponent of the two ends of the segment plus 4; unchanged where it is (a descent piece normally
% falls).  Where the exponent rises into the segment at an end, the steepest descent points out of the segment and
% DESCENT clamps it to 80 degrees from the segment, where e^(-t v^2) grows like e^(Re(t) Im(v)^2); in the Julia and C
% engines the noncentral kernel of a fixed start (which these programs do not have) climbed so to e^213 and cancelled
% against the connection (F off by 6e74 in case DK at c = 0, c0 = -0.01, g0 = 1, x = 1.5; noclimb of DFExact.jl)
while l > 1e-12 * L
    E = real(Ef(l * (1:8)' / 8));
    if any(isnan(E)) || ~(max(E) > Etop), break; end
    l = l / 2;
end
end

function y = Fnear(v, a, b, Fz)
% the integrand at v, from the nearer end of [a, b]
nb = abs(v - b) < abs(v - a); y = complex(zeros(size(v)));
if any(~nb), y(~nb) = Fz(v(~nb) - a, a); end
if any(nb), y(nb) = Fz(v(nb) - b, b); end
end

function v = tune(varargin)
% the node counts of the segment pieces: [NGROW NLEV NCON] (the least nodes of a growing piece,
% the grading levels of the connection, the least nodes of each of its sub-pieces)
persistent T
if isempty(T), T = [16 5 16]; end
if nargin > 0, T = varargin{1}; end
v = T;
end

function n = NODEMAX()
% The most nodes on one path of the general kernel's inner integrals (PIECES_NODES, PHASE_NODES; NODEMAX of DFExact.jl,
% DFL_NODEMAX of the C engine).  The paths follow the phase of the exponent -t v^2 - s_d X(v) with about 1.2 nodes per
% radian, and on a long segment between the support edges (a root b of X far out, as with c0 < 0 near zero at x > 0,
% where b ~ x/|c0|) that phase grows like Im(t) b^2: about 6e11 nodes at c0 = -1e-4, x = 1.  Beyond NODEMAX the path is
% refused (the error dfx:bound:path) before anything is evaluated on it.
n = 2e6;
end

function path_refuse(n)
error('dfx:bound:path', 'a path of the inner integrals of the general kernel needs %.3g nodes, at most 2000000', n);
end

function [y, wt] = pieces_nodes(Ef, L, first, p, nmin)
% nodes and weights for int_0^L: the first piece [0, first] (shortened, by halving, until the
% exponent falls by at most 4 and its phase changes by at most 40 rad across it) with nmin
% nodes and the substitution p, then pieces growing by the factor 4, with the node counts set
% by the change of the phase of the exponent across each piece (at most 0.85 rad per node, in
% sub-pieces of at most 160 nodes), until the exponent is below -EXPO
EXPO = 60; T = tune(); NGROW = T(1); NF = nodefac();
first = min(first, L);
for k = 1:60
    E = Ef(first * (0:8)' / 8);
    if (real(E(end)) > -4 && sum(abs(diff(imag(E)))) < 40) || first <= 1e-12 * L, break; end
    first = first / 2;
end
[y, wt] = glnodes(0, first, nmin, p);
tot = nmin * NF;                                          % the nodes of the path so far (at most NODEMAX)
a = first;
while a < L
    bn = min(L, 4 * a);
    dd = a + (bn - a) * (0:8)' / 8; E = Ef(dd);
    if all(real(E) < -EXPO), break; end
    dph = sum(abs(diff(E)));                              % the change of the exponent (phase and modulus)
    if tot + (1.2 * dph + 16) * NF > NODEMAX, path_refuse(tot + (1.2 * dph + 16) * NF); end
    n = max(NGROW, 16 * ceil((1.2 * dph + 16) / 16));
    m = ceil(n / 160); n = 16 * ceil(n / m / 16);                % at most 160 nodes per sub-piece
    tot = tot + m * n * NF;
    if tot > NODEMAX, path_refuse(tot); end
    for k = 1:m
        [y2, w2] = glnodes(a + (bn - a) * (k - 1) / m, a + (bn - a) * k / m, n, 1);
        y = [y; y2]; wt = [wt; w2]; %#ok<AGROW>
    end
    if real(E(end)) < -EXPO, break; end
    a = bn;
end
end

function [y, wt] = phase_nodes(Ef, E1, E2, EXPO)
% nodes and weights for int_0^1: Gauss-Legendre in sub-pieces graded geometrically toward both
% ends (NLEV levels, the shortest 2^-NLEV), with the node counts from the change of the
% exponent (phase and modulus, sampled at 9 points of each sub-piece, so that a stationary
% point of the phase inside is seen), skipping sub-pieces on which the exponent is below
% max(E1, E2) - EXPO (up to 2: the descent pieces end where the exponent has fallen by EXPO),
% E1 and E2 the exponents at the two ends of the segment
T = tune(); NLEV = T(2); NCON = T(3); NF = nodefac();
g = 2 .^ (-NLEV:-1); ss = [0 g 1 - g(end-1:-1:1) 1]'; ns = numel(ss) - 1;
Emax = max(real(E1), real(E2)) - EXPO + 2;
y = zeros(0, 1); wt = zeros(0, 1); tot = 0;
for k = 1:ns
    E = Ef(ss(k) + (ss(k+1) - ss(k)) * (0:8)' / 8);
    if max(real(E)) < Emax, continue; end
    dph = sum(abs(diff(E)));
    if tot + (1.2 * dph + 16) * NF > NODEMAX, path_refuse(tot + (1.2 * dph + 16) * NF); end
    n = max(NCON, 8 * ceil((1.2 * dph + 16) / 8));          % (+16: a stationary point of the phase inside)
    m = ceil(n / 160); n = 8 * ceil(n / m / 8);
    tot = tot + m * n * NF;
    if tot > NODEMAX, path_refuse(tot); end
    for i = 1:m
        [y2, w2] = glnodes(ss(k) + (ss(k+1) - ss(k)) * (i - 1) / m, ss(k) + (ss(k+1) - ss(k)) * i / m, n, 1);
        y = [y; y2]; wt = [wt; w2]; %#ok<AGROW>
    end
end
end

function [u, l] = descent(gam, alpha, dir, L, dep, sg, EXPO, deg)
% the unit direction u of steepest descent of E(v0 + z) = E(v0) + gam z - alpha z^2 at an end,
% kept within deg degrees of the direction dir (+1: into the segment from its left end, -1:
% from its right end) and, when gam = 0, vertical toward -i sg; and the length l along it:
% where the exponent has fallen by EXPO, or the turning point of its quadratic part, or the
% depth dep, or the half-length L/2 of the segment
if abs(gam) > 0, u = -conj(gam) / abs(gam); else, u = -1i * sg; end
cmax = cos(deg * pi / 180);
if real(u) * dir < cmax                                  % clamp to deg degrees from dir
    si = sign(imag(u)); if si == 0, si = -sg; end
    u = dir * cmax + 1i * si * sin(deg * pi / 180);
end
g1 = real(gam * u); g2 = -real(alpha * u * u);           % Re E(v0 + s u) - E(v0) = g1 s + g2 s^2
l = L / 2;
if abs(imag(u)) > 0, l = min(l, dep / abs(imag(u))); end
if g1 < 0
    if g2 > 0, l = min(l, -g1 / (2 * g2)); end            % the turning point of the parabola
    disc = g1 * g1 - 4 * g2 * EXPO;                       % g2 s^2 + g1 s + EXPO = 0: fallen by EXPO
    if g2 == 0
        l = min(l, -EXPO / g1);
    elseif disc >= 0
        sE = (-g1 - sqrt(disc)) / (2 * g2); if sE > 0, l = min(l, sE); end
    end
end
l = max(l, 1e-12 * L);
end

function v = lastparts(x)
% debug: the smooth part of the last OUTER_GEN value
persistent V
if nargin > 0, V = x; end
v = V;
end

function vs = saddle_loc(model, c, q, c0, tau, stat, kern, R)
% the location Re v* of the saddle of the exponent -t v^2 - s_d X(v) on the first segment,
% v* = -s_d q / (2 (t + c0 s_d)), at the reference point t = tau + i (it moves little with t
% while the saddle matters); NaN when there is no segment
vs = NaN;
if R.ne == 0, return; end
[~, s1, s2] = tparts(model, c, [tau; tau + 1i], stat);
if kern == 1 || R.sgn(1) < 0, sd = s1(2); else, sd = s2(2); end
vs = -sd * q / (2 * (tau + 1i + c0 * sd));
% relevant when it lies inside the first segment and close to the real axis (with c0 > 0 it
% can have a positive real part but sits far below the axis)
if ~(real(vs) > 0 && real(vs) < R.b(1) && abs(imag(vs)) < 0.2 * (R.b(1) - real(vs))), vs = NaN; else, vs = real(vs); end
end

function psi = saddle_phase_gen(kern, q, c0, ts, s1, s2, R)
% the phase of the saddle factor of the smooth part: the exponent -t v^2 - s_d X(v) =
% -alpha v^2 - beta v - s_d/2, alpha = t + c0 s_d, beta = s_d q, on the first segment (s_d the
% root that governs the decay of phi_t there), has its saddle at v* = -beta/(2 alpha), where it
% exceeds its value at v = 0 by beta^2/(4 alpha).  The saddle is present when it lies ahead on
% the ray from 0 (no breakpoint: Re z < 0, z = beta/(2 sqrt(alpha))) or inside the first
% segment (0 < Re v* < b_1), and is not negligible (Re beta^2/(4 alpha) > -30).  With c0 < 0
% the saddle sits inside the segment at v* = -q s_d'/(2 (1 + c0 s_d')) (s_d' = ds_d/dt) for
% large |q|, with the phase growing like q^2 Im t, which the trapezoidal rule cannot follow.
if kern == 1 || R.sgn(1) < 0, sd = s1; else, sd = s2; end
alpha = ts + c0 * sd; beta = sd * q;
z2 = beta .* beta ./ (4 * alpha);
psi = nan(size(ts));
if R.ne == 0
    z = beta ./ (2 * sqrt(alpha));
    m = real(z) < 0 & real(z2) > -30;
else
    vs = -beta ./ (2 * alpha);
    m = real(vs) > 0 & real(vs) < R.b(1) & real(z2) > -30;
    if isfield(R, 'vsad') && isnan(R.vsad), m(:) = false; end
end
psi(m) = imag(z2(m));
end

function j = lastfirst(n)
% the order of the nodes of OUTER_GEN: the last first (the paths of the inner integrals grow with Im t, so that a point
% whose paths exceed the bounds of the engine raises the error at its last node, before the other nodes are evaluated;
% tmap_last! of DFExact.jl), then the others in order; the values do not depend on the order
j = [n, 1:n-1]; if n == 0, j = []; end
end

function val = outer_gen(model, c, q, c0, R, tau, inner, weight, o, stat)
% the outer integral for the general kernel: the smooth part as in OUTER, and one Filon part
% per breakpoint, with its frequency y_k = b_k^2
WMAX = getf(o, 'wmax', wmaxc(c));
hu = getf(o, 'hu', 0.075) * hufac(c); N = getf(o, 'N', 20); cde = tau;
NF = nodefac();
kern = kernel(model, stat);
umax = asinh(WMAX / cde);
us = (0:hu:umax)';
ts = tau + 1i * cde * sinh(us);
[P, s1, s2, M0] = tparts(model, c, ts, stat);
psi = saddle_phase_gen(kern, q, c0, ts, s1, s2, R);
if max(phase_steps(psi)) > 0.25
    w0 = min(0.5, tau / 2);
    edges0 = [0, w0 * 2.^(0:ceil(log2(WMAX / w0)))]; edges0(end) = WMAX;
    K = 8; wk = [];
    for p = 1:numel(edges0)-1
        tmp = linspace(edges0(p), edges0(p+1), K + 1); wk = [wk, tmp(1:end-1)]; %#ok<AGROW>
    end
    wk(end+1) = WMAX;
    te0 = tau + 1i * wk(:);
    [~, s1e, s2e] = tparts(model, c, te0, stat);
    steps = phase_steps(saddle_phase_gen(kern, q, c0, te0, s1e, s2e, R));
    wn = []; wt = [];
    for p = 1:numel(edges0)-1
        a = edges0(p); b = edges0(p+1);
        dp = sum(steps((p-1)*K + (1:K)));
        m = max(1, ceil(dp / 50));
        n = max(24, 16 * ceil(ceil(3 * dp / m) / 16)) * NF * getf(o, 'glfac', 1);
        [u, gw] = glrule(n);
        for k = 1:m
            aa = a + (b - a) * (k - 1) / m; bb = a + (b - a) * k / m;
            wn = [wn; aa + (bb - aa) * u]; wt = [wt; (bb - aa) * gw]; %#ok<AGROW>
        end
    end
    ts = tau + 1i * [0; wn];
    [P, s1, s2, M0] = tparts(model, c, ts, stat);
    gs = zeros(numel(wn), 1);
    for j = lastfirst(numel(wn))
        Is = inner(ts(j+1), s1(j+1), s2(j+1), M0(j+1), 1);
        gs(j) = real(P(j+1) * Is * weight(ts(j+1))) * wt(j);
    end
    val = sum(gs) / pi;
else
    gs = zeros(numel(ts), 1);
    for j = lastfirst(numel(ts))
        Is = inner(ts(j), s1(j), s2(j), M0(j), 1);
        gs(j) = real(P(j) * Is * weight(ts(j))) * cde * cosh(us(j));
    end
    gs(1) = gs(1) * 0.5;
    val = sum(gs) * hu / pi;
end
lastparts(val);
ehw = getf(o, 'edge_hw', []);
if ~isempty(ehw) && any(~R.negl)
    ws = (0:ehw:WMAX)'; te = tau + 1i * ws;
    [P, s1, s2, M0] = tparts(model, c, te, stat);
    ge = zeros(numel(te), 1);
    for j = 1:numel(te)
        [~, Ie] = inner(te(j), s1(j), s2(j), M0(j), 2);
        ge(j) = real(P(j) * sum(Ie) * weight(te(j)));
    end
    ge(1) = ge(1) * 0.5;
    val = val + sum(ge) * ehw / pi;
    return
end
pfac = getf(o, 'pfac', 1);
u = glrule(N);
for k = 1:R.ne
    if R.negl(k), continue; end
    y = R.y(k);
    w0 = min([0.5, 1 / y, tau / 2]);
    edges0 = [0, w0 * 2.^(0:ceil(log2(WMAX / w0)))]; edges0(end) = WMAX;
    edges = 0;
    for p = 1:numel(edges0)-1
        a = edges0(p); b = edges0(p+1);
        m = filonpanels(a, b, y, c);
        tmp = linspace(a, b, m + 1); edges = [edges, tmp(2:end)]; %#ok<AGROW>
    end
    if pfac > 1
        e2 = [];
        for p = 1:numel(edges)-1
            tmp = linspace(edges(p), edges(p+1), pfac + 1); e2 = [e2, tmp(1:end-1)]; %#ok<AGROW>
        end
        edges = [e2, edges(end)];
    end
    np = numel(edges) - 1;
    wn = reshape(edges(1:np) + (edges(2:end) - edges(1:np)) .* u, [], 1);
    te = tau + 1i * [0; wn];
    [P, s1, s2, M0] = tparts(model, c, te, stat);
    G = complex(zeros(numel(wn), 1));
    for j = lastfirst(numel(wn))
        t = te(j+1);
        [~, Ie] = inner(t, s1(j+1), s2(j+1), M0(j+1), 10 + k);
        G(j) = P(j+1) * Ie(k) * weight(t) * exp(t * y);
    end
    acc = 0;
    for p = 1:np
        acc = acc + filon_panel(G((p-1)*N+1:p*N), edges(p), edges(p+1), y, N);
    end
    val = val + real(exp(-tau * y) * acc) / pi;
end
end

% =====================================================================================
% Cancellation-free forms for cases DU, DT at large c (the DF-Test paper)

