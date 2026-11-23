function varargout = kan_core(action, varargin)
%KAN_CORE  The devices of Bao and Kan for the limits of the classical cases.
%
%   The numerical devices of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least
%   Squares Estimator in First-Order Autoregressions: O(1) Algorithms", ported to the limit pair (U, V) of these
%   programs (a port of src/DFKan.jl and src/DFKan_ratio.jl of DFLimit.jl and of c/dfl_kan.c; the method of every
%   function here is theirs):
%     - the joint characteristic function of the limit pair in closed forms that are free of cancellation on both
%       sides of the unit root, with the branch of the square root settled at each point from factors whose arguments
%       are known (no tracking along a grid), in Taylor-mode arithmetic (jets), so that derivatives of any order come
%       out of the same assembly;
%     - the single inversion for statistics that reject on a half-plane (the coefficient statistic, the linear form
%       U - aV, the power envelope): adaptive quadrature over decades and the tail completed by repeated integration
%       by parts, the derivatives read off a jet of order eight, each value with an error estimate;
%     - the moments from one jet per node of an adaptive quadrature in log s;
%     - the distribution function of the t-ratio at zero in closed form (the chi-square representation).
%   The pair is (U, V): the t-ratio is (U - c0 V)/sqrt(V) and the coefficient statistic (U - c0 V)/V.  Bao and Kan write
%   A = U - cV and B = V, so that E exp(i t1 A + i t2 B) = M(i t1, i (t2 - c t1)).  The quadrature is the algorithm of
%   QuadGK.jl (the 15-point Kronrod rule, bisection of the segment of the largest error estimate kept in a binary heap),
%   so that the values are those of the Julia package to rounding.  A jet of order n at m points is an n x m complex
%   matrix (row k + 1 the coefficient of x^k); a number at m points is a 1 x m row.
%
%   Actions (F the transform of a case, kan_core('limcf', ...)):
%     F = kan_core('limcf', model, stat, c, g0)       the limit pair of model 1, 2, 3 (cases DK, DU, DT; stat true: SK, SU,
%                                                      ST) at c with the fixed start g0
%     L = kan_core('logcf', F, t1, t2)                log E exp(i t1 A + i t2 B) at numbers or jets (n x m)
%     L = kan_core('logmgf', F, s, t)                 log E exp(sU + tV)
%     sc = kan_core('scale', F)                       [E A, E B, sd B, Var A] from a jet of order three at the origin
%     K = kan_core('kcf', F)                          the transform interface of the single inversion (lcf, scale, tail, w0)
%     [P, e] = kan_core('half_cdf', K, alpha, beta, b, tol, w0)     P[alpha A + beta B <= b] by the single inversion
%     [I, e] = kan_core('decades', f, b, tol, rtol)   int_0^b f over decades (f vectorized)
%     [I, e] = kan_core('quadgk', f, pts, atol, rtol, maxevals)     adaptive Gauss-Kronrod quadrature (QuadGK.jl)
%     x = kan_core('illinois', g, lo, hi, xtol, grow) the root of an increasing g by the Illinois variant of regula falsi
%     [P, e] = kan_core('coef_cdf1', F, y, c0, tol)   the coefficient statistic's distribution function at y
%     [f, e, ftol] = kan_core('coef_pdf1', F, y, c0, tol)   its density (tol relative to the scale of the density)
%     [y, e] = kan_core('coef_icdf1', F, p, c0, xtol, tol)  its p-quantile
%     [P, e] = kan_core('lin_cdf1', F, b, a, tol)     P(U - aV <= b)
%     [P, e] = kan_core('envelope1', c, alpha, tol)   the power envelope of Elliott, Rothenberg and Stock at c (case DK)
%     [m, e] = kan_core('ratio_mom', F, k, x0, power, tol)  m_j = E[(A - x0 B)^j / B^(j power)], j = 1..k
%     tf = kan_core('lin_ok', model, stat, c, g0, a), kan_core('mom_ok', c), kan_core('envelope_ok', c)
%                                                      the regions where the single inversion is the default
%     P = kan_core('ray_cdf', model, stat, c, g0)     P[U <= 0] in closed form (cases DK, DU, SK, SU)
%
%   See also COEF_CDF, COEF_PDF, COEF_ICDF, COEF_MOM, LTRATIO_LINCDF, POWER_ENVELOPE, LTRATIO_MOM, KAN_TEST.
switch lower(action)
    case 'limcf', varargout{1} = limcf(varargin{:});
    case 'logcf', F = varargin{1}; varargout{1} = lcf(F, varargin{2}, varargin{3});
    case 'logmgf'
        F = varargin{1}; s = varargin{2}; t = varargin{3};
        varargout{1} = lcf(F, -1i * s, -1i * (t + F.c * s));
    case 'scale', varargout{1} = scale(varargin{1});
    case 'kcf', varargout{1} = kcf(varargin{1});
    case 'half_cdf', [varargout{1:2}] = half_cdf(varargin{:});
    case 'decades', [varargout{1:2}] = decades(varargin{:});
    case 'quadgk', [varargout{1:2}] = qgk(varargin{:});
    case 'illinois', varargout{1} = illinois(varargin{:});
    case 'coef_cdf1', [varargout{1:2}] = coef_cdf1(varargin{:});
    case 'coef_pdf1', [varargout{1:3}] = coef_pdf1(varargin{:});
    case 'coef_icdf1', [varargout{1:2}] = coef_icdf1(varargin{:});
    case 'lin_cdf1', [varargout{1:2}] = lin_cdf1(varargin{:});
    case 'envelope1', [varargout{1:2}] = envelope1(varargin{:});
    case 'ratio_pdf', [varargout{1:3}] = ratio_pdf(varargin{:});
    case 'ratio_mom', [varargout{1:2}] = ratio_mom(varargin{:});
    case 'lin_ok', [model, stat, c, g0, a] = varargin{:}; varargout{1} = a ~= 0 && c >= -100 && c <= 2 + 2 * (model == 1 && ~stat && g0 == 0);
    case 'mom_ok', varargout{1} = varargin{1} <= 20;
    case 'envelope_ok', c = varargin{1}; varargout{1} = ~(c < 0) || c <= -1e-4;
    case 'ray_cdf', varargout{1} = ray_cdf(varargin{:});
    otherwise, error('kan_core: unknown action %s', action);
end
end

% =====================================================================================================================
% The series of the closed forms (the constants of src/DFKan.jl), computed once
% =====================================================================================================================
function S = ser()
persistent SS
if isempty(SS)
    FA = exp(gammaln((0:170)' + 1));                 % (as exp(loggamma(k + 1)) in the Julia package)
    f = @(k) FA(k + 1); KMS = 26; k = (0:KMS)';
    SS.fact = FA; SS.KMS = KMS;
    SS.cos = (-1).^k ./ f(2 * k); SS.sin = (-1).^k ./ f(2 * k + 1); SS.sh = (-1).^k ./ (4.^k .* f(2 * k + 1));
    k1 = (1:KMS + 1)'; SS.cd = (-1).^(k1 + 1) .* (2 * k1) ./ (4.^k1 .* f(2 * k1 + 1));
    SS.a = 2 * (-1).^k ./ f(2 * k + 2);
    KM = 40; kk = (0:KM)'; sn = (-1).^kk ./ f(2 * kk + 1); cs = (-1).^kk ./ f(2 * kk); a = zeros(KM + 1, 1);
    for i = 0:KM, acc = sn(i + 1); for j = 1:i, acc = acc - cs(j + 1) * a(i - j + 1); end, a(i + 1) = acc; end
    SS.d1 = a(2:end) / 4;
    j = (1:42)'; SS.f1 = 4 * j ./ (4.^j .* f(2 * j + 1));
    j = (2:43)'; SS.f2 = (12 ./ f(2 * j) - 4 ./ f(2 * j - 1) - 12 ./ f(2 * j + 1)) ./ 4.^j;
    k = (0:40)'; SS.ch = 1 ./ (4.^k .* f(2 * k)); SS.shc = 1 ./ (2.^(2 * k + 1) .* f(2 * k + 1));
end
S = SS;
end

% =====================================================================================================================
% Jets: n x m complex matrices, the coefficients of x^0..x^(n-1) at m points (Jet of src/DFKan.jl)
% =====================================================================================================================
function y = jconst(n, v), y = [v; complex(zeros(n - 1, numel(v)))]; end
function y = jvar(n, x0, lam), y = jconst(n, x0); if n > 1, y(2, :) = lam; end, end
function a = addr(a, r), a(1, :) = complex(real(a(1, :)) + r, imag(a(1, :))); end     % a real number added (the real part)
function a = subr(a, r), a(1, :) = complex(real(a(1, :)) - r, imag(a(1, :))); end
function a = rsub(r, a), v = complex(r - real(a(1, :)), -imag(a(1, :))); a = -a; a(1, :) = v; end
function m = jmul(a, b)
n = size(a, 1);
if n == 1, m = a .* b; return, end
m = complex(zeros(size(a)));
for k = 1:n, acc = a(1, :) .* b(k, :); for j = 2:k, acc = acc + a(j, :) .* b(k - j + 1, :); end, m(k, :) = acc; end
end
function m = jdiv(a, b)
n = size(a, 1);
if n == 1, m = a ./ b; return, end
m = complex(zeros(size(a))); ib = 1 ./ b(1, :);
for k = 1:n, acc = a(k, :); for j = 2:k, acc = acc - b(j, :) .* m(k - j + 1, :); end, m(k, :) = acc .* ib; end
end
function m = jsqrt(a)
n = size(a, 1); m = complex(zeros(size(a))); m(1, :) = sqrt(a(1, :));
if n == 1, return, end
h = 1 ./ (2 * m(1, :));
for k = 2:n, acc = a(k, :); for j = 2:k - 1, acc = acc - m(j, :) .* m(k - j + 1, :); end, m(k, :) = acc .* h; end
end
function m = jexp(a)
n = size(a, 1); m = complex(zeros(size(a))); m(1, :) = exp(a(1, :));
for k = 1:n - 1, acc = 0; for j = 1:k, acc = acc + (j * a(j + 1, :)) .* m(k - j + 1, :); end, m(k + 1, :) = acc / k; end
end
function m = jlogv(a, v)
% the logarithm of a jet with the constant coefficient v (the branch is the caller's)
n = size(a, 1); m = complex(zeros(size(a))); m(1, :) = v;
if n == 1, return, end
ia = 1 ./ a(1, :);
for k = 1:n - 1
    acc = k * a(k + 1, :);
    for i = 1:k - 1, acc = acc - (i * m(i + 1, :)) .* a(k - i + 1, :); end
    m(k + 1, :) = (acc .* ia) / k;
end
end
function m = jlogp(a), m = jlogv(a, log(a(1, :))); end
function m = jlogc(a, ph), m = jlogv(a, complex(log(abs(a(1, :))), ph)); end
function y = jser(cf, z)
% sum_k cf(k + 1) z^k by Horner's rule
y = complex(zeros(size(z)));
for k = numel(cf):-1:1, y = addr(jmul(y, z), cf(k)); end
end
function y = pick(mask, a, b), y = b; y(:, mask) = a(:, mask); end
function y = cols(a, mask), y = a(:, mask); end

% =====================================================================================================================
% The phases: the arguments of factors that are known (Bao and Kan)
% =====================================================================================================================
function r = rhe(x)
% round half to even (Julia's round)
r = round(x); h = abs(x - fix(x)) == 0.5; r(h) = 2 * round(x(h) / 2);
end
function a = argcos(z)
% the continuous argument of cos(theta) as a function of zeta = theta^2 (Re zeta <= 0): the principal argument shifted by
% the multiple of 2 pi that brings it within pi/2 of -sign(Im theta) Re theta
th = sqrt(z); tr = real(th); ti = imag(th);
e = exp(-2 * abs(ti));
cz = complex(cos(tr) .* (1 + e), -sin(tr) .* sign(ti) .* (1 - e));
a = angle(cz); target = -sign(ti) .* tr;
a = a + 2 * pi * rhe((target - a) / (2 * pi));
end
function a = argtau4(z0, zc0, D2, T0)
% the continuous phase of the factor of the demeaned cases under the square root as the sum of the principal arguments
% of three factors, each off the closed negative real axis
S = ser();
D1 = (2 * D2 - 1) ./ z0;
sm = abs(z0) < 4;
if any(sm), D1(sm) = jser(S.d1, z0(sm) / 4); end
rho = 1 + zc0 .* D1;
Sinf = T0 .* D2 ./ (2 * rho);
W = 1 ./ (2 * D2);
a = angle(W) + angle(rho) + angle(Sinf);
end

% =====================================================================================================================
% The joint characteristic function of the limit pair: the closed forms of Bao and Kan (their AR1LIM_LOGCF)
% =====================================================================================================================
function F = limcf(model, stat, c, g0)
if nargin < 4, g0 = 0; end
if ~any(model == [1 2 3]), error('kan_core: the model must be 1, 2 or 3 (cases DK, DU, DT, SK, SU, ST)'); end
if ~(isfinite(c) && isfinite(g0)), error('kan_core: c and g0 must be finite'); end
stat = logical(stat);
if stat && c >= 0, error('kan_core: a stationary start requires c < 0'); end
par = g0; if model == 2, par = c * g0; elseif model == 3, par = c * c * g0; end
if stat, par = 0; end
w0 = 0.5; if model == 1 && ~stat, w0 = (1 + par^2) / 2; end
ag = [];
if model == 3
    S = ser(); c2 = c * c; c3 = c2 * c; c4 = c3 * c; c5 = c4 * c;
    gam = [24 * c2 * (c - 2) * (c2 - 3 * c + 6), -24 * c2 * (c - 2) * (c2 - 3 * c + 6), 8 * (c5 - 4 * c4 + 6 * c3 + 6 * c2 - 36 * c + 36), ...
           -c5 + 24 * c3 - 120 * c2 + 288 * c - 288, c4 - 8 * c3 + 32 * c2 - 72 * c + 72];
    ag = zeros(41, 1);
    for k = 4:44, s = 0; for j = 0:4, s = s + 2 * gam(j + 1) / S.fact(2 * k - j + 1); end, ag(k - 3) = s; end
end
F = struct('det', model, 'stat', stat, 'c', c, 'par', par, 'w0', w0, 'ag', ag);
end

function nu = vtmc(vt, zc, c)
% vartheta - c without cancellation
if c >= 0, nu = jdiv(-zc, addr(vt, c)); else, nu = subr(vt, c); end
end
function nu = vtmo(vt, zc, c, om, t1, tau0)
% vartheta - omega without cancellation: on the explosive side as tau_0/(vartheta + omega) where that is the larger
nu = vtmc(vt, zc, c) - t1 * 1i;
if c > 0
    vpo = vt + om;
    sel = abs(vpo(1, :)) > abs(nu(1, :));
    if any(sel), nu = pick(sel, jdiv(tau0, vpo), nu); end
end
end
function y = bpoly(u, c, sg)
% B0(vartheta) of the trend cases expanded about vartheta = c (sg = 1) or -c (sg = -1), u = vartheta -/+ c
c2 = c * c; c3 = c2 * c; c4 = c3 * c; c5 = c4 * c; c6 = c5 * c; c7 = c6 * c;
k4 = c4 - 8 * c3 + 32 * c2 - 72 * c + 72;
if sg > 0
    q = c2 - 6 * c + 12;
    b0 = c * ((c - 2) * (c - 2)) * (q * q);
    b1 = 3 * c6 - 40 * c5 + 232 * c4 - 744 * c3 + 1344 * c2 - 1152 * c + 288;
    b2 = 3 * c5 - 32 * c4 + 152 * c3 - 408 * c2 + 576 * c - 288;
else
    b0 = -c7;
    b1 = 3 * c6 - 8 * c5 + 8 * c4 + 24 * c3 - 96 * c2 + 288;
    b2 = -3 * c5 + 16 * c4 - 40 * c3 + 24 * c2 + 144 * c - 288;
end
y = addr(jmul(addr(jmul(addr(u * k4, b2), u), b1), u), b0);
end

function L = lcf(F, t1, t2)
% log E exp(i t1 A + i t2 B), A = U - cV, B = V (_logcf of src/DFKan.jl): numbers (1 x m) or jets (n x m)
if size(t1, 1) == 1 && size(t2, 1) == 1                    % (numbers: a row of points, a scalar repeated)
    if isscalar(t1) && ~isscalar(t2), t1 = t1 * ones(size(t2)); end
    if isscalar(t2) && ~isscalar(t1), t2 = t2 * ones(size(t1)); end
end
t1 = complex(t1); t2 = complex(t2);
c = F.c; c2 = c * c; n = size(t1, 1); S = ser();
om = addr(t1 * 1i, c);                                       % omega = c + i t1
zc = t2 * 2i - t1 * complex(0 * c, 2 * c);                   % theta^2 + c^2
zeta = subr(zc, c2);                                         % theta^2
tau0 = jmul(t1, t1) - t2 * 2i;
if F.det == 3, L = trendlim(F, t1, t2, om, zc, tau0); return, end
z0 = zeta(1, :); small = abs(z0) < 4;
if any(small) && ~all(small)                                 % (the points of each branch apart)
    L = complex(zeros(size(t1)));
    L(:, small) = lcf(F, cols(t1, small), cols(t2, small)); L(:, ~small) = lcf(F, cols(t1, ~small), cols(t2, ~small));
    return
end
small = all(small);
if F.det == 1
    if small
        C = jser(S.cos, zeta); Sn = jser(S.sin, zeta);
        f = C - jmul(om, Sn);
        g = rsub(1, jdiv(jmul(om, Sn), C));
        lC = jlogc(C, argcos(z0)); lg = jlogp(g);
        Kf = jdiv(-jmul(tau0, Sn), f * 2);
    else
        vt = jsqrt(-zeta); E2 = jexp(vt * -2);
        nu = vtmo(vt, zc, c, om, t1, tau0);
        Bk = nu + jmul(E2, vt + om);
        oe = addr(E2, 1);
        lC = subr(vt + jlogp(oe), log(2));                   % log cos theta: Im vartheta + Arg(1 + e^{-2 vartheta})
        lg = jlogp(jdiv(Bk, jmul(vt, oe)));
        Kf = jdiv(-jmul(tau0, rsub(1, E2)), Bk * 2);
    end
    L = (-om) / 2 - (lC + lg) / 2;
    if F.stat, L = L - jlogp(addr(Kf / c, 1)) / 2;
    elseif F.par ~= 0, L = L + Kf * (F.par * F.par); end
    return
end
% an intercept
tau1 = (jmul(t1, t1) + t1 * complex(0, c2) - t2 * complex(0 * c, 2 * c)) * 2;
tau2 = addr(t1 * (c * complex(0, c - 2)) - tau0, c2 * c);
if small
    Sn = jser(S.sin, zeta); Sh = jser(S.sh, zeta); ShmCh = jser(S.cd, zeta);
    tau10 = tau1 + zeta * (2 * c); tau20 = tau2 - zeta;
    Nz = complex(zeros(size(t1)));
    for k = S.KMS:-1:1                                       % N/zeta, N = 2 tau1 Sh^2 - 4 c^2 cos th + 4 tau2 sin th/th
        cst = -4 * c * S.a(k) - 4 * c2 * S.cos(k + 1) + 4 * S.sin(k);
        Nz = jmul(Nz, zeta) + addr(tau10 * (2 * S.a(k + 1)) + tau20 * (4 * S.sin(k + 1)), cst);
    end
    T0 = jdiv(Nz, jmul(Sh, Sh));
    lS = jlogp(Sh);
    Q = jdiv(jmul(jmul(tau0, ShmCh), Sh) * 2, Nz);
    D2 = Sh(1, :).^2 ./ (2 * Sn(1, :));
else
    vt = jsqrt(-zeta); E = jexp(-vt); E2 = jmul(E, E);
    nu = vtmo(vt, zc, c, om, t1, tau0);
    Q0 = jmul(t1, t1) * ((2 - c) * (1 + c)) + t2 * ((c - 1) * complex(0 * c, 2 * c)) + jmul(nu, tau2) - jmul(t1, tau0) * 1i;
    vt2 = jmul(vt, vt);
    Den = Q0 - jmul(E, tau1) * 2 + jmul(E2, tau1 - vt2 * c2 - jmul(tau2, vt));
    omE = rsub(1, E);
    T0 = jdiv(Den * -2, jmul(vt2, jmul(omE, omE)));
    lS = vt / 2 + jlogp(omE) - jlogp(vt);                    % sin(th/2)/(th/2) = e^{vt/2}(1 - E)/vt
    br = rsub(1, vt / 2) - jmul(E, addr(vt / 2, 1));
    Q = jdiv(jmul(jmul(tau0, br), omE), Den);
    D2 = omE(1, :) ./ ((1 + E(1, :)) .* vt(1, :));
end
lT = jlogc(T0, argtau4(z0, zc(1, :), D2, T0(1, :)));
L = rsub(log(2), om / 2) - lS - lT / 2;
if F.stat, L = L - jlogp(rsub(1, Q * c)) / 2;
elseif F.par ~= 0, L = L - Q * (F.par * F.par); end
end

function L = trendlim(F, t1, t2, om, zc, tau0)
% the cases with an intercept and a trend (the trend limit of Bao and Kan): the phase is that of the case without
% deterministic terms corrected by the principal arguments of three factors that are off the closed negative real axis
c = F.c; c2 = c * c; S = ser();
sm = abs(zc(1, :) - c2) < 16;
if any(sm) && ~all(sm)
    L = complex(zeros(size(t1)));
    L(:, sm) = trendlim(F, cols(t1, sm), cols(t2, sm), cols(om, sm), cols(zc, sm), cols(tau0, sm));
    L(:, ~sm) = trendlim(F, cols(t1, ~sm), cols(t2, ~sm), cols(om, ~sm), cols(zc, ~sm), cols(tau0, ~sm));
    return
end
Kb = tau0 * (-4 * (c2 - 3 * c + 3)) + (t1 * 1i) * (c * (c - 2) * (c2 - 6 * c + 12));
zv = rsub(c2, zc);                                           % vartheta^2
if all(sm)
    z0 = zv(1, :);
    f1 = jser(S.f1, zv); f2 = jser(S.f2, zv); Ge = jser(F.ag, zv);
    f12 = jmul(f1, f2);
    Psi = jmul(f12, Kb) + Ge;                                % e^{vartheta} D3/vartheta^8
    Lr = addr((-om) / 2 - jlogp(Psi) / 2, log(2) / 2);
    Q = jdiv(jmul(tau0, f12), Psi * 2);
    g11 = f1(1, :) ./ (2 * jser(S.ch, z0)); gll = -f2(1, :) ./ (24 * jser(S.shc, z0));
else
    vt = jsqrt(zv); E = jexp(-vt); E2 = jmul(E, E);
    if c >= 0
        vmc = jdiv(-zc, addr(vt, c)); vpc = addr(vt, c);
        Bp = bpoly(vmc, c, 1); Bm = bpoly(-vmc, c, -1);
    else
        vpc = jdiv(-zc, subr(vt, c)); vmc = subr(vt, c);
        Bp = bpoly(vpc, c, -1); Bm = bpoly(-vpc, c, 1);
    end
    B1 = subr(zv * ((c - 3) * (c2 - 4 * c + 6)), 6 * (c2 * c - 5 * c2 + 12 * c - 12));
    h1 = subr(vt, 2) + jmul(E, addr(vt, 2));
    h2 = -addr(zv - vt * 6, 12) + jmul(E, addr(zv + vt * 6, 12));
    h12 = jmul(h1, h2);
    if c > 0                                                 % the explosive side: every term of the order of D3 itself
        K1 = c * (c - 2) * (c2 - 6 * c + 12);
        Pp = addr(jmul(vt, addr(zv + vt * 8, 24)), 24);
        Qc = subr(vt * (c - 3), 3 * (c - 2));
        vQ = jmul(vmc, Qc);
        ta = jmul(vQ, vQ) * 8;
        tb = jmul(vmc, jmul(subr(zv, 6), E) * 8 + jmul(Pp, E2)) * K1;
        tc = jmul(jmul(tau0, h12), addr(jdiv(jconst(size(vt, 1), K1 * ones(1, size(vt, 2))), vt + om), 4 * (c2 - 3 * c + 3)));
        td = jmul(zc, jmul(B1, E)) * 8;
        te = jmul(vpc, jmul(Bm, E2));
        D3 = (((ta + tb) - tc) + td) - te;
    else
        D3 = jmul(h12, Kb) + (jmul(vmc, Bp) + jmul(zc, jmul(B1, E)) * 8 - jmul(vpc, jmul(Bm, E2)));
    end
    Lr = addr(-(om + vt) / 2 - (jlogp(D3) - jlogp(vt) * 8) / 2, log(2) / 2);
    Q = jdiv(jmul(tau0, h12), D3 * 2);
    v0 = vt(1, :); e0 = E(1, :);
    g11 = h1(1, :) ./ (v0.^3 .* (1 + e0)); gll = -h2(1, :) ./ (12 * v0.^4 .* (1 - e0));
end
% the continuous phase from the factors
t10 = t1(1, :); t20 = t2(1, :);
FK = F; FK.det = 1; FK.stat = false; FK.par = 0; FK.w0 = 0.5;
L0 = lcf(FK, t10, t20); z1 = zc(1, :);
tr = z1 .* (g11 + 12 * gll);
dt = (12 * z1.^2 + 144 * t10.^2) .* g11 .* gll;
sq = sqrt(tr.^2 - 4 * dt);
k1 = (tr + sq) / 2; k2 = (tr - sq) / 2;
rho = exp(2 * (L0 - Lr(1, :))) ./ (1 + tr + dt);
ph = imag(L0) - (angle(1 + k1) + angle(1 + k2) + angle(rho)) / 2;
L = Lr; L(1, :) = complex(real(Lr(1, :)), ph);
if F.stat, L = L - jlogp(addr(Q * (c2 * c), 1)) / 2;
elseif F.par ~= 0, L = L + Q * (F.par * F.par); end
end

% =====================================================================================================================
% Adaptive Gauss-Kronrod quadrature (QuadGK.jl 2.11, order 7)
% =====================================================================================================================
function [I, E] = qgk(f, pts, atol, rtol, maxevals)
% f vectorized (a row of points to a row of values); the segments between the points pts
persistent X W G
if isempty(X)
    X = [-9.9145537112081263920685469752598e-01, -9.4910791234275852452618968404809e-01, -8.6486442335976907278971278864098e-01, ...
         -7.415311855993944398638647732811e-01, -5.8608723546769113029414483825842e-01, -4.0584515137739716690660641207707e-01, ...
         -2.0778495500789846760068940377309e-01, 0.0];
    W = [2.2935322010529224963732008059913e-02, 6.3092092629978553290700663189093e-02, 1.0479001032225018383987632254189e-01, ...
         1.4065325971552591874518959051021e-01, 1.6900472663926790282658342659795e-01, 1.9035057806478540991325640242055e-01, ...
         2.0443294007529889241416199923466e-01, 2.0948214108472782801299917489173e-01];
    G = [1.2948496616886969327061143267787e-01, 2.797053914892766679014677714229e-01, 3.8183005050511894495036977548818e-01, ...
         4.1795918367346938775510204081658e-01];
end
ns = numel(pts) - 1;
H = zeros(ns, 4);                                            % the heap: [a b I E] by rows
for i = 1:ns, H(i, :) = evalrule(f, pts(i), pts(i + 1), X, W, G); end
numevals = 15 * ns;
I = H(1, 3); E = H(1, 4); for i = 2:ns, I = I + H(i, 3); E = E + H(i, 4); end
if numevals >= maxevals || E <= atol || E <= rtol * abs(I), return, end
len = ns;
for i = floor(len / 2):-1:1, H = hdown(H, i, H(i, :), len); end
while E > atol && E > rtol * abs(I) && numevals < maxevals
    s = H(1, :); y = H(len, :); len = len - 1;
    if len > 0, H = hdown(H, 1, y, len); end
    mid = (s(1) + s(2)) / 2;
    if roundoff(s(1), mid, X) || roundoff(mid, s(2), X), len = len + 1; H = hup(H, len, s); break, end
    s1 = evalrule(f, s(1), mid, X, W, G); s2 = evalrule(f, mid, s(2), X, W, G);
    I = (I - s(3)) + s1(3) + s2(3); E = (E - s(4)) + s1(4) + s2(4); numevals = numevals + 30;
    len = len + 1; H = hup(H, len, s1);
    len = len + 1; H = hup(H, len, s2);
end
I = H(1, 3); E = H(1, 4); for i = 2:len, I = I + H(i, 3); E = E + H(i, 4); end
end
function s = evalrule(f, a, b, X, W, G)
h = 0.5 * (b - a);
x = [a + h, a + (1 + X(7)) * h, a + (1 - X(7)) * h, a + (1 + X([2 4 6])) * h, a + (1 - X([2 4 6])) * h, a + (1 + X([1 3 5])) * h, a + (1 - X([1 3 5])) * h];
v = f(x);
f0 = v(1);
Ig = f0 * G(4);
Ik = f0 * W(8) + (v(2) + v(3)) * W(7);
for i = 1:3
    fg = v(3 + i) + v(6 + i); fk = v(9 + i) + v(12 + i);
    Ig = Ig + fg * G(i);
    Ik = Ik + (fg * W(2 * i) + fk * W(2 * i - 1));
end
Iks = Ik * h; Igs = Ig * h; E = abs(Iks - Igs);
if ~isfinite(E), error('kan_core:quadgk', 'the quadrature of the single inversion produced %g on (%g, %g)', E, a, b); end
s = [a b Iks E];
end
function H = hdown(H, i, x, len)
% the binary heap of the segments by their error estimate, largest first (DataStructures.jl with the reverse ordering)
while 2 * i <= len
    l = 2 * i; r = l + 1;
    if r > len || H(r, 4) < H(l, 4), j = l; else, j = r; end
    if ~(x(4) < H(j, 4)), break, end
    H(i, :) = H(j, :); i = j;
end
H(i, :) = x;
end
function H = hup(H, i, x)
while i > 1
    j = floor(i / 2);
    if ~(H(j, 4) < x(4)), break, end
    H(i, :) = H(j, :); i = j;
end
H(i, :) = x;
end
function tf = roundoff(a, b, X), c = 0.5 * (b - a); tf = a == a + (1 + X(1)) * c || b == a + (1 - X(1)) * c; end

% =====================================================================================================================
% The single inversion with tail completion (half_cdf of DFKan)
% =====================================================================================================================
function [I, e] = decades(f, b, tol, rtol)
% int_0^b f by adaptive quadrature over [0, 10], [10, 100], ..., the last up to b; with b = Inf it stops after two
% consecutive negligible decades
I = 0; e = 0; a = 0; small = 0; k = 1;
while a < b
    if isinf(b), hi = 10^k; elseif 10^(k + 0.5) >= b, hi = b; else, hi = 10^k; end
    [Ik, ek] = qgk(f, [a hi], tol / 10, rtol, 1e6);
    I = I + Ik; e = e + ek; a = hi; k = k + 1;
    if abs(Ik) + ek <= tol / 20, small = small + 1; else, small = 0; end
    if small >= 2 || k > 31, break, end
end
end
function sc = ray_scale(s4, alpha, beta)
vv = alpha^2 * s4(4) + beta^2 * s4(3)^2;
if ~(isfinite(vv) && vv > 0), vv = s4(4); end
sc = 1 / sqrt(vv);
end
function [P, err] = half_cdf(K, alpha, beta, b, tol, w0)
% P[alpha A + beta B <= b] (alpha ~= 0): 1/2 - (1/pi) int_0^inf Im[phi(u alpha, u beta) e^{-i b u}]/u du, over [0, S] by
% adaptive quadrature over decades and beyond S by the tail's expansion (S = 400/|w|, w = w0 alpha + b; doubled, at most
% five times, while the first omitted term exceeds tol/10)
if nargin < 6, w0 = K.w0; end
sc = ray_scale(K.scale(), alpha, beta);
f = @(z) hcf(K, alpha, beta, b, sc, z);
w = w0 * alpha + b; T = 0; et = 0; S = Inf;
if abs(w) > 1e-7
    S = 400 / abs(w);
    for dbl = 0:5
        [Tc, et] = K.tail(alpha, beta, b, w, S, true); T = imag(Tc);
        if et <= tol / 10 || dbl == 5, break, end
        S = S * 2;
    end
end
[I, e] = decades(f, S / sc, tol, max(tol, 1e-14));
P = 0.5 - (I + T) / pi; err = (e + et) / pi;
end
function y = hcf(K, alpha, beta, b, sc, z)
u = z * sc;
L = K.lcf(alpha * u, beta * u);
y = exp(real(L)) .* sin(imag(L) - b * u) ./ z;
y(~isfinite(y)) = 0;
end
function K = kcf(F)
% the transform interface of the classical case F
K = struct('c', F.c, 'w0', F.w0, 'lcf', @(t1, t2) lcf(F, t1, t2), 'scale', @() scale(F), 'tail', @(alpha, beta, b, w, S, cdf) lim_tail(F, alpha, beta, b, w, S));
end
function s4 = scale(F)
% (E A, E B, sd B, Var A) from jets of order three at the origin (ratio_scale)
z = jconst(3, 0); x = jvar(3, 0, 1);
Lt = lcf(F, z, x); Ls = lcf(F, x, z);
Ev = real(Lt(2) / 1i); Vv = real(-2 * Lt(3));
Eu = real(Ls(2) / 1i); Vu = real(-2 * Ls(3));
sdv = Inf; if isfinite(Vv), sdv = sqrt(max(Vv, 0)); end
s4 = [Eu, Ev, sdv, Vu];
end
function [T, e] = lim_tail(F, alpha, beta, b, w, S)
% the tail int_S^inf Im[phi(u alpha, u beta) e^{-i b u}]/u du by repeated integration by parts, the derivatives of the
% slowly varying factor from a jet of order eight at S (tail_cdf of DFKan)
N = 9; FA = ser(); fa = FA.fact;
sj = jvar(N, S, 1);
ph = jexp(lcf(F, sj * alpha, sj * beta) - sj * complex(0, b));
g = jdiv(jmul(ph, jexp(jvar(N, 0, complex(0, w)))), sj);
acc = 0; iw = 1 / complex(0, w); p = iw;
for k = 0:N - 2, acc = acc + (fa(k + 1) * g(k + 1)) * p; p = p * iw; end
T = acc; e = abs((fa(N) * g(N)) * p);
end

% ---- the density of a ratio (ratio_pdf of DFKan) ----
function [f, err, fs] = ratio_pdf(F, r, tol)
s4 = scale(F); Ev = s4(2); sdv = s4(3); Vu = s4(4);
fs = max(1, Ev / sqrt(max(Vu, realmin)));
vv = Vu + r^2 * sdv^2; if ~(isfinite(vv) && vv > 0), vv = Vu; end
sc = 1 / sqrt(vv);
fn = @(z) rpf(F, r, sc, z);
w = F.w0; S = 400 / w; T = 0; et = 0;
for dbl = 0:5
    [T, et] = tail_pdf(F, r, w, S, 1 / Ev);
    if et <= tol / 10 || dbl == 5, break, end
    S = S * 2;
end
[I, e] = decades(fn, S / sc, tol * fs, max(tol, 1e-14));
f = (I + T) / pi; err = (e + et) / pi;
end
function y = rpf(F, r, sc, z)
u = z * sc; m = numel(u);
L = lcf(F, [u; zeros(1, m)], [-r * u; ones(1, m)]);
y = sc * imag(exp(L(1, :)) .* L(2, :));
y(~isfinite(y)) = 0;
end
function [T, lst] = tail_pdf(F, r, w, S, rho)
% the tail int_S^inf Im H(u) du of the density, H(u) = d/dt2 phi(t1, t2) at (u, -r u): the derivatives of H along the ray
% are mixed derivatives of phi, read off jets of x -> phi(S + x, -r (S + x) + xi x) at N complex xi on a circle (the
% coefficient of xi in the k-th derivative is k H^(k-1)(S), extracted by a discrete Fourier transform), as in RATIO_PDF of
% Bao and Kan
N = 9; K = N - 1; FA = ser(); fa = FA.fact;
a1 = complex(zeros(N, 1));
x = jvar(N, 0, 1); sj = jvar(N, S, 1);
for m = 0:K
    ang = 2 * pi * m / N; xi = rho * complex(cos(ang), sin(ang));
    ph = jexp(lcf(F, sj, x * xi - sj * r));
    for k = 0:K, a1(k + 1) = a1(k + 1) + (fa(k + 1) * ph(k + 1)) / xi; end
end
a1 = a1 / N;
H = a1(2:K + 1) ./ (1:K)';                                   % H^(0..K-1)(S)
iwp = complex(0, w).^(0:N)';
T = 0; lst = 0;
for k = 0:K - 1
    acc = 0; bin = 1;
    for l = 0:k, acc = acc + (bin * iwp(k - l + 1)) * H(l + 1); bin = bin * (k - l) / (l + 1); end
    term = acc / iwp(k + 2);
    if k < K - 1, T = T + term; else, lst = abs(term); end
end
T = imag(T);
end

% ---- the root of an increasing function (illinois of DFKan) ----
function x = illinois(g, lo, hi, xtol, grow)
if nargin < 5 || isempty(grow) || ~(grow > 0), grow = hi - lo; end
glo = g(lo); ghi = g(hi); d = grow; k = 0;
while glo > 0 && k < 60, hi = lo; ghi = glo; lo = lo - d; d = d * 2; glo = g(lo); k = k + 1; end
while ghi < 0 && k < 120, lo = hi; glo = ghi; hi = hi + d; d = d * 2; ghi = g(hi); k = k + 1; end
if ~(isfinite(glo) && isfinite(ghi) && glo <= 0 && ghi >= 0)
    error('kan_core:bracket', 'the quantile could not be bracketed (the distribution function is not finite or not monotone)');
end
side = 0; x = (lo + hi) / 2;
for it = 1:200
    x = (lo * ghi - hi * glo) / (ghi - glo);
    if ~(lo < x && x < hi), x = (lo + hi) / 2; end
    gx = g(x);
    if gx == 0, break, end
    if gx < 0, lo = x; glo = gx; if side == -1, ghi = ghi / 2; end, side = -1;
    else, hi = x; ghi = gx; if side == 1, glo = glo / 2; end, side = 1; end
    if hi - lo < xtol, x = (lo * ghi - hi * glo) / (ghi - glo); break, end
end
end

% ---- the moments from jets (ratio_mom, mom_one of DFKan) ----
function [m, err] = ratio_mom(F, k, x0, power, tol)
if k > 31, error('kan_core: the moments from jets go up to p = 31'); end
s4 = scale(F); Ev = s4(2); FA = ser();
pts = [-60, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 45];
m = zeros(1, k); err = zeros(1, k);
for j = 1:k
    a = j * power; cst = FA.fact(j + 1) / gamma(a);
    [m(j), err(j)] = qgk(@(y) mof(F, j + 1, x0, a, cst, Ev, y), pts, tol * 1e-2, tol, 1e6);
end
end
function v = mof(F, N, x0, a, cst, Ev, y)
j = N - 1; s2 = exp(y) / Ev; mm = numel(y);
t1 = jvar(N, zeros(1, mm), complex(0, -1));
L = lcf(F, t1, jvar(N, complex(0, s2), complex(0, x0) * ones(1, mm)));
l0 = L(1, :); L(1, :) = 0;
Ex = jexp(L);
v = cst * real(exp(complex(a * log(s2) + real(l0), imag(l0))) .* Ex(j + 1, :));
v(~isfinite(v)) = 0;
end

% =====================================================================================================================
% The statistics of the classical cases by the single inversion (coef_cdf1 and so on of DFKan)
% =====================================================================================================================
function [P, e] = coef_cdf1(F, y, c0, tol)
if isnan(y), P = NaN; e = 0; return, end
if isinf(y), P = double(y > 0); e = 0; return, end
[P, e] = half_cdf(kcf(F), 1, -(y + c0 - F.c), 0, tol, F.w0);
P = min(max(P, 0), 1);
end
function [f, e, ftol] = coef_pdf1(F, y, c0, tol)
if ~isfinite(y), f = 0; if isnan(y), f = NaN; end, e = 0; ftol = tol; return, end
[f, e, fs] = ratio_pdf(F, y + c0 - F.c, tol);
f = max(f, 0); ftol = tol * fs;
end
function [y, e] = coef_icdf1(F, p, c0, xtol, tol)
s4 = scale(F); Eu = s4(1); Ev = s4(2); Vu = s4(4);
sd0 = sqrt(Vu) / Ev + abs(Eu) / Ev; m0 = Eu / Ev;
K = kcf(F); el = 0;
    function v = g(x)
        [P, ee] = half_cdf(K, 1, -x, 0, tol, F.w0); el = ee; v = P - p;
    end
r = illinois(@g, m0 - 2 * sd0, m0 + 2 * sd0, xtol);
y = r - (c0 - F.c); e = el;
end
function [P, e] = lin_cdf1(F, b, a, tol)
if isnan(b), P = NaN; e = 0; return, end
if isinf(b), P = double(b > 0); e = 0; return, end
[P, e] = half_cdf(kcf(F), 1, F.c - a, b, tol, F.w0);
P = min(max(P, 0), 1);
end
function [P, e] = envelope1(c, alpha, tol)
% the asymptotic power envelope of Elliott, Rothenberg and Stock at c < 0 without deterministic terms (and with an
% intercept, the same curve): the power at c of the test that rejects for small U - (c/2) V, its critical value the
% alpha-quantile at the unit root by the Illinois variant of regula falsi to 1e-12 from [-0.45, 0.5] (grown in steps of
% one), both by the single inversion; e the sum of the two estimates at the end
if isnan(c), P = NaN; e = 0; return, end
if ~(c < 0), P = alpha; e = 0; return, end
if c == -Inf, P = 1; e = 0; return, end
a = c / 2;
K0 = kcf(limcf(1, false, 0, 0)); Kc = kcf(limcf(1, false, c, 0)); e0 = 0;
    function v = g(x)
        [Pp, ee] = half_cdf(K0, 1, -a, x, tol, K0.w0); e0 = ee; v = Pp - alpha;
    end
b = illinois(@g, -0.45, 0.5, 1e-12, 1);
[P, e1] = half_cdf(Kc, 1, c - a, b, tol, Kc.w0);
P = min(max(P, 0), 1); e = e1 + e0;
end

% =====================================================================================================================
% The distribution function of the t-ratio at zero in closed form (AR1LIM_RAY of Bao and Kan; ray_cdf of DFKan)
% =====================================================================================================================
function v = Phi_(x), v = erfc(-x / sqrt(2)) / 2; end
function v = npdf_(x), v = exp(-x .* x / 2) / sqrt(2 * pi); end
function v = cexp1(x), if x == 0, v = 1; else, v = expm1(x) / x; end, end
function P = ray_cdf(model, stat, c, g0)
if nargin < 4, g0 = 0; end
if ~stat && model == 1
    c0 = cexp1(2 * c); dl = g0^2 * exp(2 * c) / c0; a = sqrt((1 + g0^2) / c0);
    P = Phi_(a - sqrt(dl)) - Phi_(-a - sqrt(dl)); return
end
if ~stat && model == 2, P = ray_u(c, c * g0); return, end
if ~(c < 0), error('kan_core: ray_cdf: a stationary start requires c < 0'); end
if model == 1
    lc = sqrt(-expm1(2 * c)) / (-4 * c); a = 1 / (2 * lc);
    P = qgk(@(z) erf(sqrt(max(a + z .* z, 0) / 2)) .* npdf_(z), [-12 0 12], 1e-14, 1e-12, 1e7); return
end
if model == 2
    sg = sqrt(-c / 2);
    P = qgk(@(z) arrayfun(@(zz) ray_u(c, sg * zz), z) .* npdf_(z), [-12 0 12], 1e-12, 1e-10, 1e7); return
end
error('kan_core: ray_cdf: cases DK, DU, SK, SU');
end
function P = ray_u(c, kap)
% an intercept: U = lambda_1 X_1 - lambda_2 X_2 - 1/2 with independent noncentral chi-square(1) variables
C5 = [1/12, 1/12, 17/360, 7/360, 43/6720, 107/60480, 769/1814400, 163/1814400, 4097/239500800, 709/239500800, ...
      6827/14529715200, 15019/217945728000, 19661/2092278988800, 1139/951035904000];
c0 = cexp1(2 * c); c1 = cexp1(c);
if abs(c) < 0.25, q = 0; for k = 14:-1:1, q = q * c + C5(k); end, c2 = c^2 * q;
else, c2 = c0 - c1^2; q = c2 / c^2; end
vs = sqrt(c2^2 + 4 * q);
l1 = (vs + c2) / 4; l2 = q / (vs + c2);
if c == 0, d1 = kap^2 / 2; else, d1 = kap^2 / 2 * (1 - (c - 4) * c2 / (c * vs)); end
d2 = kap^2 - d1;
sd1 = sqrt(max(d1, 0)); sd2 = sqrt(max(d2, 0));
wp = [-12, 12];
if sd2 < 12, wp(end + 1) = -sd2; end
ys = (l1 * d1 - 0.5) / l2;
if ys > 0, for sg = [-1 1], zz = sg * sqrt(ys) - sd2; if abs(zz) < 12, wp(end + 1) = zz; end, end, end
wp = unique(wp);
F1 = @(y) Phi_(sqrt(y) - sd1) - Phi_(-sqrt(y) - sd1);
P = qgk(@(z) F1((0.5 + l2 * (sd2 + z).^2) / l1) .* npdf_(z), wp, 1e-14, 1e-12, 1e7);
end
