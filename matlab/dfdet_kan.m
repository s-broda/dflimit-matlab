function varargout = dfdet_kan(action, varargin)
%DFDET_KAN  The devices of Bao and Kan for the general engine: the transform with the branch from the pivots.
%
%   The joint transform of the pair (U, V) of a model of the general engine (a design with a zero or stationary start,
%   any detrending) at c, with the branch of the square root of its determinant fixed at each point from the pivots of
%   the precision matrix, and the single inversion of the point-optimal statistic (a port of src/DFKan_det.jl of
%   DFLimit.jl and of c/dfl_det_kan.c).  The method is that of Yong Bao and Raymond Kan, "Exact Finite-Sample
%   Distribution and Moments of the Least Squares Estimator in First-Order Autoregressions: O(1) Algorithms": det(I -
%   2 Q Sigma) = det(K) det(Sigma) with K = Sigma^{-1} - 2 Q the precision matrix of the functionals under the tilt;
%   both have a positive definite Hermitian part for imaginary (s, t), for any design, start and detrending, so that
%   elimination without pivoting has every pivot in the open right half-plane and the logarithm of the determinant, with
%   its continuous argument, is the sum of the principal logarithms of the pivots (with a stationary start det Sigma_x =
%   v det Sigma, v = -1/(2c) > 0).  The fixed grid of the inversion of Gil-Pelaez (DFDET_QGRID) continued the phase of
%   the determinant along u instead.
%
%   K = DFDET_KAN('cf', Mo, c)          the transform of the model Mo (DFDET_STARTMODEL, DFDET_MODEL) at c, the interface
%                                       of the single inversion of KAN_CORE (fields lcf, scale, tail, w0, c)
%   K = DFDET_KAN('cf', Mo, c, true)    the same with the values along the ray kept (the quantile searches)
%   L = DFDET_KAN('logcf', K, t1, t2)   log E exp(i t1 A + i t2 B), A = U - cV, B = V, at real points (rows)
%   [P, e] = DFDET_KAN('pt_cdf1', K, cbar, q, tol)    P(Q <= q), Q = U + 1/2 + cbar^2 V = A + (cbar^2 + c) B + 1/2, by the
%                                       single inversion with the tail completed by integration by parts
%   [f, e] = DFDET_KAN('pt_pdf1', K, cbar, q, tol)    its density
%   [q, e] = DFDET_KAN('pt_quantile1', K, cbar, p, tol)  its p-quantile: regula falsi from E Q -/+ its standard deviation
%                                       to 1e-12 times it, with the values along the ray kept across the search
%
%   See also PT_CDF, PT_PDF, PT_QUANTILE, PT_POWER, PT_ENVELOPE, KAN_CORE.
switch lower(action)
    case 'cf'
        memo = false; if numel(varargin) > 2, memo = varargin{3}; end
        varargout{1} = mkcf(varargin{1}, varargin{2}, memo);
    case 'logcf', varargout{1} = lcfv(varargin{:});
    case 'pt_cdf1', [varargout{1:2}] = pt_cdf1(varargin{:});
    case 'pt_pdf1', [varargout{1:2}] = pt_pdf1(varargin{:});
    case 'pt_quantile1', [varargout{1:2}] = pt_quantile1(varargin{:});
    otherwise, error('dfdet_kan: unknown action %s', action);
end
end

function K = mkcf(Mo, c, memo)
if strcmp(Mo.start, 'fixed') && Mo.g0 ~= 0, error('dfdet:args', 'dfdet_kan: a zero or stationary start (the fixed start''s transform has a further term)'); end
K = struct('Mo', Mo, 'c', c, 'w0', Mo.a, 'memo', []);
if memo, K.memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
s4 = scale(K);
K.lcf = @(t1, t2) lcfv(K, t1, t2);
K.scale = @() s4;
K.tail = @(alpha, beta, b, w, S, cdf) tail_int(K, alpha, beta, b, w, S, cdf);
end

function L = lcfv(K, t1, t2)
% the transform at the points (t1(k), t2(k)), from the memo where it is kept
L = complex(zeros(size(t1)));
for k = 1:numel(t1)
    if isempty(K.memo), L(k) = raw(K, t1(k), t2(k)); continue, end
    key = [num2hex(t1(k)) num2hex(t2(k))];
    if isKey(K.memo, key), L(k) = K.memo(key); else, L(k) = raw(K, t1(k), t2(k)); K.memo(key) = L(k); end
end
end

function L = raw(K, t1, t2)
% log M(i t1, i (t2 - c t1)) at one point with the logarithm of the determinant from the pivots
Mo = K.Mo; P = Mo.P; c = K.c;
s = 1i * complex(t1); t = 1i * (complex(t2) - c * complex(t1));
theta = sqrt(complex(c * c - 2 * t)); gam = (c + theta) / 2;
if strcmp(Mo.start, 'stationary')
    [S, h] = dfdet_sigma(P, theta);
    v = -1 / (2 * c); n = P.N; N = n + 1;
    Sx = complex(zeros(N)); Sx(1:n, 1:n) = S + v * (h * h.'); Sx(N, 1:n) = v * h.'; Sx(1:n, N) = v * h; Sx(N, N) = v;
    Q = s * Mo.A + t * Mo.B; Q(1, 1) = Q(1, 1) + gam; Q(N, N) = Q(N, N) - gam;
    M = eye(N) - 2 * (Q * Sx);
    ld = pivlogs(Sx \ M.') + pivlogs(S) + log(-1 / (2 * c));
else
    S = dfdet_sigma(P, theta);
    M = buildM(P, S, c, s, t, theta);
    ld = pivlogs(S \ M.') + pivlogs(S);
end
L = -s / 2 - (c + theta) / 2 - ld / 2;
end

function s = pivlogs(A)
% the sum of the principal logarithms of the pivots of Gaussian elimination without pivoting
n = size(A, 1); s = 0;
for k = 1:n
    p = A(k, k); s = s + log(p);
    if k < n
        l = A(k+1:n, k) * (1 / p);
        A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - l * A(k, k+1:n);
    end
end
end

function M = buildM(P, S, c, s, t, theta)
% I - 2 (s A + t B + gamma E11) Sigma with its (1,1) entry without cancellation (buildM of DFDET_MGF)
N = P.N; gam = (c + theta) / 2;
Q = -2 * ((s * P.A + t * P.B) * S);
M = eye(N) + Q;
M(1, :) = M(1, :) - 2 * gam * S(1, :);
M(1, 1) = m011(c, t, theta) + Q(1, 1);
end
function v = m011(c, t, theta)
% 1 - (c + theta) e_0(2 theta) = 1 - 2 gamma Var J(1), without the cancellation of theta - c for c > 0
if abs(theta) < 0.5
    z = 2 * theta;
    if abs(z) < 1e-3, e = 1 - z / 2 + z * z / 6 - z^3 / 24; else, e = (1 - exp(-z)) / z; end
    v = 1 - (c + theta) * e; return
end
if c > 0, d = -2 * t / (theta + c); else, d = theta - c; end
v = (d + (c + theta) * exp(-2 * theta)) / (2 * theta);
end

function s4 = scale(K)
% means and variances of (A, B) from the log characteristic function next to the origin (no jets for this transform)
m = zeros(2);
for q = 1:2
    a = double(q == 1); b = 1 - a; h = 1e-2;
    L = raw(K, a * h, b * h);
    while abs(L) > 0.05 && h > 1e-12, h = h / 4; L = raw(K, a * h, b * h); end
    m(q, :) = [imag(L) / h, -2 * real(L) / h^2];
end
s4 = [m(1, 1), m(2, 1), sqrt(max(m(2, 2), 0)), m(1, 2)];
end

function [T, e] = tail_int(K, alpha, beta, b, w, S, cdf)
% the tail int_S^inf psi(u) e^{-i w u} du for a transform without jets, psi(u) = phi(u alpha, u beta) e^{-i b u} e^{i w u}/u
% (cdf) or without the 1/u (the density): the derivatives of the slowly varying factor at S from the polynomial through
% 9 Chebyshev points of the real axis about S, four terms of the expansion (tail_int of DFKan_det.jl)
m = 9; L = S / 16; xk = cos(pi * (0:m-1)' / (m - 1));
g = complex(zeros(m, 1));
for k = 1:m
    u = S + L * xk(k);
    z = raw(K, alpha * u, beta * u) - complex(0, b) * u + complex(0, w) * (u - S);
    g(k) = exp(z); if cdf, g(k) = g(k) / u; end
end
cf = (xk .^ (0:m-1)) \ g;                                   % g(S + L x) = sum_j cf(j) x^(j-1)
acc = 0; iw = 1 / complex(0, w); p = iw; KT = 4; fa = exp(gammaln((0:KT)' + 1));
for k = 0:KT - 1, acc = acc + (fa(k + 1) * cf(k + 1)) / L^k * p; p = p * iw; end
T = acc; e = abs((fa(KT + 1) * cf(KT + 1)) / L^KT * p);
end

function [f, err] = half_pdf(K, alpha, beta, b, tol, w0)
% the density of alpha A + beta B at b: (1/pi) int_0^inf Re[phi(u alpha, u beta) e^{-i b u}] du by adaptive quadrature over
% decades and beyond S = 400/|w| by the tail's expansion (half_pdf of DFKan_det.jl)
s4 = K.scale(); vv = alpha^2 * s4(4) + beta^2 * s4(3)^2; if ~(isfinite(vv) && vv > 0), vv = s4(4); end
sc = 1 / sqrt(vv);
fn = @(z) hpf(K, alpha, beta, b, sc, z);
w = w0 * alpha + b; T = 0; et = 0; S = Inf;
if abs(w) > 1e-7
    S = 400 / abs(w);
    for dbl = 0:5
        [Tc, et] = tail_int(K, alpha, beta, b, w, S, false); T = real(Tc);
        if et <= tol / 10 || dbl == 5, break, end
        S = S * 2;
    end
end
[I, e] = kan_core('decades', fn, S / sc, tol, max(tol, 1e-14));
f = (I + T) / pi; err = (e + et) / pi;
end
function y = hpf(K, alpha, beta, b, sc, z)
u = z * sc; L = complex(zeros(size(u)));
for k = 1:numel(u), L(k) = raw(K, alpha * u(k), beta * u(k)); end
y = sc * exp(real(L)) .* cos(imag(L) - b * u);
y(~isfinite(y)) = 0;
end

function [P, e] = pt_cdf1(K, cbar, q, tol)
if isnan(q), P = NaN; e = 0; return, end
if isinf(q), P = double(q > 0); e = 0; return, end
[P, e] = kan_core('half_cdf', K, 1, cbar^2 + K.c, q - 0.5, tol, 0.5);
P = min(max(P, 0), 1);
end
function [f, e] = pt_pdf1(K, cbar, q, tol)
if ~isfinite(q), f = 0; if isnan(q), f = NaN; end, e = 0; return, end
[f, e] = half_pdf(K, 1, cbar^2 + K.c, q - 0.5, tol, 0.5);
f = max(f, 0);
end
function [q, e] = pt_quantile1(K, cbar, p, tol)
M = K; M.memo = containers.Map('KeyType', 'char', 'ValueType', 'any');
M.lcf = @(t1, t2) lcfv(M, t1, t2);
s4 = M.scale(); bq = cbar^2 + K.c;
m = s4(1) + bq * s4(2) + 0.5; s = sqrt(s4(4) + bq^2 * s4(3)^2);
el = 0;
    function v = g(x)
        [Pp, ee] = pt_cdf1(M, cbar, x, tol); el = ee; v = Pp - p;
    end
q = kan_core('illinois', @g, max(m - s, 0), m + s, 1e-12 * s, s);
e = el;
end
