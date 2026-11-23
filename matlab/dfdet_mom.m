function [m,err] = dfdet_mom(p,D,c,opts)
%DFDET_MOM  Moments of the limiting t-ratio for any deterministic terms (the general engine).
%
%   [m,err] = DFDET_MOM(p,D,c,opts): m(j) = E tau^j, j = 1..p, of
%   tau = (U - c0 V)/sqrt(V) (arguments and options c0, start, g0, gls as in
%   DFDET_CDF).  By the gamma integral V^(-j/2) = Gamma(j/2)^-1 int_0^Inf
%   s^(j/2-1) e^(-sV) ds, E tau^j = Gamma(j/2)^-1 int_0^Inf s^(j/2-1) m_j(s) ds,
%   m_j(s) = d^j/du^j M(u, -s - c0 u) at u = 0.  The derivatives are Cauchy
%   integrals on a circle of radius r(s) (trapezoidal rule, opts.Nu = 64
%   nodes) with the logarithm of the determinant continued around the circle
%   from its real point u = r (a circle on which it does not close encloses a
%   zero, and its radius is halved); r(s) = 0.4 min |s_k(-s)| (the roots at
%   the real point t = -s), at most 0.4 (s + t_*)/|c0| with c0 ~= 0, at most 0.25
%   with a fixed or stationary start, and at most 2 for c < -20 and with GLS
%   detrending.  The integral over s is the trapezoidal rule in log s with
%   step opts.h = 0.05 on [min(-60, -2c - 40), max(8, log(80|c|))]; err is the
%   difference to the rule with step 2h.  With c0 = 0 every point of a circle
%   has t = -s, and the determinant there is det M0 prod_k (1 - u mu_k) from the
%   roots at t = -s (DFDET_MGF); with c0 ~= 0 each point has its own t and its
%   own covariance, at Nu times the cost.  Not for the statistic that detrends
%   first.
%
%   Example: dfdet_mom(2, dfdet_perron('C', 0.3))   % the first two moments with a break
%
%   See also DFDET_CDF, LTRATIO_MOM.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[c,c0,start,g0,gls] = dfdet_args(D,c,opts);
if strcmp(D.detrend,'first'), error('dfdet:args','dfdet_mom: %s.',dfdet_msg('first_mom')); end
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 1), error('dfdet:args','dfdet_mom: p must be a positive integer.'); end
Nu = getf(opts,'Nu',64); h = getf(opts,'h',0.05);
Mo = dfdet_startmodel(D,c,start,g0,gls);
tlo = min(-60,-2*c - 40); if c < -20, thi = max(8,log(80*abs(c))); else, thi = 8; end
% (the grid grows with c above the unit root: 2c/h points, 4e5 at c = 1e4; beyond 1e6 an error)
if ~((thi - tlo)/h <= 1e6), error('dfdet:args','dfdet_mom: the grid of the integral in s would have %.3g points at c = %g (more than 1e6): c is beyond the range of the moments.',(thi - tlo)/h,c); end
sv = exp(tlo:h:thi); ns = numel(sv);
smax = dfdet_tstar(Mo,c);
wk = exp(2i*pi*(0:Nu-1)'/Nu);
acc = zeros(p,ns);
fixed = strcmp(Mo.start,'fixed');
for i = 1:ns
    s = sv(i);
    [th0,ld0,mu0,fx0] = dfdet_mgf(Mo,c,complex(-s));
    r = 0.4/max(abs(mu0));
    if c0 ~= 0, r = min(r,0.4*(s + smax)/abs(c0)); end
    if ~strcmp(Mo.start,'zero'), r = min(r,0.25); end
    if c < -20, r = min(r,2); end
    if ~isnan(Mo.P.cbar), r = min(r,2); end
    ok = false;
    for it = 1:8
        U = r*wk;
        if c0 == 0
            % (every point has t = -s: the determinant from the roots there)
            L0 = -U/2 - (c + th0)/2;
            if fixed, L0 = L0 + Mo.g0^2*(fx0.cf + dfdet_fixe1(fx0,U) - Mo.kappa*U); end
            LD = (ld0 - c - th0) + sum(log(1 - U*mu0.'),2);
        else
            L0 = complex(zeros(Nu,1)); LD = L0;
            cut = -s + abs(c0)*r > c*c/2;              % (whether the circle in t = -s - c0 u reaches the cut of theta)
            for k = 1:Nu
                tk = -s - c0*U(k);
                [L0(k),LD(k)] = dfdet_logmgf(Mo,c,U(k),tk);
                % (where the circle reaches the cut, the logarithm continued around it is that of the entire function
                % e^{c + theta} det(...), real and positive at the real point u = r, and not that of the determinant
                % alone, whose phase jumps with theta across the cut: where -s - c0 r > c^2/2, as at c = 0 with c0 < 0
                % for every s below |c0| r, theta is imaginary at the real point, the check of the real start failed on
                % every circle, and the moments refused (case DT at c = 0, c0 = -2); elsewhere as before)
                if cut
                    e = c + sqrt(complex(c*c - 2*tk));
                    LD(k) = LD(k) + e; L0(k) = L0(k) + e/2;
                end
            end
        end
        [ld,ok] = circlog(exp(LD - real(LD(1))));
        L = L0 - (ld + real(LD(1)))/2;
        if ok, break, end
        r = r/2;
    end
    if ~ok, error('dfdet:moments','dfdet_mom: no circle in u free of zeros of the determinant at s = %g.',s); end
    f = exp(L);
    for j = 1:p
        Mj = factorial(j)/(Nu*r^j)*real(sum(f.*wk.^(-j)));
        acc(j,i) = s^(j/2)*Mj;
    end
end
m = zeros(1,p); err = m;
for j = 1:p
    m(j) = h*sum(acc(j,:))/gamma(j/2);
    err(j) = abs(h*sum(acc(j,:)) - 2*h*sum(acc(j,1:2:ns)))/gamma(j/2);
end
end

function [l,ok] = circlog(Z)
% the logarithm along a closed circle, continued from its first point (where it is real); ok: the steps are resolved and
% the circle closes
lz = log(Z); n = numel(Z);
ok = abs(imag(lz(1))) < 1e-8;
ph = zeros(n,1); ph(1) = imag(lz(1));
d = angle(exp(1i*(imag(lz(2:n)) - imag(lz(1:n-1)))));
if any(abs(d) >= pi/3), ok = false; end
ph(2:n) = ph(1) + cumsum(d);
ok = ok && abs(ph(n) + angle(Z(1)/Z(n)) - ph(1)) < 1e-6;
l = complex(real(lz),ph);
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
