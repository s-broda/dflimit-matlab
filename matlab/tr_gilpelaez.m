function [P,phi] = tr_gilpelaez(b,a,cas,c,opts,phi)
%TR_GILPELAEZ  The distribution of a linear form of (U, V) with a zero start, by the inversion formula of Gil-Pelaez.
%
%   P = TR_GILPELAEZ(b,a,cas,c,opts) is P(U - aV < b) (elementwise in b; a NaN: P(V < b)) under phi = 1 + c/n, zero
%   start, no drift, in case 'DK' (also the constant case with GLS demeaning, whose limit is the same) or 'DG' (the
%   GLS-detrended trend case at the detrending parameter opts.cbar, default -13.5; a finite only).  By
%   P(X < b) = 1/2 - (1/pi) int_0^Inf Im[phi(u) e^(-iub)]/u du from the joint moment generating function with
%   a zero start, M(s,t) = e^(-(s+c)/2) (cosh mu - (c+s) sinh(mu)/mu)^(-1/2), mu^2 = c^2 - 2t, in case DK (with the
%   growth e^mu taken out and the phase unwrapped along u), and e^(-(s+c)/2) (cosh th |I - Q Sig|)^(-1/2)
%   (1 - s/s_1)^(-1/2) from the Girsanov forms of case T in case DG (GLS_PARTS); the trapezoidal rule in u = sinh(v)
%   with opts.N (60000) steps up to opts.umax (4000), from u = 1e-7.  [P,phi] returns the characteristic function on
%   the grid, which a further call with the same a, cas, c accepts as its sixth argument (the bisections of
%   POWER_ENVELOPE).  As tr_gilpelaez of julia/DFExact.jl.
%
%   The transform of U - aV decays like a power of u for a = 0 (R = U + 1/2 is a gamma(1/2) variable, whose density
%   is infinite at its edge) and like exp(-(|a| u)^(1/2)) only from u ~ 1/|a| on, so that the inversion truncated at
%   umax is off by up to 7e-2 near the edge b = -1/2 and 1e-3 in the bulk for |a| up to about 0.3 (more at very
%   negative c, where R and V are concentrated), and for c > 0, where U - aV is wide.  Where |phi| exceeds 5e-9 on the
%   last 101 points of the grid, or c > 0, the probability comes from a conditional integral instead (TR_LINEXACT: the
%   region {U <= aV + b} in case DK, the conditional integral given R in case DG; to about 1e-10).  P(V < b) comes from
%   a longer grid where its transform has not decayed by umax (c <= -50) and from a finer one for c > 2, where V spreads
%   like 1/t_* (the step h0 min(1, (t_*(c)/t_*(2))^(1/2)), h0 = asinh(4000)/60000, up to |phi| < 1e-16; at most 4e6
%   nodes, c up to about 7).  The values are clamped to [0, 1], and b = -Inf, Inf give 0, 1.
%
%   See also LTRATIO_LINCDF, POWER_ENVELOPE, LTRATIO_REGION.
if nargin < 5 || isempty(opts), opts = struct; end
N = getf(opts,'N',60000); umax = getf(opts,'umax',4000);
vmax = asinh(umax); h = vmax/N; v = (0:N)'*h; u = sinh(v); u(1) = 1e-7; w = cosh(v)*h; w(1) = w(1)/2;
if nargin < 6 || isempty(phi), phi = gpphi(u,a,cas,c,opts); end
tail = max(abs(phi(max(1,end-100):end)));
if isnan(a)
    [Nv,uv] = vgrid(c,tail);                          % (the law of V: a longer or finer grid where needed)
    if ~isempty(Nv)
        vmax = asinh(uv); h = vmax/Nv; v = (0:Nv)'*h; u = sinh(v); u(1) = 1e-7; w = cosh(v)*h; w(1) = w(1)/2;
        phi = gpphi(u,a,cas,c,opts);
    end
elseif ~(c <= 0 && tail <= 5e-9)
    % the inversion is not trusted (a transform that has not decayed by umax, or c > 0): a conditional integral
    P = tr_linexact(b,a,upper(cas),c,opts);
    return
end
P = zeros(size(b));
for k = 1:numel(b)
    if ~isfinite(b(k)), if isnan(b(k)), P(k) = NaN; else, P(k) = double(b(k) > 0); end, continue, end
    P(k) = min(1,max(0,0.5 - sum(imag(phi.*exp(-1i*u*b(k)))./u.*w)/pi));
end
end
function phi = gpphi(u,a,cas,c,opts)
% the characteristic function of U - aV (a NaN: of V) on the grid u
switch upper(cas)
    case 'DK'
        if isnan(a), S = zeros(size(u)); T = 1i*u; else, S = 1i*u; T = -a*1i*u; end
        mu = sqrt(c^2 - 2*T); E = exp(-2*mu);
        q = -expm1(-2*mu)./mu; sm = abs(mu) < 1e-6; q(sm) = 2 - 2*mu(sm);
        lm = mu - log(2); br = (1 + E) - (c + S).*q;
        ph = unwrap(angle(br));
        phi = exp(-(S + c)/2 - (lm + log(abs(br)) + 1i*ph)/2);
    case 'DG'
        if isnan(a), error('tr_gilpelaez:args','tr_gilpelaez: P(V < b) is available in case DK only.'); end
        s = 1i*u; t = -a*s;
        [lch,D,s1] = gls_parts(t,c,getf(opts,'cbar',-13.5));
        Dt = exp(lch).*D; F = 1 - s./s1;
        pd = unwrap(angle(Dt)); pf = unwrap(angle(F));
        phi = exp(-(s + c)/2).*abs(Dt).^(-0.5).*exp(-1i*pd/2).*abs(F).^(-0.5).*exp(-1i*pf/2);
    otherwise
        error('tr_gilpelaez:args','tr_gilpelaez: the case must be DK or DG.');
end
end
function [N,umax] = vgrid(c,tail)
% the grid of P(V < b) where the default one does not serve: c <= 2 and |phi| < 1e-9 at its end ([] then); as vgrid
% of julia/DFExact.jl (dfl_vgrid of the C engine)
N = []; umax = [];
if c <= 2 && tail <= 1e-9, return, end
h0 = asinh(4000)/60000; h = h0;
if c > 2, h = h0*min(1,sqrt(dfx_tstar(1,c)/dfx_tstar(1,2))); end
U = 4000;
while absphiV(c,U) > 1e-16 && U < 1e9, U = U*1.5; end
N = ceil(asinh(U)/h);
if N > 4e6, error('ltratio:args','ltratio_lincdf: P(V <= b) at c = %g needs more than 4000000 nodes in the inversion of Gil-Pelaez (case DK is available for c up to about 7).',c); end
umax = sinh(N*h);
end
function p = absphiV(c,u)
mu = sqrt(complex(c^2,-2*u)); q = -expm1(-2*mu)/mu; if abs(mu) < 1e-6, q = 2 - 2*mu; end
br = (1 + exp(-2*mu)) - c*q;
p = exp(-c/2 - (real(mu) - log(2) + log(abs(br)))/2);
end
function v = getf(o,f,d), if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end, end
