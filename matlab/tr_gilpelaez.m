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
%   See also LTRATIO_LINCDF, POWER_ENVELOPE, LTRATIO_REGION.
if nargin < 5 || isempty(opts), opts = struct; end
N = getf(opts,'N',60000); umax = getf(opts,'umax',4000);
vmax = asinh(umax); h = vmax/N; v = (0:N)'*h; u = sinh(v); u(1) = 1e-7; w = cosh(v)*h; w(1) = w(1)/2;
if nargin < 6 || isempty(phi)
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
P = zeros(size(b));
for k = 1:numel(b)
    P(k) = 0.5 - sum(imag(phi.*exp(-1i*u*b(k)))./u.*w)/pi;
end
end
function v = getf(o,f,d), if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end, end
