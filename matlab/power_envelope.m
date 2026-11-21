function P = power_envelope(c,cas,opts)
%POWER_ENVELOPE  The asymptotic power envelope of Elliott, Rothenberg and Stock (1996).
%
%   P = POWER_ENVELOPE(c,cas,opts) is the power at c (an array) of the point-optimal level-alpha test against c
%   itself (opts.alpha, default 0.05), the envelope of all tests of the unit root, zero start, no drift: case 'DK'
%   (no deterministic terms; also the constant case, case 'DU', where the point-optimal test uses the GLS-demeaned
%   data and the envelope is the same curve) and case 'DG' (the trend case, the point-optimal test against c with the
%   data GLS-detrended at c).  The test against cbar rejects for small U - aV (LTRATIO_LINCDF), a = cbar/2 in case
%   'DK' and a = cbar^2/(2 (cbar - 1)) in case 'DG' (where the limit of the invariant statistic of Elliott,
%   Rothenberg and Stock is cbar^2 V + (1 - cbar) J(1)^2, J(1)^2 = 2U + 1): its critical value is the
%   alpha-quantile at c = 0 (bisection on [-50, 50] to 1e-9) and the envelope its power at cbar = c; alpha for
%   c >= 0.  To about 1e-10; 0.5 s per value in case DK, 3 s in case DG.  Near c = 0, where a is small and the
%   inversion of Gil-Pelaez is not trusted (TR_GILPELAEZ), the critical value and the power come from the
%   conditional integral (TR_LINEXACT), the critical value by bisection in x = b + 1/2 to 1e-12 of the scale
%   s = E(U + 1/2) + |a| E V of U + 1/2 - aV (in case DG detrended at cbar = c the law of U + 1/2 shrinks like c^4,
%   so that a tolerance in b of 1e-9 does not resolve it).  As power_envelope of
%   julia/DFExact.jl (the data of the paper's figures of local power and of the weighted-average-power test).
%
%   Example: power_envelope(-7,'DK')    % 0.4994: at cbar = -7 the envelope has power one half
%
%   See also LTRATIO_LINCDF, LTRATIO_CDF, COEF_CDF.
if nargin < 2 || isempty(cas), cas = 'DK'; end
if nargin < 3 || isempty(opts), opts = struct; end
cas = upper(char(cas)); if strcmp(cas,'DU'), cas = 'DK'; end
if ~any(strcmp(cas,{'DK','DG'})), error('ltratio:args','power_envelope: the case must be DK (or DU) or DG.'); end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
P = zeros(size(c));
for k = 1:numel(c)
    if ~(c(k) < 0), P(k) = alpha; continue, end
    if strcmp(cas,'DK'), a = c(k)/2; else, a = c(k)^2/(2*(c(k) - 1)); end   % the point-optimal weight
    o = opts; o.cbar = c(k);
    [~,phi0] = tr_gilpelaez(0,a,cas,0,o);
    if isempty(phi0), tail = gptail(a,cas,o); else, tail = max(abs(phi0(end-100:end))); end
    if ~(tail <= 5e-9)
        % near c = 0: the conditional integral
        if isempty(phi0)                               % (the C engine, through the drop-in TR_GILPELAEZ of c/mex)
            P(k) = ltratio_mex('envelope',c(k),1 + 3*strcmp(cas,'DG'),alpha,[getf(o,'N',NaN) getf(o,'umax',NaN)]);
        else
            P(k) = tr_linexact(lin_quantile_x(alpha,a,cas,o),a,cas,c(k),o,true);
        end
        continue
    end
    b = ltratio_bisect(@(b) tr_gilpelaez(b,a,cas,0,o,phi0) - alpha,-50,50,1e-9);
    P(k) = tr_gilpelaez(b,a,cas,c(k),o);
end
end

function x = lin_quantile_x(p,a,cas,o)
% x = b + 1/2 with P(U - aV <= b) = p at c = 0 by the conditional integral (lin_quantile_x of julia/DFExact.jl)
if strcmp(cas,'DK'), kap = 1; ev = 0.5; else, [~,L] = tr_lincond([],a,'DG',0,o.cbar); kap = L.kap; ev = L.ev; end
s = 1/(2*kap) + abs(a)*ev;
F = @(x) tr_linexact(x,a,cas,0,o,true) - p;
lo = -s; if a < 0, lo = 0; end
hi = s;
for it = 1:200, if F(hi) >= 0, break, end, hi = hi*2; end
for it = 1:200, if a < 0 || F(lo) <= 0, break, end, lo = lo*2; end
x = ltratio_bisect(F,lo,hi,1e-12*s);
end

function m = gptail(a,cas,o)
% max |phi| of U - aV at c = 0 on the last 101 points of the grid of TR_GILPELAEZ (the test of its trust)
N = getf(o,'N',60000); umax = getf(o,'umax',4000); h = asinh(umax)/N; u = sinh((N-100:N)'*h);
if strcmp(cas,'DK')
    S = 1i*u; T = -a*1i*u; mu = sqrt(-2*T);
    q = -expm1(-2*mu)./mu; sm = abs(mu) < 1e-6; q(sm) = 2 - 2*mu(sm);
    br = (1 + exp(-2*mu)) - S.*q;
    m = max(exp(-(real(mu) - log(2) + log(abs(br)))/2));
else
    s = 1i*u; [lch,D,s1] = gls_parts(-a*s,0,o.cbar);
    m = max(abs(exp(lch).*D).^(-0.5).*abs(1 - s./s1).^(-0.5));
end
end
function v = getf(o,f,d), if isfield(o,f) && ~isempty(o.(f)), v = double(o.(f)); else, v = d; end, end
