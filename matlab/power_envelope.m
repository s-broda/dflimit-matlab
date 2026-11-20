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
%   c >= 0.  To about 1e-10; 0.5 s per value in case DK, 3 s in case DG.  As power_envelope of
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
    b = ltratio_bisect(@(b) tr_gilpelaez(b,a,cas,0,o,phi0) - alpha,-50,50,1e-9);
    P(k) = tr_gilpelaez(b,a,cas,c(k),o);
end
end
