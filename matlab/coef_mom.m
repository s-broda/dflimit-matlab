function [m,err] = coef_mom(p,cas,c,opts)
%COEF_MOM  Moments of the limiting Dickey-Fuller coefficient statistic.
%
%   m = COEF_MOM(p,cas,c,opts): m(j) = E[{(U - c0 V)/V}^j], j = 1..p, the moments of the limit of n(phi^ - phi0) (the
%   coefficient statistic of COEF_CDF; options c0, g0 as there), in the cases DK, DU, DT (zero or fixed start) and SK,
%   SU, ST (a design that is one of them): Sawa's representation m_j = int_0^inf [d^j M(s1, -s2)/ds1^j at s1 = 0]
%   s2^(j-1)/Gamma(j) ds2 with the derivative from one jet of order j per node, the method of Yong Bao and Raymond Kan,
%   "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in First-Order Autoregressions: O(1)
%   Algorithms" (KAN_CORE), adaptive quadrature in log s2 to the relative tolerance opts.tol (1e-10); [m,err] =
%   COEF_MOM(...) returns the estimates.  Up to c = 20 (opts.method = 'inversion' or 'jets' beyond) and p = 31; the
%   engines have no moments of the coefficient statistic (an error for case DG, other designs and opts.engine =
%   'general').
%
%   Example: m = coef_mom(1,'DK'); m(1)   % -1.7814: the mean of n(phi^ - 1) at the unit root
%
%   See also COEF_CDF, COEF_PDF, LTRATIO_MOM, KAN_CORE.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 1), error('ltratio:args','coef_mom: p must be a positive integer.'); end
if ltratio_croute('use',cas,opts,{'tol'},true), [m,err] = ltratio_croute('coefmom',p,cas,c,opts); return, end   % (the C router)
mt = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), mt = lower(char(opts.method)); end
if ~any(strcmp(mt,{'auto','inversion','jets'})), error('ltratio:args','coef_mom: opts.method must be ''auto'' (the jets of Bao and Kan).'); end
[engine,oc] = ltratio_route('engine',opts);
if strcmp(engine,'general'), error('ltratio:args','coef_mom: the general engine has no moments of the coefficient statistic.'); end
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    if isempty(S.case) || strcmp(S.case,'DG'), error('ltratio:args','coef_mom: the moments of the coefficient statistic are available for the cases DK, DU, DT, SK, SU, ST (dfdet_trend(-1), dfdet_trend(0), dfdet_trend(1) with a zero, fixed or stationary start).'); end
    oc.g0 = S.g0; [m,err] = coef_mom(p,S.case,c,oc); return
end
[cs,c,c0,g0] = ltratio_args(cas,c,rmf(oc,{'method','tol','verified'}),false);
if ~any(strcmp(cs,{'DK','DU','DT','SK','SU','ST'})), error('ltratio:args','coef_mom: the moments of the coefficient statistic are available for the cases DK, DU, DT, SK, SU, ST.'); end
if ~(c <= 20) && strcmp(mt,'auto'), error('ltratio:args','coef_mom: c must be at most 20 (opts.method = ''jets'' beyond).'); end
if p > 31, error('ltratio:args','coef_mom: the moments from jets go up to p = 31.'); end
tol = kan_route('tol',opts,1e-10);
[model,stat] = kan_route('model',cs);
[m,err] = kan_core('ratio_mom',kan_core('limcf',model,stat,c,g0),p,c0 - c,1,tol);
end

function o = rmf(o,f)
if ischar(f), f = {f}; end
for k = 1:numel(f), if isfield(o,f{k}), o = rmfield(o,f{k}); end, end
end
