function [y,err] = coef_icdf(p,cas,c,opts)
%COEF_ICDF  Quantiles of the limiting Dickey-Fuller coefficient statistic.
%
%   y = COEF_ICDF(p,cas,c,opts): the p-quantiles (p an array in (0, 1)) of the limit of n(phi^ - phi0) (COEF_CDF),
%   by the Illinois variant of regula falsi to the width opts.tol (1e-9).  By default (in the cases and the range of the
%   single inversion of COEF_CDF, the method of Bao and Kan, "Exact Finite-Sample Distribution and Moments of the Least
%   Squares Estimator in First-Order Autoregressions: O(1) Algorithms") on the single inversion's distribution function
%   (tolerance 1e-12), from the bracket E[U - cV]/E[V] -/+ twice the scale of the ratio; [y,err] = COEF_ICDF(...) returns
%   the error estimate of the distribution function at the quantile divided by the density there, and the warning
%   ltratio:inversion names the points at which the estimate of the distribution function exceeds 1e-12.  opts.method =
%   'engine' (and the t-ratio's representations, opts.engine 'classical' or 'general', case DG, options of the engines)
%   is the search of before on COEF_CDF with opts.method = 'engine' from a bracket grown from [-40, 10], about ten
%   evaluations; opts.verified = true that search on the checked computation where the single inversion's estimate
%   exceeds its tolerance.
%
%   cas may be a design of the general engine (LTRATIO_SPEC; a design that is
%   not a classical case: DFDET_COEF_ICDF).
%
%   Example: coef_icdf(0.05,'DK')         % -8.039, the 5% critical value of rho^
%
%   See also COEF_CDF, COEF_PDF, LTRATIO_ICDF, DFDET_COEF_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts,{'tol'},true), [y,err] = ltratio_croute('coeficdf',p,cas,c,opts); return, end   % (the C router)
m = kan_route('method',opts,'coef_icdf',true);
if ~(ischar(cas) || isstring(cas))
    [engine,oc] = ltratio_route('engine',opts); S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && isempty(S.case), error('ltratio:args','coef_icdf: engine = ''classical'': the classical programs have the cases DK, DU, DT, SK, SU, ST and DG only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && (strcmp(engine,'classical') || ltratio_route('valid',S,c,rmf(oc,{'tol','method'})))
        oc.g0 = S.g0; oc.engine = engine; [y,err] = coef_icdf(p,S.case,c,oc); return
    end
    if strcmp(m,'inversion'), kan_route('take','coef_icdf',m,engine,'',c,true); end
    go = ltratio_route('gopts',S,rmf(oc,{'method'})); if isfield(oc,'tol'), go.tol = oc.tol; end
    y = dfdet_coef_icdf(p,S.D,c,rmf(go,'method')); err = NaN(size(p)); return
end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','coef_icdf: p must lie in (0,1).'); end
xtol = 1e-9; if isfield(opts,'tol') && ~isempty(opts.tol), xtol = opts.tol; end
if ~(isnumeric(xtol) && isscalar(xtol) && xtol > 0 && isfinite(xtol)), error('ltratio:args','coef_icdf: opts.tol must be a positive finite number.'); end
[engine,oc] = ltratio_route('engine',opts);
o = rmf(opts,{'method','tol'}); if any(strcmp(m,{'single','double','contour'})), o.method = m; end
[cs,c,c0,g0] = ltratio_args(cas,c,rmf(oc,{'method','tol','verified'}),~any(strcmp(m,{'auto','inversion'})));
if kan_route('take','coef_icdf',m,engine,cs,c,kan_route('plain',oc))
    [model,stat] = kan_route('model',cs); K = kan_core('limcf',model,stat,c,g0);
    y = zeros(size(p)); e = zeros(size(p));
    for k = 1:numel(p), [y(k),e(k)] = kan_core('coef_icdf1',K,p(k),c0,xtol,1e-12); end
    if isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified)
        bad = ~(e <= 1e-12);
        if any(bad(:)), o.verified = true; o.method = 'engine'; y(bad) = coef_icdf(p(bad),cas,c,setfield(o,'tol',xtol)); e(bad) = NaN; end %#ok<SFLD>
    else
        kan_route('warn','coef_icdf',e,1e-12);
    end
    if nargout > 1                                  % (in y: the estimate of F divided by the density)
        err = NaN(size(p));
        for k = 1:numel(p), if ~isnan(e(k)), err(k) = e(k) / kan_core('coef_pdf1',K,y(k),c0,1e-12); end, end
    end
    return
end
% the search of before on the engines' distribution function (COEF_CDF with opts.method = 'engine' or a representation)
if ~isfield(o,'method'), o.method = 'engine'; end
y = zeros(size(p)); err = NaN(size(p));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (each brittle regime warned once)
for k = 1:numel(p)
    f = @(x) coef_cdf(x,cas,c,o) - p(k);
    lo = -40; hi = 10; flo = f(lo); fhi = f(hi);
    while flo > 0, hi = lo; fhi = flo; lo = 2*lo - 10; flo = f(lo); end
    while fhi < 0, lo = hi; flo = fhi; hi = 2*hi + 10; fhi = f(hi); end
    side = 0; x = (lo + hi)/2;
    for it = 1:200
        x = (lo*fhi - hi*flo)/(fhi - flo); fx = f(x);
        if fx == 0, break, end
        if fx < 0, lo = x; flo = fx; if side == -1, fhi = fhi/2; end, side = -1;
        else, hi = x; fhi = fx; if side == 1, flo = flo/2; end, side = 1; end
        if hi - lo < xtol, x = (lo*fhi - hi*flo)/(fhi - flo); break, end
    end
    y(k) = x;
end
end

function o = rmf(o,f)
if ischar(f), f = {f}; end
for k = 1:numel(f), if isfield(o,f{k}), o = rmfield(o,f{k}); end, end
end
