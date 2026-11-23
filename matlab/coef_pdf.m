function [f,err] = coef_pdf(y,cas,c,opts)
%COEF_PDF  Limiting density of the Dickey-Fuller coefficient statistic.
%
%   f = COEF_PDF(y,cas,c,opts) is the limiting density of n(phi^ - phi0), the coefficient statistic (U - c0 V)/V of
%   COEF_CDF (y an array; options c0, g0 as there).  In the cases DK, DU, DT (zero or fixed start) and SK, SU, ST up to
%   c = 20 it is the single inversion of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the
%   Least Squares Estimator in First-Order Autoregressions: O(1) Algorithms" (KAN_CORE): the density of a ratio by the
%   formula of Geary, one integral of the derivative of the joint characteristic function along the ray, from jets,
%   its tail completed by integration by parts.  opts.tol (1e-12) is its tolerance relative to the scale
%   max(1, E[V]/sd(U)) of the density, [f,err] = COEF_PDF(...) returns the estimates, and the warning ltratio:inversion
%   names the points at which they exceed the tolerance.  opts.method = 'inversion' forces it; opts.method = 'engine'
%   (and opts.engine = 'general', case DG, a design beyond the classical cases) is the general engine's density, the
%   kernel 2 v^3 of its contour integral (DFDET_COEF_PDF; err NaN).  opts.verified = true recomputes the values whose
%   estimate exceeds the tolerance by the general engine's checked computation.  The classical programs have no
%   density of the coefficient statistic: with opts.engine = 'classical' the single inversion computes.
%
%   Example: coef_pdf(-8.039,'DK')        % the density at the 5% critical value of rho^
%
%   See also COEF_CDF, COEF_ICDF, COEF_MOM, DFDET_COEF_PDF, KAN_CORE.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts,{'tol'},true), [f,err] = ltratio_croute('coefpdf',y,cas,c,opts); return, end   % (the C router)
m = kan_route('method',opts,'coef_pdf',false);
[engine,oc] = ltratio_route('engine',opts);
S = ltratio_spec(cas,opts);
cs = S.case;
if ~(ischar(cas) || isstring(cas))
    if strcmp(engine,'classical') && isempty(cs), error('ltratio:args','coef_pdf: engine = ''classical'': the classical programs have the cases DK, DU, DT, SK, SU, ST and DG only.'); end
    if ~strcmp(engine,'general') && ~isempty(cs) && (strcmp(engine,'classical') || ltratio_route('valid',S,c,rmf(oc,{'method','tol','verified'})))
        oc.g0 = S.g0; oc.engine = engine; [f,err] = coef_pdf(y,cs,c,oc); return
    end
    if strcmp(m,'inversion'), kan_route('take','coef_pdf',m,engine,'',c,true); end
    f = gen(y,S,c,oc); err = NaN(size(y)); return
end
if ~(isnumeric(y) && isreal(y)), error('ltratio:args','coef_pdf: y must be real.'); end
c0 = 0; g0 = 0;
if ~strcmp(engine,'general'), [cs,c,c0,g0] = ltratio_args(cas,c,rmf(oc,{'method','tol','verified'}),false); end
eng = engine; if strcmp(eng,'classical'), eng = 'auto'; end
if kan_route('take','coef_pdf',m,eng,cs,c,kan_route('plain',oc))
    tol = kan_route('tol',opts,1e-12);
    [model,stat] = kan_route('model',cs); K = kan_core('limcf',model,stat,c,g0);
    f = zeros(size(y)); err = zeros(size(y)); T = zeros(size(y));
    for k = 1:numel(y), [f(k),err(k),T(k)] = kan_core('coef_pdf1',K,y(k),c0,tol); end
    bad = ~(err <= T);
    if isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified)
        if any(bad(:)), o = oc; o.verified = true; f(bad) = gen(y(bad),S,c,o); end
    else
        kan_route('warn','coef_pdf',err,T);
    end
    return
end
f = gen(y,S,c,oc); err = NaN(size(y));
end

function f = gen(y,S,c,oc)
% the general engine's density (DFDET_COEF_PDF), a negative value set to zero
f = dfdet_coef_pdf(y,S.D,c,ltratio_route('gopts',S,rmf(oc,{'method','tol'})));
f(f < 0) = 0;
end

function o = rmf(o,f)
if ischar(f), f = {f}; end
for k = 1:numel(f), if isfield(o,f{k}), o = rmfield(o,f{k}); end, end
end
