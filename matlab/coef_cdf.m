function [F,err] = coef_cdf(y,cas,c,opts)
%COEF_CDF  Limiting distribution function of the Dickey-Fuller coefficient statistic.
%
%   F = COEF_CDF(y,cas,c,opts) is the limit as n -> Inf of P[n(phi^ - phi0) <= y] (y an array) for the least squares
%   estimate phi^ of phi in the regression of LTRATIO_CDF, under phi = 1 + c/n, with phi0 = 1 + c0/n (opts.c0,
%   default 0: the coefficient test of Dickey and Fuller, rho^ = n(phi^ - 1)).  The limit is (U - c0 V)/V, so that
%   P[(U - c0 V)/V <= y] = P[U - (c0 + y)V <= 0].
%
%   In the cases DK, DU, DT (zero or fixed start opts.g0) and SK, SU, ST the default (opts.method 'auto' with
%   opts.engine 'auto', up to c = 20) is the single inversion of the joint characteristic function along the ray of the
%   half-plane, the method of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least
%   Squares Estimator in First-Order Autoregressions: O(1) Algorithms" (KAN_CORE): adaptive quadrature over decades,
%   the tail completed by integration by parts with the derivatives from jets, the closed-form transform with the
%   branch fixed at each point, one evaluation with an error estimate per value; opts.tol (1e-12) its absolute
%   tolerance, [F,err] = COEF_CDF(...) returns the estimates, and the warning ltratio:inversion names the points at
%   which they exceed it.  opts.method = 'inversion' forces it (also above c = 20); opts.method = 'engine' is the
%   computation of before, the distribution function of the t-ratio at 0 with the local null value c0 + y,
%   LTRATIO_CDF(0,cas,c,opts) with c0 + y in place of c0, with its choice of engine and representation (err NaN); the
%   t-ratio's representations ('single' or 'line', 'double' or 'conditional', 'contour'), opts.engine 'classical' or
%   'general', options of the engines and case DG take that computation too.  opts.verified = true keeps each value of
%   the single inversion whose estimate meets opts.tol and recomputes the others by the engines' checked computation.
%
%   opts.engine is passed to LTRATIO_CDF (each value routed there); cas may be
%   a design of the general engine (LTRATIO_SPEC; a design that is not a
%   classical case: DFDET_COEF_CDF).
%
%   Example: coef_cdf(-14.094,'DU',0)    % 0.05: the 5% critical value of rho^_mu
%
%   See also COEF_ICDF, COEF_PDF, COEF_MOM, LTRATIO_CDF, DFDET_COEF_CDF, KAN_CORE.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts,{'tol'},true), [F,err] = ltratio_croute('coef',y,cas,c,opts); return, end   % (the C router)
m = kan_route('method',opts,'coef_cdf',true);
if ~(ischar(cas) || isstring(cas))
    [engine,oc] = ltratio_route('engine',opts); S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && isempty(S.case), error('ltratio:args','coef_cdf: engine = ''classical'': the classical programs have the cases DK, DU, DT, SK, SU, ST and DG only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && (strcmp(engine,'classical') || ltratio_route('valid',S,c,rmf(oc,{'method','tol'})))
        oc.g0 = S.g0; oc.engine = engine; [F,err] = coef_cdf(y,S.case,c,oc); return
    end
    if strcmp(m,'inversion'), kan_route('take','coef_cdf',m,engine,'',c,true); end
    F = dfdet_coef_cdf(y,S.D,c,ltratio_route('gopts',S,rmf(oc,{'method','tol'}))); err = NaN(size(y)); return
end
if ~(isnumeric(y) && isreal(y)), error('ltratio:args','coef_cdf: y must be real.'); end
[engine,oc] = ltratio_route('engine',opts);
o = rmf(opts,{'method','tol'}); if any(strcmp(m,{'single','double','contour'})), o.method = m; end
[cs,c,c0,g0] = ltratio_args(cas,c,rmf(oc,{'method','tol','verified'}),~any(strcmp(m,{'auto','inversion'})));
if kan_route('take','coef_cdf',m,engine,cs,c,kan_route('plain',oc))
    % the single inversion of Bao and Kan
    tol = kan_route('tol',opts,1e-12);
    [model,stat] = kan_route('model',cs); K = kan_core('limcf',model,stat,c,g0);
    F = zeros(size(y)); err = zeros(size(y));
    for k = 1:numel(y), [F(k),err(k)] = kan_core('coef_cdf1',K,y(k),c0,tol); end
    if isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified)
        bad = ~(err <= tol);
        if any(bad(:)), o.verified = true; F(bad) = coef_engine(y(bad),cs,c,c0,g0,o); end
    else
        kan_route('warn','coef_cdf',err,tol);
    end
    return
end
F = coef_engine(y,cas,c,c0,g0,o); err = NaN(size(y));
end

function F = coef_engine(y,cas,c,c0,g0,opts)
% the engines' computation: the t-ratio at 0 with the null value c0 + y (with its final check, named as this call)
F = zeros(size(y));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (the regimes c0 + y passes, each warned once)
for k = 1:numel(y)
    if isnan(y(k)), F(k) = NaN; elseif y(k) == Inf, F(k) = 1; elseif y(k) == -Inf, F(k) = 0;
    else, o = opts; o.c0 = c0 + y(k); F(k) = ltratio_dist(0,cas,c,o,'cdf'); end
end
a = sprintf('case %s, c = %g',upper(char(cas)),c); if c0 ~= 0, a = [a sprintf(', c0 = %g',c0)]; end
if g0 ~= 0, a = [a sprintf(', g0 = %g',g0)]; end
F = ltratio_rangecheck('coef_cdf',a,'y',y,F);
end

function o = rmf(o,f)
if ischar(f), f = {f}; end
for k = 1:numel(f), if isfield(o,f{k}), o = rmfield(o,f{k}); end, end
end
