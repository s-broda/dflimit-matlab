function [m,info] = ltratio_mom(p,cas,c,opts)
%LTRATIO_MOM  Moments of the limiting Dickey-Fuller t-ratio.
%
%   m = LTRATIO_MOM(p,cas,c,opts): m(j) = E[tau^j], j = 1..p, for the limit
%   tau = (U - c0 V)/sqrt(V) of the t-ratio under phi = 1 + c/n (arguments
%   and options c0, g0 as in LTRATIO_CDF).  By the gamma integral
%   V^(-j/2) = Gamma(j/2)^(-1) int_0^Inf s^(j/2-1) exp(-s V) ds the moments
%   are one-dimensional integrals of the joint moment generating function of
%   (U, V) at real arguments (the moments of the limit in the Supplementary Material), where it is
%   smooth and positive; the derivatives in the argument of U are Cauchy
%   integrals.  Engine: TR_LIMMOM (the zero-start moment generating function
%   from the Girsanov closed forms, free of cancellation for c > 0; the
%   initial-value terms of a fixed or a stationary start in closed form,
%   TR_NSSTART).  DFX_MOMENTS is an independent implementation for c0 = 0
%   with a zero or a stationary start, used by LTRATIO_TEST as a check.
%   info: the error estimate of the integral over s per moment, and the time.
%   For c < -20 the circles have a radius of at most 2 and the integral over
%   s extends to 80|c| (TR_LIMMOM).  Far from the unit root with c0 = c the
%   third and fourth moments converge slowly in the step of the integral over
%   s: about 1e-9 relative with a zero start, 1e-8 to 2e-7 with a stationary
%   start at c = -100 to -480; the first two moments to 1e-11.
%
%   By default (opts.method 'auto' with opts.engine 'auto', the six classical
%   cases, any c0 and fixed start, c up to 20 and p up to 31) the derivative
%   comes from one Taylor-mode jet of order j at each node of an adaptive
%   quadrature in log s, the method of Yong Bao and Raymond Kan, "Exact
%   Finite-Sample Distribution and Moments of the Least Squares Estimator in
%   First-Order Autoregressions: O(1) Algorithms" (KAN_CORE), with the
%   closed-form transform of the six cases: opts.tol (1e-10) is the relative
%   tolerance of its quadrature, info.err the estimates, and the warning
%   ltratio:inversion names the moments whose estimate exceeds opts.tol times
%   max(|m_j|, m_2^(j/2)); opts.method = 'jets' forces it, opts.method =
%   'engine' is the computation of before (the Cauchy circles below, or the
%   general engine by the rule of LTRATIO_ENGINE), which reproduces those
%   values; opts.engine 'classical' or 'general' and options of the engines
%   take that computation too.  opts.verified = true recomputes the moments by
%   the engines where an estimate exceeds the tolerance.
%
%   opts.engine: 'classical' TR_LIMMOM; 'general' the moments of the general
%   engine (DFDET_MOM); 'auto' (the default) the engine of LTRATIO_ENGINE.  cas
%   may be a design of the general engine (LTRATIO_SPEC).
%
%   Example: m = ltratio_mom(4,'DU');  m(2)   % 3.0559953450 (Nabeya, 1999)
%
%   See also LTRATIO_CDF, LTRATIO_ICDF, TR_LIMMOM, DFX_MOMENTS, DFDET_MOM.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts,{'tol'},true), [m,info] = ltratio_croute('mom',p,cas,c,opts); return, end   % (the C router)
[engine,oc] = ltratio_route('engine',opts);
mt = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), mt = lower(char(opts.method)); end
if strcmp(mt,'inversion'), mt = 'jets'; end
if ~any(strcmp(mt,{'auto','jets','engine'})), error('ltratio:args','ltratio_mom: opts.method must be ''auto'', ''jets'' (one jet per node, after Bao and Kan) or ''engine'' (the engines'' circles).'); end
if strcmp(mt,'jets') && strcmp(engine,'general'), error('ltratio:args','ltratio_mom: opts.method = ''jets'' (the jets of Bao and Kan) does not apply with engine = ''general''.'); end
oc = rmf(oc,{'method','tol'});
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && (isempty(S.case) || strcmp(S.case,'DG')), error('ltratio:args','ltratio_mom: engine = ''classical'': the classical programs have the moments of the cases DK, DU, DT, SK, SU, ST only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && ~strcmp(S.case,'DG') && (strcmp(engine,'classical') || ltratio_route('valid',S,c,rmf(oc,'verified')))
        o2 = opts; o2.g0 = S.g0; o2.engine = engine; o2 = rmf(o2,{'start','gls'});
        [m,info] = ltratio_mom(p,S.case,c,o2);
        return
    end
    if strcmp(mt,'jets'), error('ltratio:args','ltratio_mom: opts.method = ''jets'' applies to the cases DK, DU, DT, SK, SU, ST.'); end
    [m,err] = dfdet_mom(p,S.D,c,ltratio_route('gopts',S,oc)); info = struct('err',err,'engine','general');
    return
end
if strcmp(engine,'general')
    S = ltratio_spec(cas,oc);
    [m,err] = dfdet_mom(p,S.D,c,ltratio_route('gopts',S,oc)); info = struct('err',err,'engine','general');
    return
end
[cas,c,c0,g0] = ltratio_args(cas,c,rmf(oc,'verified'),false);
if strcmp(cas,'DG'), error('ltratio:args','ltratio_mom: the moments are not available for case DG.'); end
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 1), error('ltratio:args','ltratio_mom: p must be a positive integer.'); end
ic = any(strcmp(cas,{'DK','DU','DT','SK','SU','ST'}));
if strcmp(mt,'jets') || (strcmp(mt,'auto') && strcmp(engine,'auto') && kan_route('plain',opts) && ic && c <= 20 && p <= 31)
    % the moments from jets (Bao and Kan; tratio_mom1 of DFLimit.jl)
    if p > 31, error('ltratio:args','ltratio_mom: the moments from jets go up to p = 31 (opts.method = ''engine'', the engines'' circles, beyond).'); end
    tol = kan_route('tol',opts,1e-10); t0 = tic;
    [model,stat] = kan_route('model',cas);
    [mm,ee] = kan_core('ratio_mom',kan_core('limcf',model,stat,c,g0),max(p,2),c0 - c,0.5,tol);
    % (the tolerance relative to max(|m_j|, m_2^(j/2)): a mean near zero, as with c0 = c far above the unit root, has no
    % relative accuracy to speak of)
    sc = max(abs(mm),abs(mm(2)).^((1:numel(mm))/2));
    m = mm(1:p); e = ee(1:p); bad = find(~(e <= tol*sc(1:p)));
    if ~isempty(bad)
        if isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified)
            o2 = rmf(opts,{'method','tol','verified'}); o2.method = 'engine'; m = ltratio_mom(p,cas,c,o2);
        else
            warning('ltratio:inversion','ltratio_mom: the error estimate of the moments from jets exceeds the tolerance %g (relative) for the moments %s; opts.verified = true recomputes them by the engines.',tol,strjoin(arrayfun(@num2str,bad,'UniformOutput',false),', '));
        end
    end
    info = struct('err',e,'engine','jets','time',toc(t0));
    return
end
opts = oc;
[cas,c,c0,g0] = ltratio_args(cas,c,opts);
if strcmp(engine,'auto') && ltratio_route('plain',opts) && ~strcmp(cas,'DG') && ltratio_engine('mom',ltratio_spec(cas,opts),c,c0,0)
    S = ltratio_spec(cas,opts);
    [m,err] = dfdet_mom(p,S.D,c,ltratio_route('gopts',S,opts)); info = struct('err',err,'engine','general');
    return
end
if strcmp(cas,'DG'), error('ltratio:args','ltratio_mom: the moments are not available for case DG.'); end
if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 1), error('ltratio:args','ltratio_mom: p must be a positive integer.'); end
o = opts; if isfield(o,'method'), o = rmfield(o,'method'); end
o.c0 = c0; o.g0 = g0;
[m,info] = tr_limmom(p,cas,c,o);
end

function o = rmf(o,f)
if ischar(f), f = {f}; end
for k = 1:numel(f), if isfield(o,f{k}), o = rmfield(o,f{k}); end, end
end
