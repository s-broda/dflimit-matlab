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
if ltratio_croute('use',cas,opts), [m,info] = ltratio_croute('mom',p,cas,c,opts); return, end   % (the C router)
[engine,oc] = ltratio_route('engine',opts);
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && (isempty(S.case) || strcmp(S.case,'DG')), error('ltratio:args','ltratio_mom: engine = ''classical'': the classical programs have the moments of the cases DK, DU, DT, SK, SU, ST only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && ~strcmp(S.case,'DG') && (strcmp(engine,'classical') || ltratio_route('valid',S,c,oc))
        oc.g0 = S.g0; oc.engine = engine;
        [m,info] = ltratio_mom(p,S.case,c,oc);
        return
    end
    [m,err] = dfdet_mom(p,S.D,c,ltratio_route('gopts',S,oc)); info = struct('err',err,'engine','general');
    return
end
if strcmp(engine,'general')
    S = ltratio_spec(cas,oc);
    [m,err] = dfdet_mom(p,S.D,c,ltratio_route('gopts',S,oc)); info = struct('err',err,'engine','general');
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
