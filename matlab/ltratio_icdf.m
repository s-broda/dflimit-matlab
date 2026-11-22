function [q,info] = ltratio_icdf(p,cas,c,opts)
%LTRATIO_ICDF  Quantiles of the limiting Dickey-Fuller t-ratio.
%
%   q = LTRATIO_ICDF(p,cas,c,opts): the limit of the p-quantile of the
%   t-ratio, p in (0,1) (an array), by Newton's method on LTRATIO_CDF with
%   LTRATIO_PDF as the derivative, safeguarded by bisection, started at the
%   normal approximation from the first two moments (LTRATIO_MOM).  Arguments
%   and options as in LTRATIO_CDF; opts.tol (1e-10) is the tolerance on q.
%   info: F, the distribution function at q; iter; conv (logical).
%
%   opts.engine: 'classical' the search above; 'general' the search of the
%   general engine (DFDET_ICDF, which evaluates F and f in one pass on one
%   path); 'auto' (the default) the engine of LTRATIO_ENGINE (and in the
%   classical search each evaluation routed by LTRATIO_CDF).  cas may be a
%   design of the general engine (LTRATIO_SPEC; any design that is not a
%   classical case goes to DFDET_ICDF).
%
%   Example: ltratio_icdf([0.01 0.05 0.10], 'DU')   % the asymptotic critical values of tau_mu
%
%   See also LTRATIO_CDF, LTRATIO_PDF, LTRATIO_MOM, DFDET_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts,{'tol'}) || ltratio_croute('usemethod',cas,opts,{'tol'})                % (the C router; info, the distribution function at q, only if asked for)
    if nargout > 1, [q,info] = ltratio_croute('icdf',p,cas,c,opts); else, q = ltratio_croute('icdf',p,cas,c,opts); end
    return
end
[engine,oc] = ltratio_route('engine',opts);
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && isempty(S.case), error('ltratio:args','ltratio_icdf: engine = ''classical'': the classical programs have the cases DK, DU, DT, SK, SU, ST and DG only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && (strcmp(engine,'classical') || ltratio_route('valid',S,c,oc))
        oc.g0 = S.g0; oc.engine = engine;
        [q,info] = ltratio_icdf(p,S.case,c,oc);
        return
    end
    ltratio_route('nomethod',opts);
    [q,info] = general_icdf(p,S,c,oc);
    return
end
if strcmp(engine,'general')
    ltratio_route('nomethod',opts);
    [q,info] = general_icdf(p,ltratio_spec(cas,oc),c,oc);
    return
end
[cas,c,c0] = ltratio_args(cas,c,oc);
if strcmp(engine,'auto') && ltratio_route('plain',oc,{'tol'})
    S = ltratio_spec(cas,oc);
    if ltratio_engine('icdf',S,c,c0,0)
        if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','ltratio_icdf: p must lie in (0,1).'); end
        [q,info] = general_icdf(p,S,c,oc);
        return
    end
end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','ltratio_icdf: p must lie in (0,1).'); end
tol = 1e-10; if isfield(opts,'tol') && ~isempty(opts.tol), tol = opts.tol; o = rmfield(opts,'tol'); else, o = opts; end
if strcmp(upper(char(cas)),'DG')
    mu = -1.9 + 0.1*min(c,0); sd = 0.8;                % (no moments for case DG: a rough start, the search is safeguarded)
else
    m = ltratio_mom(2,cas,c,o); mu = m(1); sd = sqrt(max(m(2) - m(1)^2,1e-6));
end
q = zeros(size(p)); Fq = q; it = q; conv = false(size(p));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (a brittle regime warned once, not at every step)
for k = 1:numel(p)
    pk = p(k);
    x = mu + sd*sqrt(2)*erfinv(2*pk - 1);
    [q(k),Fq(k),it(k),conv(k)] = ltratio_newton(@(y) ltratio_cdf(y,cas,c,o), @(y) ltratio_pdf(y,cas,c,o), x, sd, tol, pk);
end
info = struct('F',Fq,'iter',it,'conv',conv);
if ~all(conv(:)), warning('ltratio:icdf','ltratio_icdf: Newton''s method did not converge at %d of %d points.',nnz(~conv),numel(p)); end
end

function [q,info] = general_icdf(p,S,c,oc)
% the quantiles of the general engine, started at the normal approximation from the moments of the engine the rule takes
% for them (the classical moments where it takes the classical engine)
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','ltratio_icdf: p must lie in (0,1).'); end
o = ltratio_route('gopts',S,oc);
if isfield(oc,'tol') && ~isempty(oc.tol), o.tol = oc.tol; end
c0 = 0; if isfield(oc,'c0') && ~isempty(oc.c0), c0 = oc.c0; end
m = [];
if ~isempty(S.case) && ~strcmp(S.case,'DG') && ltratio_route('valid',S,c,oc) && ~ltratio_engine('mom',S,c,c0,0)
    m = ltratio_mom(2,S.case,c,struct('c0',c0,'g0',S.g0,'engine','classical'));
end
q = zeros(size(p));
for k = 1:numel(p)
    ok = o;
    if ~isempty(m), ok.x0 = m(1) + sqrt(max(m(2) - m(1)^2,1e-6))*sqrt(2)*erfinv(2*p(k) - 1); end
    q(k) = dfdet_icdf(p(k),S.D,c,ok);
end
o = rmfield(o,intersect(fieldnames(o),{'tol','x0'}));
Fq = dfdet_cdf(q,S.D,c,o);                              % (the distribution function at the quantiles, as the classical search)
info = struct('F',Fq,'iter',NaN(size(p)),'conv',abs(Fq - p) < 1e-6,'engine','general');
end
