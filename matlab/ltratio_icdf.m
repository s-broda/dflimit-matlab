function [q,info] = ltratio_icdf(p,cas,c,opts)
%LTRATIO_ICDF  Quantiles of the limiting Dickey-Fuller t-ratio.
%
%   q = LTRATIO_ICDF(p,cas,c,opts): the limit of the p-quantile of the
%   t-ratio, p in (0,1) (an array), by Newton's method on LTRATIO_CDF with
%   LTRATIO_PDF as the derivative, safeguarded by bisection, started at the
%   normal approximation from the first two moments: from jets, the device of
%   Yong Bao and Raymond Kan (KAN_CORE), where they hold (c <= 20, an error
%   estimate within 1e-6 relative; mom_start of DFLimit.jl), otherwise from
%   LTRATIO_MOM.  Where the rule gives a Newton step's point to the general
%   engine (opts.engine 'auto' without further options), F and f come from one
%   pass of its paths (route_dist2 of DFLimit.jl), the values of LTRATIO_CDF
%   and LTRATIO_PDF bit for bit.  LTRATIO_SWITCH('mom_start',false) and
%   LTRATIO_SWITCH('dist2',false) restore the computation of before.  Arguments
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
[cas,c,c0,g0] = ltratio_args(cas,c,oc);
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
    % (the start: the moments from jets of Bao and Kan where they hold, 0.2 to 1 ms in Julia; otherwise LTRATIO_MOM)
    m = mom_start(2,cas,c,c0,g0);
    if isempty(m), m = ltratio_mom(2,cas,c,o); end
    mu = m(1); sd = sqrt(max(m(2) - m(1)^2,1e-6));
end
% (F and f in one evaluation where the rule gives the point to the general engine's paths: route_dist2)
Fb = [];
if strcmp(engine,'auto') && ltratio_route('plain',oc,{'tol'}) && ltratio_switch('dist2'), S = ltratio_spec(cas,oc); Fb = @(y) dist2(y,S,cas,c,c0,g0,oc); end
q = zeros(size(p)); Fq = q; it = q; conv = false(size(p));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (a brittle regime warned once, not at every step)
for k = 1:numel(p)
    pk = p(k);
    x = mu + sd*sqrt(2)*erfinv(2*pk - 1);
    [q(k),Fq(k),it(k),conv(k)] = ltratio_newton(@(y) ltratio_cdf(y,cas,c,o), @(y) ltratio_pdf(y,cas,c,o), x, sd, tol, pk, Fb);
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
% (the six classical cases: the moments from jets, mom_start, where they hold; otherwise the engine of the rule)
cv = ~isempty(S.case) && ~strcmp(S.case,'DG') && ltratio_route('valid',S,c,oc);
if cv, m = mom_start(2,S.case,c,c0,S.g0); end
if isempty(m) && cv && ~ltratio_engine('mom',S,c,c0,0)
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

function m = mom_start(p,cas,c,c0,g0)
% The first p moments of the t-ratio of a classical case for the start of a quantile search (mom_start of DFLimit.jl):
% from one jet per node, the method of Bao and Kan (KAN_CORE 'ratio_mom', where the circles of the engines take 10 to
% 150 ms in Julia), or [] outside its region (c > 20; case DG), where it fails or where its error estimate exceeds 1e-6
% relative (to max(|m_j|, m_2^(j/2))); the search then starts from the engines' circles.  The start only places the
% first step and the bracket of the search: the quantile is the root of F(x) - p to the tolerance of the search either way
m = [];
cas = upper(char(cas));
if ~ltratio_switch('mom_start') || ~any(strcmp(cas,{'DK','DU','DT','SK','SU','ST'})) || ~(c <= 20), return, end
[model,stat] = kan_route('model',cas);
try
    [mm,ee] = kan_core('ratio_mom',kan_core('limcf',model,stat,c,g0),max(p,2),c0 - c,0.5,1e-10);
catch ME
    if strcmp(ME.identifier,'MATLAB:nomem'), rethrow(ME); end
    return
end
sc = max(abs(mm),abs(mm(2)).^((1:numel(mm))/2));
if all(isfinite(mm(1:p)) & ee(1:p) <= 1e-6*sc(1:p)), m = mm(1:p); end
end

function Ff = dist2(x,S,cas,c,c0,g0,oc)
% F and f at a point x of the Newton steps in one evaluation (route_dist2 of DFLimit.jl): where the rule gives the point
% to the general engine, its paths return both from one pass (DFDET_DIST, mode 'both'), with the values of LTRATIO_CDF
% and LTRATIO_PDF bit for bit; [] elsewhere (the classical programs, the tails that go back to them), and the search takes
% the two evaluations
Ff = [];
if ~isfinite(x) || ~ltratio_engine('dist',S,c,c0,x), return, end
o = struct('start',S.start,'g0',S.g0,'gls',S.gls);
for f = {'c0','verified','vtol','phi','hfac','n','ntail','alpha','filon','warn','tol','method'}
    if isfield(oc,f{1}), o.(f{1}) = oc.(f{1}); end
end
o = rmfield_(o,'tol');                                    % (the tolerance of the search, not of the engine)
[F,f] = dfdet_dist(x,S.D,c,o,'both');
if f < 0, f = 0; end
% (the tails below 1e-3 where the classical programs would take the conditional integral go back to them, as in
% LTRATIO_DIST, for F and f each: the search then takes the two evaluations)
if (F < 1e-3 || f < 1e-3) && ltratio_method(cas,c,c0,g0,x,false), return, end
a = sprintf('case %s, c = %g',cas,c);
if c0 ~= 0, a = [a sprintf(', c0 = %g',c0)]; end
if g0 ~= 0, a = [a sprintf(', g0 = %g',g0)]; end
F = ltratio_rangecheck('ltratio_cdf',a,'x',x,F);
Ff = [F f];
end

function o = rmfield_(o,f)
if isfield(o,f), o = rmfield(o,f); end
end
