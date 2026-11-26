function [P,err] = ltratio_lincdf(b,a,cas,c,opts)
%LTRATIO_LINCDF  Distribution function of the linear form U - aV of the limiting functionals.
%
%   P = LTRATIO_LINCDF(b,a,cas,c,opts) is the limit of P(U - aV <= b) under phi = 1 + c/n (b an array), zero start,
%   no drift, with U = int J dJ and V = int J^2 as in LTRATIO_CDF: case 'DK' (no deterministic terms, and the
%   constant case with GLS demeaning, whose limit is the same) or 'DG' (the GLS-detrended trend case of Elliott,
%   Rothenberg and Stock, detrended at opts.cbar, default -13.5).  The point-optimal test of Elliott, Rothenberg
%   and Stock against cbar < 0 rejects for U - aV < b: in case 'DK' with a = cbar/2 (its log likelihood ratio is
%   cbar U - cbar^2 V/2), in case 'DG', with the data detrended at cbar, with a = cbar^2/(2 (cbar - 1)) (it rejects
%   for small cbar^2 V + (1 - cbar) J(1)^2, J(1)^2 = 2U + 1); see POWER_ENVELOPE.  a = NaN gives P(V <= b) (case
%   DK).  By the inversion formula of Gil-Pelaez (TR_GILPELAEZ) where its transform has decayed by the end of its
%   grid and c <= 0, to about 1e-10, about 0.008 s per call in case DK, 0.02 to 0.04 s in case DG (one thread, November
%   2026; the default, the single inversion of Bao and Kan, 0.004 to 0.007 s in case DK); otherwise (|a| below about
%   0.3 to 3, depending on c, which includes a = 0, the law of U alone, and c > 0) by a conditional integral: in case
%   DK that of the region {U <= aV + b} (LTRATIO_REGION), in case DG the conditional integral given U (TR_LINCOND),
%   to about 1e-10 and in [0, 1].  c must be at least -550; in case DG at most about 3.5, and for a = NaN about 7
%   (an error beyond).
%
%   In case DK the default (opts.method 'auto' with opts.engine 'auto', a ~= 0, -100 <= c <= 4) is the single
%   inversion of the joint characteristic function along the ray of the half-plane, the method of Yong Bao and
%   Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in First-Order
%   Autoregressions: O(1) Algorithms" (KAN_CORE), with the tail beyond S = 400/|w| (w = 1/2 + b) completed by
%   integration by parts; opts.tol (1e-12) its absolute tolerance, [P,err] = LTRATIO_LINCDF(...) returns the
%   estimates, and the warning ltratio:inversion names the points at which they exceed it.  The six classical cases
%   through their designs (dfdet_trend(-1), dfdet_trend(0), dfdet_trend(1) with a zero, fixed or stationary start) take
%   it too for -100 <= c <= 2.  opts.method = 'inversion' forces it; opts.method = 'engine' (and opts.engine
%   'classical' or 'general', case DG, a = 0 and a = NaN, options of the engines) is the computation of before above
%   (err NaN); opts.verified = true recomputes the values whose estimate exceeds the tolerance by it.
%
%   opts.engine: 'classical' the programs above; 'general' the general engine
%   (DFDET_LIN_CDF, DFDET_V_CDF); 'auto' (the default) the engine of
%   LTRATIO_ENGINE.  cas may be a design of the general engine (LTRATIO_SPEC;
%   with opts.start, opts.g0, opts.gls): the cases DK and DG with a zero start by
%   the classical programs, any other design by the general engine.
%
%   Example: ltratio_lincdf(-0.2639,-3.5,'DK',0)   % 0.05: the size of the point-optimal test at cbar = -7
%
%   See also POWER_ENVELOPE, LTRATIO_REGION, TR_GILPELAEZ, DFDET_LIN_CDF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts); ltratio_route('nostart','ltratio_lincdf',cas,opts);
if ltratio_croute('use',cas,opts,{'cbar','tol'},true), [P,err] = ltratio_croute('lincdf',b,a,cas,c,opts); return, end   % (the C router)
m = kan_route('method',opts,'ltratio_lincdf',false);
[engine,oc] = ltratio_route('engine',opts);
oe = rmf(rmf(oc,'method'),'tol');                   % (the options of the engines' computation)
err = NaN(size(b));
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    % the single inversion of Bao and Kan for the six classical cases with any start, within the classical checks of the
    % arguments; opts.method = 'inversion' forces it
    inv = ~strcmp(engine,'general') && any(strcmp(S.case,{'DK','DU','DT','SK','SU','ST'})) && isnumeric(a) && isscalar(a) && isfinite(a) ...
          && cvalid(S,c);
    if strcmp(m,'inversion') && ~inv, error('ltratio:args','ltratio_lincdf: opts.method = ''inversion'' applies to the cases DK, DU, DT, SK, SU, ST (dfdet_trend(-1), dfdet_trend(0), dfdet_trend(1)) within the classical checks, a finite a and engine ''auto'' or ''classical''.'); end
    if inv
        [P,err,done] = lin_inv(b,a,S,c,engine,m,opts,oe);
        if done, return, end
    end
    % (the classical linear form has the zero start only: no deterministic terms with a fixed start, case 'DK' too, take
    % the general engine, as in the C and Julia routers)
    if ~strcmp(engine,'general') && any(strcmp(S.case,{'DK','DG'})) && strcmp(S.start,'zero')
        oe.engine = engine; oe.method = 'engine'; P = ltratio_lincdf(b,a,S.case,c,oe); return
    end
    if strcmp(engine,'classical'), error('ltratio:args','ltratio_lincdf: engine = ''classical'': the linear form of the classical programs has the cases DK and DG only, with a zero start.'); end
    P = dfdet_lin_cdf(b,a,S.D,c,ltratio_route('gopts',S,oe)); return
end
cb = -13.5; if isfield(oc,'cbar') && ~isempty(oc.cbar), cb = oc.cbar; end
% (a non-finite cbar did not return in the classical programs)
if ~(isnumeric(cb) && isscalar(cb) && isreal(cb) && isfinite(cb) && cb < 0), error('ltratio:args','ltratio_lincdf: opts.cbar must be a finite real scalar < 0.'); end
if strcmpi(cas,'DG'), S = struct('D',dfdet_trend(1),'start','zero','g0',0,'gls',cb,'case','DG');
elseif any(strcmpi(cas,{'DK','DU'})), S = struct('D',dfdet_trend(-1),'start','zero','g0',0,'gls',[],'case','DK');
else, error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
if strcmp(engine,'general') && strcmp(m,'inversion'), error('ltratio:args','ltratio_lincdf: opts.method = ''inversion'' (the single inversion of Bao and Kan) does not apply with engine = ''general''.'); end
if ~strcmp(engine,'general'), checks(b,a,cas,c); end     % (the checks of the classical programs come first)
if strcmp(m,'inversion') && isnan(a), error('ltratio:args','ltratio_lincdf: opts.method = ''inversion'' needs a finite a (a = NaN, the law of V, is the general engine''s).'); end
if strcmp(m,'inversion') && strcmp(S.case,'DG'), error('ltratio:args','ltratio_lincdf: opts.method = ''inversion'' applies to case DK (the designs of the six classical cases), not to case DG.'); end
if strcmp(S.case,'DK') && ~isnan(a) && ~strcmp(engine,'general')
    S.model = 1; S.stat = false;
    [P,err,done] = lin_inv(b,a,S,c,engine,m,opts,oe);
    if done, return, end
end
if strcmp(engine,'general') || (strcmp(engine,'auto') && ltratio_route('plain',rmf(oe,'cbar')) && isnumeric(a) && isscalar(a) && isreal(a) ...
        && isnumeric(c) && isscalar(c) && ltratio_engine('lincdf',S,c,0,a))
    P = dfdet_lin_cdf(b,a,S.D,c,struct('gls',S.gls)); return
end
P = tr_gilpelaez(b,a,upper(char(cas)),c,oe);
end

function [P,err,done] = lin_inv(b,a,S,c,engine,m,opts,oe)
% the linear form of a classical case by the single inversion of Bao and Kan where KAN_ROUTE takes it and KAN_CORE's
% region lin_ok holds (all the points with opts.method = 'inversion'; a ~= 0, -100 <= c <= 4 in case DK with a zero start,
% -100 <= c <= 2 in the other cases), the engines at the other points (lin_engine of merged.jl: case DK with a zero start
% by the classical programs, the other cases by the general engine)
P = NaN(size(b)); err = NaN(size(b)); done = false;
if ~kan_route('take','ltratio_lincdf',m,engine,S.case,c,kan_route('plain',rmf(oe,'cbar'))), return, end
[model,stat] = kan_route('model',S.case); g0 = S.g0;
if strcmp(m,'inversion'), ok = true(size(b)); else, ok = repmat(kan_core('lin_ok',model,stat,c,g0,a),size(b)); end
if ~any(ok(:)), return, end
tol = kan_route('tol',opts,1e-12);
K = kan_core('limcf',model,stat,c,g0);
for k = find(ok(:))', [P(k),err(k)] = kan_core('lin_cdf1',K,b(k),a,tol); end
if any(~ok(:)), P(~ok) = lin_engine(b(~ok),a,S,c); end
if isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified)
    bad = ok & ~(err <= tol);
    if any(bad(:)), P(bad) = lin_engine(b(bad),a,S,c); end
else
    kan_route('warn','ltratio_lincdf',err(ok),tol);
end
done = true;
end

function ok = cvalid(S,c)
% whether the classical checks of the C and Julia engines accept the case of S at c (LTRATIO_ARGS without the restriction of
% the MATLAB programs' conditional integral, which the single inversion does not take)
try
    ltratio_args(S.case,c,struct('g0',S.g0),false); ok = true;
catch ME
    if ~strcmp(ME.identifier,'ltratio:args'), rethrow(ME); end
    ok = false;
end
end

function P = lin_engine(b,a,S,c)
if strcmp(S.case,'DK') && strcmp(S.start,'zero'), P = ltratio_lincdf(b,a,'DK',c,struct('method','engine')); return, end
P = dfdet_lin_cdf(b,a,S.D,c,ltratio_route('gopts',S,struct));
end

function checks(b,a,cas,c)
if ~(ischar(cas) || isstring(cas)) || ~any(strcmpi(cas,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_lincdf: c must be a finite real scalar.'); end
if c < -550, error('ltratio:args','ltratio_lincdf: c must be at least -550 (the range in which the programs are checked; the moments fail from about c = -600, the line integral from about c = -650).'); end
if ~(isnumeric(a) && isscalar(a) && isreal(a) && ~isinf(a)), error('ltratio:args','ltratio_lincdf: a must be a finite real scalar (NaN: the distribution of V).'); end
if ~(isnumeric(b) && isreal(b)), error('ltratio:args','ltratio_lincdf: b must be real.'); end
end

function o = rmf(o,f)
if isfield(o,f), o = rmfield(o,f); end
end
