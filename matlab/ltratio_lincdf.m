function P = ltratio_lincdf(b,a,cas,c,opts)
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
%   grid and c <= 0, to about 1e-10, about 0.02 s per call in case DK, 0.1 s in case DG; otherwise (|a| below about
%   0.3 to 3, depending on c, which includes a = 0, the law of U alone, and c > 0) by a conditional integral: in case
%   DK that of the region {U <= aV + b} (LTRATIO_REGION), in case DG the conditional integral given U (TR_LINCOND),
%   to about 1e-10 and in [0, 1].  c must be at least -550; in case DG at most about 3.5, and for a = NaN about 7
%   (an error beyond).
%
%   opts.engine: 'classical' the programs above; 'general' the general engine
%   (DFDET_LIN_CDF, DFDET_V_CDF); 'auto' (the default) the engine of
%   LTRATIO_ENGINE.  cas may be a design of the general engine (LTRATIO_SPEC;
%   with opts.start, opts.g0, opts.gls): the cases DK and DG by the classical
%   programs, any other design by the general engine.
%
%   Example: ltratio_lincdf(-0.2639,-3.5,'DK',0)   % 0.05: the size of the point-optimal test at cbar = -7
%
%   See also POWER_ENVELOPE, LTRATIO_REGION, TR_GILPELAEZ, DFDET_LIN_CDF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts); ltratio_route('nostart','ltratio_lincdf',cas,opts);
[engine,oc] = ltratio_route('engine',opts);
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    % (the classical linear form has the zero start only: no deterministic terms with a fixed start, case 'DK' too, take
    % the general engine, as in the C and Julia routers)
    if ~strcmp(engine,'general') && any(strcmp(S.case,{'DK','DG'})) && strcmp(S.start,'zero')
        oc.engine = engine; P = ltratio_lincdf(b,a,S.case,c,oc); return
    end
    if strcmp(engine,'classical'), error('ltratio:args','ltratio_lincdf: engine = ''classical'': the linear form of the classical programs has the cases DK and DG only, with a zero start.'); end
    P = dfdet_lin_cdf(b,a,S.D,c,ltratio_route('gopts',S,oc)); return
end
if ischar(cas) || isstring(cas)
    cb = -13.5; if isfield(oc,'cbar') && ~isempty(oc.cbar), cb = oc.cbar; end
    if strcmpi(cas,'DG'), S = struct('D',dfdet_trend(1),'start','zero','g0',0,'gls',cb,'case','DG');
    elseif any(strcmpi(cas,{'DK','DU'})), S = struct('D',dfdet_trend(-1),'start','zero','g0',0,'gls',[],'case','DK');
    else, error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
    if ~strcmp(engine,'general'), checks(b,a,cas,c); end     % (the checks of the classical programs come first)
    if strcmp(engine,'general') || (strcmp(engine,'auto') && ltratio_route('plain',rmf(oc,'cbar')) && isnumeric(a) && isscalar(a) && isreal(a) ...
            && isnumeric(c) && isscalar(c) && ltratio_engine('lincdf',S,c,0,a))
        P = dfdet_lin_cdf(b,a,S.D,c,struct('gls',S.gls)); return
    end
end
opts = oc;
P = tr_gilpelaez(b,a,upper(char(cas)),c,opts);
end

function checks(b,a,cas,c)
if ~(ischar(cas) || isstring(cas)) || ~any(strcmpi(cas,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_lincdf: c must be a finite real scalar.'); end
if c < -550, error('ltratio:args','ltratio_lincdf: c must be at least -550 (the range in which the programs are checked; the moments fail from about c = -600, the line integral from about c = -650).'); end
if ~(isnumeric(a) && isscalar(a) && isreal(a)), error('ltratio:args','ltratio_lincdf: a must be a real scalar (NaN: the distribution of V).'); end
if ~(isnumeric(b) && isreal(b)), error('ltratio:args','ltratio_lincdf: b must be real.'); end
end

function o = rmf(o,f)
if isfield(o,f), o = rmfield(o,f); end
end
