function [q,info] = adfmax_icdf(p,cas,c,opts)
%ADFMAX_ICDF  Quantiles of the limiting ADFmax statistic of Leybourne (1995).
%
%   q = ADFMAX_ICDF(p,cas,c,opts): the p-quantiles (p an array in (0, 1)) of the limit of ADFMAX_CDF, by Newton's
%   method with ADFMAX_PDF as the derivative, safeguarded by bisection on a bracket grown from [x0 - 1, x0 + 1]
%   (x0 = opts.x0, default -2) to opts.tol (1e-10), as adfmax_icdf of julia/DFExact.jl.  info: F at q.
%
%   cas may be a design of the general engine (LTRATIO_SPEC; a design other than
%   the cases DU, DT, SU, ST: DFDET_ADFMAX_ICDF).
%
%   Example: adfmax_icdf([0.01 0.05 0.10],'DT')   % -3.67949 -3.11891 -2.82890 (the quantiles of the ADFmax statistic in the paper)
%
%   See also ADFMAX_CDF, ADFMAX_PDF, DFDET_ADFMAX_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if isfield(opts,'c0') && ~isempty(opts.c0) && opts.c0 ~= 0, error('ltratio:args','adfmax_icdf: ADFmax is the statistic of the null value phi = 1 (opts.c0 does not apply).'); end
fixed = isfield(opts,'g0') && ~isempty(opts.g0) && opts.g0 ~= 0;
if fixed && (ischar(cas) || isstring(cas)) && strncmpi(char(cas),'S',1), error('ltratio:args','adfmax: g0 (a fixed start) and a stationary start (cases SU, ST) exclude each other.'); end
if ~(ischar(cas) || isstring(cas)) || fixed
    % (a design, or a case string with a fixed start, which the classical ADFmax does not have: the general engine)
    [engine,oc] = ltratio_route('engine',opts); S = ltratio_spec(cas,opts);
    if ~isempty(S.gls), error('ltratio:args','adfmax_icdf: ADFmax has least-squares detrending only (no opts.gls).'); end
    if ~strcmp(engine,'general') && any(strcmp(S.case,{'DU','DT','SU','ST'})) && ~strcmp(S.start,'fixed')
        oc.engine = engine; [q,info] = adfmax_icdf(p,S.case,c,oc); return
    end
    if strcmp(engine,'classical'), error('ltratio:args','adfmax_icdf: engine = ''classical'': ADFmax of the classical programs has the cases DU, DT, SU, ST only, with a zero or a stationary start.'); end
    go = ltratio_route('gopts',S,oc); if isfield(oc,'tol'), go.tol = oc.tol; end; if isfield(oc,'x0'), go.x0 = oc.x0; end
    q = dfdet_adfmax_icdf(p,S.D,c,go); info = struct('F',NaN(size(p)),'engine','general'); return
end
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','adfmax_icdf: p must lie in (0,1).'); end
tol = 1e-10; x0 = -2; o = opts;
if isfield(o,'tol'), if ~isempty(o.tol), tol = o.tol; end, o = rmfield(o,'tol'); end
if isfield(o,'x0'), if ~isempty(o.x0), x0 = o.x0; end, o = rmfield(o,'x0'); end
q = zeros(size(p)); Fq = q;
for k = 1:numel(p)
    pk = p(k); x = x0; lo = x - 1; hi = x + 1;
    while adfmax_cdf(lo,cas,c,o) > pk, lo = lo - 1; end
    while adfmax_cdf(hi,cas,c,o) < pk, hi = hi + 1; end
    x = (lo + hi)/2; conv = false;
    for it = 1:100
        F = adfmax_cdf(x,cas,c,o); f = adfmax_pdf(x,cas,c,o);
        if F > pk, hi = x; else, lo = x; end
        xn = x - (F - pk)/f;
        if ~(lo <= xn && xn <= hi) || ~isfinite(xn), xn = (lo + hi)/2; end   % (a step below half an ulp lands on the end of the bracket: converged)
        if abs(xn - x) < tol, x = xn; conv = true; break, end
        x = xn;
    end
    if ~conv, warning('ltratio:icdf','adfmax_icdf: no convergence at p = %g.',pk); end
    q(k) = x; Fq(k) = adfmax_cdf(x,cas,c,o);
end
info = struct('F',Fq);
end
