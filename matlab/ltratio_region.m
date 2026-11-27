function [P,info] = ltratio_region(bnd,cas,c,opts)
%LTRATIO_REGION  The limiting probability of a region {U <= g(V)} of the plane.
%
%   P = LTRATIO_REGION(bnd,cas,c,opts) is P{U <= g(V)} under phi = 1 + c/n, no drift, zero start, case 'DK' (no
%   deterministic terms; also the constant case with GLS demeaning), with U = int J dJ and V = int J^2 as in
%   LTRATIO_CDF.  Every test of the unit root in the limit experiment of case DK rejects on such a region (the companion
%   paper's section on the weighted-average-power test), and the first part of the conditional integral representation, the
%   probability of a region of the plane, gives its probability for every boundary g, by one outer integral over V
%   (TR_REGION).
%   bnd is a structure with the field type and its parameters:
%      'line'      g(v) = b + a v           (fields a, b: the point-optimal test, U - aV <= b)
%      'parabola'  g(v) = x sqrt(v) + c0 v  (fields x, c0: the t-ratio, P = LTRATIO_CDF(x,'DK',c) with c0)
%      'wap'       the weighted-average-power test with the half-normal weight of scale omega on c < 0,
%                  (1 + omega^2 V)^(-1/2) erfcx(U/sqrt(2 (V + 1/omega^2))) >= k   (fields k, omega)
%      'fun'       g a function handle, elementwise in v (field g; the MATLAB engine only)
%   and optionally v0 (g = +Inf on (0, v0), whose probability P(V < v0) is added) and brks (kinks of g, as
%   breakpoints of the outer rule).  Accurate to about 1e-11; 0.01 to 0.03 s per probability (one thread, November 2026).
%
%   Example: ltratio_region(struct('type','wap','k',3.6658,'omega',10),'DK',0)   % 0.05 (the weighted-average-power test of the companion paper)
%
%   See also LTRATIO_LINCDF, POWER_ENVELOPE, LTRATIO_CDF, TR_REGION.
if nargin < 2 || isempty(cas), cas = 'DK'; end
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
% (opts.engine: 'general' the general engine, DFDET_REGION, for the lines and parabolas; 'auto' and 'classical' the
% programs above; a design other than case DK goes to DFDET_REGION)
ltratio_route('case',cas,opts); ltratio_route('nostart','ltratio_region',cas,opts);
if ltratio_croute('use',cas,opts) && ltratio_croute('bnd',bnd), [P,info] = ltratio_croute('region',bnd,cas,c,opts); return, end   % (the C router)
[engine,oc] = ltratio_route('engine',opts);
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    % (the classical regions have the zero start only: no deterministic terms with a fixed start, whose case is also 'DK',
    % take the general engine, as in the C and Julia routers)
    if ~strcmp(engine,'general') && strcmp(S.case,'DK') && strcmp(S.start,'zero'), oc = rmfield_(oc,'g0'); [P,info] = ltratio_region(bnd,'DK',c,oc); return, end
    if strcmp(engine,'classical'), error('ltratio:args','ltratio_region: engine = ''classical'': the regions of the classical programs have case DK only, with a zero start.'); end
    P = dfdet_region(bnd,S.D,c,ltratio_route('gopts',S,oc)); info = struct('engine','general'); return
end
if strcmp(engine,'general'), P = dfdet_region(bnd,dfdet_trend(-1),c); info = struct('engine','general'); return, end
opts = oc;
if ~(ischar(cas) || isstring(cas)) || ~strcmpi(cas,'DK'), error('ltratio:args','ltratio_region: the case must be DK (zero start).'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_region: c must be a finite real scalar.'); end
if c < -550, error('ltratio:args','ltratio_region: c must be at least -550 (the range in which the programs are checked).'); end
if ~isstruct(bnd) || ~isfield(bnd,'type'), error('ltratio:args','ltratio_region: bnd must be a structure with the field type.'); end
% (the line and the parabola at a non-finite b or x: 1, 0, NaN, as the C router gives them; the classical programs gave 0,
% also at NaN; the slope a and c0 must be finite)
ty = lower(char(bnd.type));
if any(strcmp(ty,{'line','parabola'})) && ~any(isfield(bnd,{'v0','brks','g'}))
    if strcmp(ty,'line'), s = bnd.a; v = bnd.b; else, s = 0; if isfield(bnd,'c0') && ~isempty(bnd.c0), s = bnd.c0; end, v = bnd.x; end
    if ~(isnumeric(s) && isscalar(s) && isreal(s) && isfinite(s))
        if strcmp(ty,'line'), error('ltratio:args','ltratio_region: the slope a of the line must be a finite real scalar.'); end
        error('ltratio:args','ltratio_region: c0 of the parabola must be a finite real scalar.');
    end
    if isnumeric(v) && isscalar(v) && ~isfinite(v)
        if isnan(v), P = NaN; elseif v > 0, P = 1; else, P = 0; end
        info = struct('engine','classical','err',0,'status',0); return
    end
end
[P,info] = tr_region(bnd,c,opts);
end

function o = rmfield_(o,f)
if isfield(o,f), o = rmfield(o,f); end
end
