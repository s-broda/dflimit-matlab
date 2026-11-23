function q = dfdet_icdf(p,D,c,opts)
%DFDET_ICDF  Quantiles of the limiting t-ratio for any deterministic terms (the general engine).
%
%   q = DFDET_ICDF(p,D,c,opts): the limit of the p-quantile of the t-ratio,
%   p in (0,1) (an array), by Newton's method on F(x) - p with the density
%   from the same evaluation of the contour integral (DFDET_CDF, DFDET_PDF on
%   one path), each step limited to 1, at most 30 steps, and where that does
%   not converge the safeguarded Newton's method (a bracket, bisection where a
%   step leaves it).  The start is opts.x0 or the normal approximation from
%   the first two moments (DFDET_MOM; at c0 = 0, where they are cheap, when
%   c0 ~= 0; -2 for the statistic that detrends first).  opts.tol (1e-10) is
%   the tolerance on q; the other options as in DFDET_CDF.
%
%   Examples
%      dfdet_icdf([0.01 0.05 0.10], dfdet_perron('C', 0.3))   % critical values with a break at 0.3
%      dfdet_icdf(0.05, dfdet_detrend_first(dfdet_perron('B', 0.5)))
%      dfdet_icdf(0.05, dfdet_trend(1), 0, struct('gls', -13.5))   % the DF-GLS test
%
%   See also DFDET_CDF, DFDET_PDF, LTRATIO_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[c,c0,start,g0,gls,o] = dfdet_args(D,c,opts,'the quantiles of the t-ratio (dfdet_icdf)');
dfdet_method(opts,'dfdet_icdf');
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','dfdet_icdf: p must lie in (0,1).'); end
tol = getf(opts,'tol',1e-10);
x0 = getf(opts,'x0',[]); mu = -2; sd = 1;
if isempty(x0) && ~strcmp(D.detrend,'first')
    try
        mo = opts; mo.c0 = 0; m = dfdet_mom(2,D,c,mo);
        mu = m(1); sd = sqrt(max(m(2) - m(1)^2,1e-6));
    catch ME
        if ~strncmp(ME.identifier,'dfdet:',6), rethrow(ME); end
    end
end
Mo = dfdet_startmodel(D,c,start,g0,gls);
CP = dfdet_cpath(Mo,c,struct('phi',getf(o,'phi',[]),'hfac',getf(o,'hfac',1)));
co = struct('c0',c0);
for f = {'n','ntail','alpha','filon','warn'}, if isfield(o,f{1}), co.(f{1}) = o.(f{1}); end, end
FF = @(x) ev(CP,x,co);
q = zeros(size(p));
for k = 1:numel(p)
    if isempty(x0), xs = mu + sd*sqrt(2)*erfinv(2*p(k) - 1); else, xs = x0; end
    q(k) = dfdet_newton(FF,p(k),xs,1,tol,100,30);
end
end

function v = ev(CP,x,co)
CPr = dfdet_contour('refine',CP,x);
[F,f] = dfdet_contour('tratio',CPr,x,co);
v = [F f];
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
