function F = coef_cdf(y,cas,c,opts)
%COEF_CDF  Limiting distribution function of the Dickey-Fuller coefficient statistic.
%
%   F = COEF_CDF(y,cas,c,opts) is the limit as n -> Inf of P[n(phi^ - phi0) <= y] (y an array) for the least squares
%   estimate phi^ of phi in the regression of LTRATIO_CDF, under phi = 1 + c/n, with phi0 = 1 + c0/n (opts.c0,
%   default 0: the coefficient test of Dickey and Fuller, rho^ = n(phi^ - 1)).  The limit is (U - c0 V)/V, so that
%   P[(U - c0 V)/V <= y] = P[U - (c0 + y)V <= 0], the distribution function of the t-ratio at 0 with the local null
%   value c0 + y: LTRATIO_CDF(0,cas,c,opts) with c0 + y in place of c0, in every case and with every start and
%   option of LTRATIO_CDF (method, g0).
%
%   opts.engine is passed to LTRATIO_CDF (each value routed there); cas may be
%   a design of the general engine (LTRATIO_SPEC; a design that is not a
%   classical case: DFDET_COEF_CDF).
%
%   Example: coef_cdf(-14.094,'DU',0)    % 0.05: the 5% critical value of rho^_mu
%
%   See also COEF_ICDF, LTRATIO_CDF, DFDET_COEF_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
ltratio_route('case',cas,opts);
if ltratio_croute('use',cas,opts), F = ltratio_croute('coef',y,cas,c,opts); return, end   % (the C router)
if ~(ischar(cas) || isstring(cas))
    [engine,oc] = ltratio_route('engine',opts); S = ltratio_spec(cas,opts);
    if strcmp(engine,'classical') && isempty(S.case), error('ltratio:args','coef_cdf: engine = ''classical'': the classical programs have the cases DK, DU, DT, SK, SU, ST and DG only.'); end
    if ~strcmp(engine,'general') && ~isempty(S.case) && (strcmp(engine,'classical') || ltratio_route('valid',S,c,oc))
        oc.g0 = S.g0; oc.engine = engine; F = coef_cdf(y,S.case,c,oc); return
    end
    F = dfdet_coef_cdf(y,S.D,c,ltratio_route('gopts',S,oc)); return
end
if ~(isnumeric(y) && isreal(y)), error('ltratio:args','coef_cdf: y must be real.'); end
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
F = zeros(size(y));
ltratio_brittle('begin'); cl = onCleanup(@() ltratio_brittle('end'));   % (the regimes c0 + y passes, each warned once)
for k = 1:numel(y)
    if isnan(y(k)), F(k) = NaN; elseif y(k) == Inf, F(k) = 1; elseif y(k) == -Inf, F(k) = 0;
    else, o = opts; o.c0 = c0 + y(k); F(k) = ltratio_dist(0,cas,c,o,'cdf'); end
end
% (the final check of the values, named as this call)
a = sprintf('case %s, c = %g',upper(char(cas)),c); if c0 ~= 0, a = [a sprintf(', c0 = %g',c0)]; end
if isfield(opts,'g0') && ~isempty(opts.g0) && opts.g0 ~= 0, a = [a sprintf(', g0 = %g',opts.g0)]; end
F = ltratio_rangecheck('coef_cdf',a,'y',y,F);
end
