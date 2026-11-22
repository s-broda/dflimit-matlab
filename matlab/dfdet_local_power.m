function P = dfdet_local_power(c,D,opts)
%DFDET_LOCAL_POWER  Asymptotic local power of the left-tailed tests of the general engine.
%
%   P = DFDET_LOCAL_POWER(c,D,opts): the power of the left-tailed level-alpha
%   test (opts.alpha, 0.05) of phi = 1 + c0/n against phi = 1 + c/n for each c in
%   the array c: P_c(stat <= q), q the alpha-quantile under c = c0 with the same
%   start (a stationary start needs c0 < 0), stat the t-ratio (opts.statistic
%   'tratio', the default) or n(phi^ - phi0) ('coef').  Options c0, start, g0,
%   gls as in DFDET_CDF.
%
%   Example: dfdet_local_power([-5 -10 -20], dfdet_perron('C', 0.3))
%
%   See also DFDET_CDF, DFDET_ICDF, DFDET_COEF_CDF.
if nargin < 3 || isempty(opts), opts = struct; end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; o = rmfield(opts,'alpha'); else, o = opts; end
st = 'tratio'; if isfield(o,'statistic') && ~isempty(o.statistic), st = lower(char(o.statistic)); o = rmfield(o,'statistic'); end
c0 = 0; if isfield(o,'c0') && ~isempty(o.c0), c0 = o.c0; end
P = zeros(size(c));
switch st
    case 'tratio'
        q0 = o; if isfield(o,'gls') && ~isempty(o.gls), q0.x0 = -2.5; else, q0.x0 = -3.5; end
        q0.tol = 1e-11;
        q = dfdet_icdf(alpha,D,c0,q0);
        for k = 1:numel(c), P(k) = dfdet_cdf(q,D,c(k),o); end
    case 'coef'
        q0 = o; q0.y0 = -15; q0.tol = 1e-11;
        q = dfdet_coef_icdf(alpha,D,c0,q0);
        for k = 1:numel(c), P(k) = dfdet_coef_cdf(q,D,c(k),o); end
    otherwise
        error('dfdet:args','dfdet_local_power: opts.statistic must be ''tratio'' or ''coef''.');
end
end
