function P = dfdet_power_envelope(c,D,opts)
%DFDET_POWER_ENVELOPE  The asymptotic power envelope of Elliott, Rothenberg and Stock by the general engine.
%
%   P = DFDET_POWER_ENVELOPE(c,D,opts): at each c < 0 (an array) the power of the
%   point-optimal level-alpha test (opts.alpha, 0.05) against c itself, which
%   rejects for small U - aV: without deterministic terms and with an
%   intercept (whose point-optimal test GLS-demeans the data, which leaves the
%   limit of the case without deterministic terms) a = c/2 and the limit of no
%   deterministic terms; with a linear trend the data GLS-detrended at
%   cbar = c and a = c^2/(2 (c - 1)).  The critical value is the alpha-quantile
%   of U - aV at c = 0 (DFDET_LIN_ICDF), the envelope its power at c
%   (DFDET_LIN_CDF); alpha for c >= 0.  D: dfdet_trend(-1), dfdet_trend(0) or
%   dfdet_trend(1); for other designs the envelope of the point-optimal
%   statistic of the design is PT_ENVELOPE.
%
%   Example: dfdet_power_envelope(-7, dfdet_trend(-1))   % 0.4994
%
%   See also PT_ENVELOPE, DFDET_LIN_CDF, POWER_ENVELOPE.
if nargin < 3 || isempty(opts), opts = struct; end
alpha = 0.05; if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
if ~(isnumeric(alpha) && isscalar(alpha) && alpha > 0 && alpha < 1), error('dfdet:args','dfdet_power_envelope: the level opts.alpha must lie in (0, 1).'); end
K = size(D.C{1},2);
trendcase = K == 2 && isequal(D.br,[0 1]) && isequal(D.C{1},eye(2)) && strcmp(D.detrend,'onestep') && ~D.loc;
nodet = K == 0 || (K == 1 && isequal(D.br,[0 1]) && isequal(D.C{1},1) && ~D.loc);
if ~(trendcase || nodet), error('dfdet:args','dfdet_power_envelope: the design must be dfdet_trend(-1), dfdet_trend(0) or dfdet_trend(1) (the limit of the point-optimal statistic of other designs: pt_envelope).'); end
P = zeros(size(c));
for k = 1:numel(c)
    ck = double(c(k));
    if ~(ck < 0), if isnan(ck), P(k) = NaN; else, P(k) = alpha; end, continue, end   % (NaN at NaN, as POWER_ENVELOPE)
    if ck == -Inf, P(k) = 1; continue, end                                            % (the limit)
    if nodet
        a = ck/2; D0 = dfdet_trend(-1);
        b = dfdet_lin_icdf(alpha,a,D0,0,struct('b0',0));
        P(k) = dfdet_lin_cdf(b,a,D0,ck);
    else
        a = ck*ck/(2*(ck - 1));
        b = dfdet_lin_icdf(alpha,a,D,0,struct('b0',0,'gls',ck));
        P(k) = dfdet_lin_cdf(b,a,D,ck,struct('gls',ck));
    end
end
end
