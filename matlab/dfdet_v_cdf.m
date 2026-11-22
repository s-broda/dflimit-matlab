function P = dfdet_v_cdf(y,D,c,opts)
%DFDET_V_CDF  Distribution function of the denominator V = int Jt^2 of the t-ratio for any deterministic terms.
%
%   P = DFDET_V_CDF(y,D,c,opts) is the limit of P(V <= y) (y an array) by the
%   single inversion of M_V(t) = M(0,t) on the outer contour:
%   P(V > y) = (1/2 pi i) int_C M_V(t) e^{-t y} dt/t, M_V(t) from the roots and
%   log P(t) along the contour (the first term of the inversion on paths, a
%   probability for V alone), each logarithm continued from node to node.
%   Where the t-ratio takes the vertical line (c > 0, GLS detrending with a
%   small |cbar|) the law of V takes the hyperbola, on which e^{-t y} decays.
%   Options start, g0, gls as in DFDET_CDF.
%
%   Example: dfdet_v_cdf(0.05, dfdet_trend(1))
%
%   See also DFDET_V_PDF, DFDET_LIN_CDF, MSB_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
P = dfdet_stat('v','cdf',y,D,c,opts);
end
