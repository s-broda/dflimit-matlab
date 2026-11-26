function [F,info] = adfmax_cdf(x,cas,c,opts)
%ADFMAX_CDF  Limiting distribution function of the ADFmax test of Leybourne (1995).
%
%   F = ADFMAX_CDF(x,cas,c,opts) is the limit as n -> Inf of P[ADFmax <= x] under phi = 1 + c/n, no drift, for the
%   larger of the Dickey-Fuller t-ratios of the forward regression (y_t on y_{t-1}) and of the reverse regression
%   (y_t on y_{t+1}), with a constant ('DU', 'SU') or a constant and a trend ('DT', 'ST'), zero start ('D') or
%   stationary start ('S', c <= -1e-4).  In the limit both regressions have the denominator V and the numerators U and
%   -U - 1, so that ADFmax => (|R| - 1/2)/sqrt(V), R = U + 1/2, and {ADFmax <= x} = {|R| <= 1/2 + x sqrt(V)}, a
%   region between two curves of the plane (the ADFmax test of the companion paper), by the conditional integral (TR_ADFMAX).  x may
%   be an array; F = 0 at -Inf, 1 at Inf, NaN at NaN.  Accurate to about 1e-11; 0.1 to 0.5 s per value with a zero
%   start and 1 to 3 s with a stationary start on the MATLAB engines (opts are passed to TR_LIMDIST).  info: the
%   error estimates and statuses of the engine; a nonzero status gives the warning tr_limdist:status.  A port of
%   adfmax_cdf of julia/DFExact.jl.
%
%   Example: adfmax_cdf(-2.43542,'DU')    % 0.05 (the quantiles of the ADFmax statistic in the companion paper)
%
%   See also ADFMAX_PDF, ADFMAX_ICDF, LTRATIO_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[F,info] = adfmax_dist(x,cas,c,opts,'cdf');
end
