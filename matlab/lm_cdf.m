function F = lm_cdf(x,D,c)
%LM_CDF  Limiting distribution function of the LM unit root statistic (Schmidt and Phillips; with breaks, Lee and Strazicich).
%
%   F = LM_CDF(x,D,c) is the limit of P(tau_LM <= x) (x an array) for the LM
%   statistic with the trend functions of D (the constant and a linear trend
%   in their span), under phi = 1 + c/n (c defaults to 0): the trend
%   coefficients estimated from the first differences, S_t = y_t - psi -
%   Z_t delta, and the t-ratio of S_{t-1}
%   in the regression of the differences on the differenced trend functions
%   and S_{t-1}.  Its limit is -1/(2 sqrt V) (DFDET_PREP, kind 'lm'), so that
%   F(x) = P(V <= 1/(4 x^2)) for x < 0 and 1 for x >= 0, one inversion of the
%   transform of V.  The limit holds only for trend functions whose span
%   contains the constant (the LM regression estimates the level from the
%   first observation): for any other design, e^r alone or r alone, say,
%   LM_CDF, LM_QUANTILE and LM_POWER raise an error (decided on the values of
%   the piecewise polynomial design at Chebyshev points of each piece: exactly
%   for polynomial trends, breaks and dummies, to 1e-6 for the interpolants of
%   DFDET_FDESIGN and DFDET_DESIGN_FROM_MATRIX).
%
%   Example: lm_cdf(-3.0230, dfdet_trend(1))   % 0.05
%
%   See also LM_QUANTILE, LM_POWER.
if nargin < 3 || isempty(c), c = 0; end
CP = dfdet_lmpath(D,c);
F = zeros(size(x));
for k = 1:numel(x)
    if isnan(x(k)), F(k) = NaN; elseif x(k) >= 0, F(k) = 1; elseif x(k) == -Inf, F(k) = 0;
    else, F(k) = dfdet_lmpath(CP,1/(4*x(k)^2)); end
end
end
