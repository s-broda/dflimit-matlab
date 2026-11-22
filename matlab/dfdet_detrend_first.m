function D = dfdet_detrend_first(D)
%DFDET_DETREND_FIRST  The design for the statistic that detrends first (Perron and Vogelsang).
%
%   D = DFDET_DETREND_FIRST(D): y_t is regressed on d(t/n), t = 1..n, and the
%   t-ratio is that of ytilde_{t-1} in the regression of Delta ytilde_t on
%   ytilde_{t-1}, t = 2..n, ytilde the residual, with a dummy for the single
%   observation after each jump of d.  The one-step regression of y_t on
%   d(t/n) and y_{t-1} is invariant to the coefficients of the trend only when
%   the span of d is closed under differentiation (polynomial trends, Perron's
%   models A and C); for a break in the slope alone (model B, the trend joined
%   at the break) it is not, and the test that detrends first is the valid
%   one.  Its limit is U_d/sqrt(V) with the same V and U_d = U - m' G^{-1}
%   int d' Jtilde (d' the piecewise derivative, Jtilde the projection
%   residual): U_d + 1/2 is a quadratic form in the same functionals with the
%   block of z sym((Pm Phi_n)(Pm Phi_d)) in place of sym(Pm Phi_n).  Zero
%   start only; the distribution function, the density, the quantiles, a local
%   null value and the coefficient statistic are supported; the moments, a
%   fixed and a stationary start and GLS detrending (which detrends first
%   already) are refused.
%
%   Example: the 5% critical value with a break in the slope at 0.5
%      dfdet_icdf(0.05, dfdet_detrend_first(dfdet_perron('B', 0.5)))
%
%   See also DFDET_PERRON, DFDET_CDF.
if ~isstruct(D) || ~isfield(D,'br'), error('dfdet:args','dfdet_detrend_first: D must be a design (DFDET_DETERM).'); end
D.name = [D.name ' (detrended first)'];
D.detrend = 'first';
end
