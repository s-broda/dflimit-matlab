function f = dfdet_pdf(x,D,c,opts)
%DFDET_PDF  Limiting density of the t-ratio for any deterministic terms (the general engine).
%
%   f = DFDET_PDF(x,D,c,opts): the density of the limit of the t-ratio;
%   arguments and options as in DFDET_CDF.  For a design whose numerator has
%   rank one (no deterministic terms, GLS detrending) with c0 > 0 the density
%   is infinite at the tangency x = -sqrt(4 c0 a) of the parabola
%   X(v) = c0 v^2 + x v + a (a = 1/2, or (1 + g0^2)/2 with a fixed start), and
%   Inf is returned there; elsewhere it is continuous.  Far in the tails the
%   computed density carries an absolute error of order 1e-14 and can be a
%   tiny negative number (it is returned as it is).
%
%   Example: dfdet_pdf(-3.5, dfdet_perron('A', [0.3 0.6]))
%
%   See also DFDET_CDF, DFDET_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
f = dfdet_dist(x,D,c,opts,'pdf');
end
