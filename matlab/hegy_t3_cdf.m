function F = hegy_t3_cdf(x,D,c,opts)
%HEGY_T3_CDF  Limiting distribution function of the first seasonal t-ratio of Hylleberg, Engle, Granger and Yoo.
%
%   F = HEGY_T3_CDF(x,D,c,opts) is the limit of P(t3 <= x) at the points x (an
%   array) for the first t-ratio of the regression of Hylleberg, Engle,
%   Granger and Yoo at a harmonic frequency (strictly between 0 and pi: pi/2
%   for quarterly data, where t3 is the t-ratio of pi_3, the coefficient of
%   y_{3,t-2}; each pair of harmonic frequencies of monthly data), the one of
%   Dickey-Fuller shape: t3 = (R - c0 V)/sqrt(V), R = Re int Z~* dZ, V = int
%   |Z~|^2 of the complex Ornstein-Uhlenbeck process dZ = c Z dr + dW (W =
%   W_a + i W_b), Z(0) = 0, and Z~ its residual on the deterministic terms of
%   the design D (the real trend functions with complex coefficients):
%   dfdet_trend(-1) for a regression without seasonal terms (a constant or a
%   linear trend leave the harmonic statistics alone), dfdet_trend(0) for
%   seasonal dummies (with or without a trend), dfdet_trend(1) for seasonal
%   dummies and seasonal trends, and any design of the general engine (breaks
%   in the seasonal means, ...).  c (<= 0, default 0) and opts.c0 (0) in the
%   units of the roots' modulus 1 + c/n, n the number of years; opts.vol a
%   volatility path (VOLPATH; c may then be a vector of one value per regime).
%
%   The moment generating function of (R, V) is e^{-s}/D(s,t) with the general
%   engine's determinant D, whose singularities in s are poles: the inner
%   integral of the inversion on paths is the sum of the residues at the roots,
%   and the distribution function is one integral over the engine's outer
%   contour of closed-form functions.  Where two roots of a group come within
%   1e-4 (relative) of each other along the contour, the engine's inversion on
%   paths of the doubled transform (every root twice) computes instead.
%   opts.verified = true: both, the latter returned with a warning where they
%   differ by more than opts.vtol (1e-9).  Zero start, least squares
%   detrending.  0, 1 and NaN at -Inf, Inf and NaN.
%
%   Example: hegy_t3_cdf(-3.3442, dfdet_trend(0))        % 0.05: the 5% critical value with seasonal dummies
%
%   See also HEGY_T3_PDF, HEGY_T3_QUANTILE, VOLPATH, DFDET_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
F = t3_values('hegy_t3_cdf',x,D,c,opts,1);
cm = min(c(:)); c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
F = ltratio_rangecheck('hegy_t3_cdf',sprintf('the first seasonal t-ratio, c = %g, c0 = %g',cm,c0),'x',x,F);
end

function v = t3_values(fn,x,D,c,opts,which)
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','%s: c0 must be a finite real scalar',fn); end
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','%s: x must be real.',fn); end
CP = dfdet_t3('path',fn,D,c,opts);
x = double(x); v = zeros(size(x)); fin = isfinite(x);
v(isnan(x)) = NaN; if which == 1, v(x == Inf) = 1; end
if ~any(fin(:)), return, end
xf = x(fin); xf = xf(:).';
[F,f] = dfdet_t3('eval',CP,xf,c0);
if isfield(opts,'verified') && ~isempty(opts.verified) && opts.verified
    vtol = 1e-9; if isfield(opts,'vtol') && ~isempty(opts.vtol), vtol = opts.vtol; end
    [F2,f2] = dfdet_t3('eval',CP,xf,c0,40,64);
    if which == 1, d = max([0 abs(F - F2)]); else, d = max([0 abs(f - f2)]); end
    if ~(d <= vtol)
        warning('dfdet:verified','%s, verified: the residues and the inversion on paths of the doubled transform differ by %.2g (above %g); the latter is returned.',fn,d,vtol);
        F = F2; f = f2;
    end
end
if which == 1, v(fin) = F; else, v(fin) = f; end
end
