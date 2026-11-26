function V = volpath(breaks,omega)
%VOLPATH  Volatility (and persistence) regimes between known dates, for the general engine (opts.vol).
%
%   V = VOLPATH(breaks,omega): the regime boundaries breaks (strictly
%   increasing in (0, 1), fractions of the sample) and the standard
%   deviations omega of the innovations on the numel(breaks) + 1 regimes, of
%   any positive scale (they are scaled so that the integral of omega^2 is
%   one; only the ratios matter; default all equal).  VOLPATH([]) is one
%   regime, the engine itself.  The structure V has the fields lam
%   ([0 breaks 1]) and om (the scaled standard deviations).
%
%   It is the option opts.vol of DFDET_CDF, DFDET_PDF, DFDET_ICDF, LTRATIO_CDF,
%   LTRATIO_PDF, LTRATIO_ICDF and HEGY_T3_CDF, HEGY_T3_PDF, HEGY_T3_QUANTILE:
%   the limit (U - c0 V)/sqrt(V) of the t-ratio when the process is
%   dJ = c(r) J dr + omega(r) dW with omega and, through a vector c of one value
%   per regime, the local-to-unity parameter piecewise constant (zero start,
%   least squares detrending, the one-step regression, every c in [-480, 0]).
%
%   Examples
%      V = volpath(0.3, [1 3]);          % the standard deviation triples at 30% of the sample
%      V = volpath(0.5);                 % two regimes of equal variance (changing persistence: c a vector)
%      dfdet_cdf(-2.86, dfdet_trend(0), 0, struct('vol', volpath(0.1, [1 5])))
%
%   See also DFDET_CDF, LTRATIO_CDF, HEGY_T3_CDF.
if nargin < 1, breaks = []; end
br = double(breaks(:).');
if nargin < 2 || isempty(omega), omega = ones(1,numel(br) + 1); end
if ~(isnumeric(breaks) && isreal(breaks)) || ~all(isfinite(br) & br > 0 & br < 1), error('dfdet:args','volpath: the break dates must lie in (0, 1)'); end
if any(diff(br) <= 0), error('dfdet:args','volpath: the break dates must be strictly increasing'); end
if numel(omega) ~= numel(br) + 1
    error('dfdet:args','volpath: one standard deviation omega per regime: %d for %d break dates, %d given',numel(br) + 1,numel(br),numel(omega));
end
om = double(omega(:).');
if ~(isnumeric(omega) && isreal(omega)) || ~all(isfinite(om) & om > 0), error('dfdet:args','volpath: the standard deviations omega must be finite and positive'); end
lam = [0 br 1];
s = 0; for i = 1:numel(om), s = s + om(i)^2*(lam(i+1) - lam(i)); end
om = om/sqrt(s);
V = struct('lam',lam,'om',om,'breaks',br,'omega',double(omega(:).'));
end
