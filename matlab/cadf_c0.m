function c0 = cadf_c0(c,rho,R)
%CADF_C0  The local null value of the covariate-augmented Dickey-Fuller t-ratio of Hansen (1995).
%
%   c0 = CADF_C0(c,rho,R) = c (1 - 1/(rho R)), 0 at c = 0 whatever rho and R:
%   rho the long-run correlation of the innovation of the series and the
%   error of the augmented regression and R the ratio of their long-run
%   standard deviations (R = |rho|, the default, for the well-specified
%   regression; Hansen's rho^2 and R^2 are the squares); rho in [-1, 1] and
%   not 0.  c is the local-to-unity parameter of this package (phi = 1 + c/n;
%   Hansen's c with the opposite sign).  With opts.rho = rho and opts.c0 =
%   CADF_C0(c,rho,R), DFDET_CDF and LTRATIO_CDF give the limit of the
%   statistic, rho tau_c0 + sqrt(1 - rho^2) Z.
%
%   Example: cadf_c0(-5, sqrt(0.5))      % 5: -c (1 - rho^2)/rho^2 when R = rho
%
%   See also PREDICTIVE_C0, DFDET_CDF, LTRATIO_CDF.
if nargin < 3 || isempty(R), R = abs(rho); end
if ~(isnumeric(rho) && isscalar(rho) && isreal(rho) && isfinite(rho) && rho >= -1 && rho <= 1 && rho ~= 0), error('dfdet:args','cadf_c0: rho must lie in [-1, 1] and not be 0'); end
if ~(isnumeric(R) && isscalar(R) && isreal(R) && isfinite(R) && R > 0), error('dfdet:args','cadf_c0: R must be finite and positive'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','cadf_c0: c must be finite'); end
if c == 0, c0 = 0; else, c0 = double(c)*(1 - 1/(double(rho)*double(R))); end
end
