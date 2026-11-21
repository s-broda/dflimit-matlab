function P = ltratio_lincdf(b,a,cas,c,opts)
%LTRATIO_LINCDF  Distribution function of the linear form U - aV of the limiting functionals.
%
%   P = LTRATIO_LINCDF(b,a,cas,c,opts) is the limit of P(U - aV <= b) under phi = 1 + c/n (b an array), zero start,
%   no drift, with U = int J dJ and V = int J^2 as in LTRATIO_CDF: case 'DK' (no deterministic terms, and the
%   constant case with GLS demeaning, whose limit is the same) or 'DG' (the GLS-detrended trend case of Elliott,
%   Rothenberg and Stock, detrended at opts.cbar, default -13.5).  The point-optimal test of Elliott, Rothenberg
%   and Stock against cbar < 0 rejects for U - aV < b: in case 'DK' with a = cbar/2 (its log likelihood ratio is
%   cbar U - cbar^2 V/2), in case 'DG', with the data detrended at cbar, with a = cbar^2/(2 (cbar - 1)) (it rejects
%   for small cbar^2 V + (1 - cbar) J(1)^2, J(1)^2 = 2U + 1); see POWER_ENVELOPE.  a = NaN gives P(V <= b) (case
%   DK).  By the inversion formula of Gil-Pelaez (TR_GILPELAEZ) where its transform has decayed by the end of its
%   grid and c <= 0, to about 1e-10, about 0.02 s per call in case DK, 0.1 s in case DG; otherwise (|a| below about
%   0.3 to 3, depending on c, which includes a = 0, the law of U alone, and c > 0) by a conditional integral: in case
%   DK that of the region {U <= aV + b} (LTRATIO_REGION), in case DG the conditional integral given U (TR_LINCOND),
%   to about 1e-10 and in [0, 1].  c must be at least -550; in case DG at most about 3.5, and for a = NaN about 7
%   (an error beyond).
%
%   Example: ltratio_lincdf(-0.2639,-3.5,'DK',0)   % 0.05: the size of the point-optimal test at cbar = -7
%
%   See also POWER_ENVELOPE, LTRATIO_REGION, TR_GILPELAEZ.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
if ~(ischar(cas) || isstring(cas)) || ~any(strcmpi(cas,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_lincdf: c must be a finite real scalar.'); end
if c < -550, error('ltratio:args','ltratio_lincdf: c must be at least -550 (the range in which the programs are checked; the moments fail from about c = -600, the line integral from about c = -650).'); end
if ~(isnumeric(a) && isscalar(a) && isreal(a)), error('ltratio:args','ltratio_lincdf: a must be a real scalar (NaN: the distribution of V).'); end
if ~(isnumeric(b) && isreal(b)), error('ltratio:args','ltratio_lincdf: b must be real.'); end
P = tr_gilpelaez(b,a,upper(char(cas)),c,opts);
end
