function [F,info] = adfmax_dist(x,cas,c,opts,mode)
%ADFMAX_DIST  Argument checks and non-finite x of ADFMAX_CDF and ADFMAX_PDF (mode 'cdf' or 'pdf').
%
%   See also ADFMAX_CDF, TR_ADFMAX.
if ~(ischar(cas) || isstring(cas)) || ~any(strcmpi(cas,{'DU','DT','SU','ST'}))
    error('ltratio:args','adfmax: the case must be one of DU, DT, SU, ST.');
end
cas = upper(char(cas));
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','adfmax: c must be a finite real scalar.'); end
if cas(1) == 'S' && ~(c <= -1e-4), error('ltratio:args','ltratio: a stationary start (cases SK, SU, ST) requires c <= -1e-4: closer to zero the variance -1/(2c) of the start is too large for the quadrature.'); end   % (LTRATIO_ARGS)
if c < -550, error('ltratio:args','adfmax: c must be at least -550 (the range in which the programs are checked).'); end
if ~(isnumeric(x) && isreal(x)), error('ltratio:args','adfmax: x must be real.'); end
F = zeros(size(x)); F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
fin = isfinite(x); info = struct('err',NaN(size(x)),'status',zeros(size(x)));
if any(fin(:))
    [Fv,inf1] = tr_adfmax(x(fin),cas,c,opts,mode);
    F(fin) = Fv; info.err(fin) = inf1.err; info.status(fin) = inf1.status;
end
if strcmp(mode,'pdf'), F(F < 0) = 0; end
end
