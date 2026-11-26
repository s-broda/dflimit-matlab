function f = hegy_t3_pdf(x,D,c,opts)
%HEGY_T3_PDF  Limiting density of the first seasonal t-ratio of Hylleberg, Engle, Granger and Yoo.
%
%   f = HEGY_T3_PDF(x,D,c,opts) is the density of the limit of t3 (HEGY_T3_CDF,
%   the same arguments and options).
%
%   Example: hegy_t3_pdf([-3 -1 0], dfdet_trend(0))
%
%   See also HEGY_T3_CDF, HEGY_T3_QUANTILE.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
f = t3_values('hegy_t3_pdf',x,D,c,opts,2);
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
