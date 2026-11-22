function F = dfdet_dist(x,D,c,opts,mode)
%DFDET_DIST  The engine dispatch of DFDET_CDF and DFDET_PDF (mode 'cdf' or 'pdf').
%
%   Validates the arguments (DFDET_ARGS), builds the model of the design with
%   its start and its outer contour (cached), sets the non-finite points
%   (F = 0 at -Inf, 1 (cdf) or 0 (pdf) at +Inf, NaN at NaN) and evaluates the
%   finite ones by one evaluation of the contour integral (DFDET_CONTOUR) on
%   the path, refined on the vertical line to the oscillation at the largest
%   |x| (opts.verified = false, the default) or by the checked computation
%   (DFDET_VERIFY, opts.verified = true, tolerance opts.vtol = 1e-9).  For a
%   design whose numerator has rank one (no deterministic terms, GLS
%   detrending) with c0 > 0 the density is infinite at the tangency
%   x = -sqrt(4 c0 a) of the parabola X(v) (a the shift of Mt), and Inf is
%   returned there.
%
%   See also DFDET_CDF, DFDET_PDF, DFDET_CONTOUR.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[c,c0,start,g0,gls,o] = dfdet_args(D,c,opts);
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','dfdet_%s: x must be real.',mode); end
x = double(x);
F = zeros(size(x)); fin = isfinite(x);
F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
if ~any(fin(:)), return, end
Mo = dfdet_startmodel(D,c,start,g0,gls);
CP = dfdet_cpath(Mo,c,struct('phi',getf(o,'phi',[]),'hfac',getf(o,'hfac',1)));
xf = x(fin); xf = xf(:).';
if strcmp(mode,'cdf'), which = 1; else, which = 2; end
co = struct('c0',c0);
for f = {'n','ntail','alpha','filon','warn','tol'}, if isfield(o,f{1}), co.(f{1}) = o.(f{1}); end, end
if logical(getf(o,'verified',false))
    v = dfdet_verify(@(CPv,oo) ev(CPv,xf,merge(co,oo),which),CP,getf(o,'vtol',1e-9));
else
    v = ev(CP,xf,co,which);
end
if which == 2 && c0 > 0 && Mo.r == 1 && Mo.rm == 0
    xs_ = -sqrt(4*c0*Mo.a);
    v(abs(xf - xs_) <= 1e-12*max(1,abs(xf))) = Inf;
end
F(fin) = v;
end

function v = ev(CP,xf,co,which)
CPr = dfdet_contour('refine',CP,xf);
[F,f] = dfdet_contour('tratio',CPr,xf,co);
if which == 1, v = F; else, v = f; end
end

function o = merge(o,p)
for f = fieldnames(p)', o.(f{1}) = p.(f{1}); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
