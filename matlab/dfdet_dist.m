function [F,f] = dfdet_dist(x,D,c,opts,mode)
%DFDET_DIST  The engine dispatch of DFDET_CDF and DFDET_PDF (mode 'cdf' or 'pdf'; 'both': [F,f] from one pass).
%
%   Validates the arguments (DFDET_ARGS), builds the model of the design with
%   its start and its outer contour (cached), sets the non-finite points
%   (F = 0 at -Inf, 1 (cdf) or 0 (pdf) at +Inf, NaN at NaN) and evaluates the
%   finite ones by one evaluation of the contour integral (DFDET_CONTOUR) on
%   the path, refined on the vertical line to the oscillation at the largest
%   |x| (opts.verified = false, the default; by the rule of one evaluation of
%   DFDET_RULES, opts.rule, 'auto' by default, with opts.hfac, opts.n,
%   opts.ntail in its place where given) or by the checked computation
%   (DFDET_VERIFY, opts.verified = true, tolerance opts.vtol = 1e-9, on the
%   present rule's contour).  Mode 'both' returns the distribution function
%   and the density from one pass (dist_both of DFDet.jl: the values of the
%   modes 'cdf' and 'pdf' bit for bit; the quantile search of LTRATIO_ICDF).  For a
%   design whose numerator has rank one (no deterministic terms, GLS
%   detrending) with c0 > 0 the density is infinite at the tangency
%   x = -sqrt(4 c0 a) of the parabola X(v) (a the shift of Mt), and Inf is
%   returned there.
%
%   See also DFDET_CDF, DFDET_PDF, DFDET_CONTOUR.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if isstruct(opts) && dfdet_ext('has',c,opts), [F,f] = dfdet_ext('dist',x,D,c,opts,mode); return, end   % (vol, rho: DFDET_EXT)
[c,c0,start,g0,gls,o] = dfdet_args(D,c,opts,sprintf('the t-ratio (dfdet_%s)',mode));
dfdet_method(opts,['dfdet_' mode]);
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','dfdet_%s: x must be real.',mode); end
x = double(x);
F = zeros(size(x)); fin = isfinite(x);
F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
f = zeros(size(x)); f(isnan(x)) = NaN;
if ~any(fin(:)), return, end
Mo = dfdet_startmodel(D,c,start,g0,gls);
verified = logical(getf(o,'verified',false));
% (the rule of one evaluation and its outer contour, DFDET_RULES; verified: the present rule)
[hf,n,nt,rl] = dfdet_rules('rule',getf(o,'rule','auto'),verified,Mo,c,c0,o);
CP = dfdet_rules('path',Mo,c,getf(o,'phi',[]),hf,rl,c0);
if ~verified, dfdet_farwarn(sprintf('dfdet_%s',mode),D,c,start,g0,gls); end
xf = x(fin); xf = xf(:).';
switch mode, case 'cdf', which = 1; case 'pdf', which = 2; otherwise, which = 3; end
co = struct('c0',c0,'n',n,'ntail',nt);
for fl = {'alpha','filon','warn','tol'}, if isfield(o,fl{1}), co.(fl{1}) = o.(fl{1}); end, end
if verified
    v = dfdet_verify(@(CPv,oo) ev(CPv,xf,merge(co,oo),which),CP,getf(o,'vtol',1e-9));
else
    v = ev(CP,xf,co,which);
end
if which == 3, w = v(2,:); v = v(1,:); end
if which >= 2 && c0 > 0 && Mo.r == 1 && Mo.rm == 0
    xs_ = -sqrt(4*c0*Mo.a); tg = abs(xf - xs_) <= 1e-12*max(1,abs(xf));
    if which == 2, v(tg) = Inf; else, w(tg) = Inf; end
end
F(fin) = v;
if which == 3, f(fin) = w; end
end

function v = ev(CP,xf,co,which)
CPr = dfdet_contour('refine',CP,xf,[],co.c0);
[F,f] = dfdet_contour('tratio',CPr,xf,co);
if which == 1, v = F; elseif which == 2, v = f; else, v = [F(:).'; f(:).']; end
end

function o = merge(o,p)
for f = fieldnames(p)', o.(f{1}) = p.(f{1}); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
