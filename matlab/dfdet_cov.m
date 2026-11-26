function varargout = dfdet_cov(op,varargin)
%DFDET_COV  A normal component in the general engine: T = rho tau_c0 + sqrt(1 - rho^2) Z (opts.rho).
%
%   [F,f] = DFDET_COV('dist',x,D,c,opts,mode) and q = DFDET_COV('icdf',p,D,c,
%   opts): the distribution function, density and quantiles of the limit
%   T = rho tau_c0 + sqrt(1 - rho^2) Z, tau_c0 = (U - c0 V)/sqrt(V) the
%   t-ratio's limit of the design D at the local null value c0 and Z standard
%   normal and independent (the covariate-augmented t-ratio of Hansen (1995)
%   at c0 = CADF_C0(c,rho,R), the t-ratio of the slope in a predictive
%   regression at c0 = PREDICTIVE_C0(c,rho,b)).  rho in [-1, 1]: rho = 0 the
%   standard normal law; negative rho by the reflection P(T <= x) = 1 -
%   P(|rho| tau + rb Z <= -x).  opts.method: 'auto' (the default), 'lines'
%   (the inversion on two lines with the normal component in the kernel, on
%   the engine's outer contour; a zero or stationary start, c <= 0) or
%   'mixture' (the normal average of the t-ratio's distribution function,
%   graded toward its kink when c0 > 0); 'auto' takes the lines for
%   0.1 <= rho^2 <= 0.99 with a zero or stationary start at c <= 0 and the
%   mixture otherwise.  opts.verified: both routes, where they differ by more
%   than opts.vtol (1e-9) the mixture with half its step on the checked
%   distribution function, with a warning.  Least squares detrending and the
%   one-step regression only (DFDet_cov.jl).
%
%   See also CADF_C0, PREDICTIVE_C0, DFDET_CDF.
switch op
    case 'dist', [F,f] = cov_eval(varargin{:}); varargout = {F,f};
    case 'icdf', varargout = {cov_icdf(varargin{:})};
    otherwise, error('dfdet:args','dfdet_cov: unknown operation %s.',op);
end
end

function [rho,c0,start,g0,meth] = cov_args(fn,D,opts,extra)
% the checks of the arguments with rho (cov_args and the options that do not apply)
rho = getf(opts,'rho',1);
if ~(isnumeric(rho) && isscalar(rho) && isreal(rho) && isfinite(rho) && rho >= -1 && rho <= 1), error('dfdet:args','%s: rho must lie in [-1, 1]',fn); end
rho = double(rho);
meth = lower(char(getf(opts,'method','auto')));
if ~any(strcmp(meth,{'auto','lines','mixture'}))
    error('dfdet:args','%s: with rho ~= 1 the method must be ''auto'', ''lines'' (the inversion on two lines with the normal component in the kernel) or ''mixture'' (the normal average of the t-ratio''s distribution function)',fn);
end
if ~isempty(getf(opts,'gls',[])), error('dfdet:args','%s: rho ~= 1 (a normal component: the covariate-augmented t-ratio, the predictive regression) is implemented for least squares detrending; GLS detrending (gls) is not supported with it',fn); end
if ~strcmp(D.detrend,'onestep'), error('dfdet:args','%s: rho ~= 1 is implemented for the one-step regression; the statistic that detrends first (detrend_first) is not supported with it',fn); end
other = {};
for f = [{'n','ntail','alpha','filon','warn'} extra], if isfield(opts,f{1}) && ~isempty(opts.(f{1})), other{end+1} = f{1}; end, end %#ok<AGROW>
if ~isempty(other), error('dfdet:args','%s: with rho ~= 1 the options %s do not apply',fn,strjoin(other,', ')); end
c0 = getf(opts,'c0',0);
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','%s: c0 must be a finite real scalar',fn); end
g0 = getf(opts,'g0',0);
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('dfdet:args','dfdet: opts.g0 must be a finite real scalar.'); end
start = getf(opts,'start','');
if isempty(start), if g0 ~= 0, start = 'fixed'; else, start = 'zero'; end, end
start = lower(char(start));
if ~any(strcmp(start,{'zero','fixed','stationary'})), error('dfdet:args','dfdet: opts.start must be ''zero'', ''fixed'' or ''stationary''.'); end
if strcmp(start,'fixed') && g0 == 0, start = 'zero'; end
end

function m = cov_method(fn,meth,rho,c,fixed)
% the lines for 0.1 <= rho^2 <= 0.99 with a zero or stationary start at c <= 0, the mixture otherwise (near |rho| = 1 the
% kernel of the lines oscillates, its phase reaching x^2/(2 (1 - rho^2)), and the mixture is the cheaper; the lines do not
% cover a fixed start or the vertical line above the unit root)
if strcmp(meth,'lines')
    if fixed, error('dfdet:args','%s: method = ''lines'' (the inversion on two lines) needs a zero or a stationary start; method = ''mixture'' covers a fixed start',fn); end
    if c > 0, error('dfdet:args','%s: method = ''lines'' (the inversion on two lines) is implemented below and at the unit root (c <= 0); method = ''mixture'' covers c > 0',fn); end
    m = 'lines'; return
end
if strcmp(meth,'mixture'), m = 'mixture'; return, end
if fixed || c > 0 || ~(0.1 <= rho^2 && rho^2 <= 0.99), m = 'mixture'; else, m = 'lines'; end
end

function [F,f] = cov_eval(x,D,c,opts,mode)
fn = ['dfdet_' mode]; if strcmp(mode,'both'), fn = 'dfdet_cdf'; end
[rho,c0,start,g0,meth] = cov_args(fn,D,opts,{'tol'});
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','%s: x must be real.',fn); end
x = double(x);
if rho == 0, F = (erfc(-x/sqrt(2)))/2; f = exp(-x.*x/2)/sqrt(2*pi); return, end
if rho < 0
    % (the finite points by the reflection, the others set here: 0 at -Inf, 1 (cdf) or 0 (pdf) at Inf, NaN at NaN)
    F = zeros(size(x)); f = F; fin = isfinite(x);
    F(isnan(x)) = NaN; f(isnan(x)) = NaN; F(x == Inf) = 1;
    if ~any(fin(:)), return, end
    xf = x(fin);
    if rho == -1
        o1 = struct('c0',c0,'start',start,'g0',g0,'verified',getf(opts,'verified',false));
        [Fv,fv] = dfdet_dist(-xf,D,c,o1,'both');
    else
        o1 = opts; o1.rho = -rho; [Fv,fv] = cov_eval(-xf,D,c,o1,'both');
    end
    F(fin) = 1 - Fv; f(fin) = fv;
    return
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','dfdet: c must be a finite real scalar.'); end
fixed = strcmp(start,'fixed') && g0 ~= 0;
m = cov_method(fn,meth,rho,c,fixed);
F = zeros(size(x)); f = F; fin = isfinite(x);
F(isnan(x)) = NaN; f(isnan(x)) = NaN; F(x == Inf) = 1;
if ~any(fin(:)), return, end
xf = x(fin); xf = xf(:).';
if logical(getf(opts,'verified',false))
    [Fv,fv] = cov_verified(D,xf,c,rho,c0,start,g0,m,getf(opts,'vtol',1e-9));
else
    [Fv,fv] = cov_both(D,xf,c,rho,c0,start,g0,m,false,[]);
end
F(fin) = Fv; f(fin) = fv;
end

function [F,f] = cov_both(D,x,c,rho,c0,start,g0,m,verified,h)
if strcmp(m,'lines')
    Mo = dfdet_startmodel(D,c,start,g0,[]);
    CP = dfdet_cpath(Mo,c);
    [F,f] = dfdet_contour('lines',CP,x,struct('rho',rho,'c0',c0));
    return
end
[F,f] = mixture_both(D,x,c,rho,c0,start,g0,h,verified);
end

function [F,f] = mixture_both(D,x,c,rho,c0,start,g0,h,verified)
% P(rho tau_c0 + rb Z <= x) and its density from the t-ratio's distribution function (mixture_both of DFDet_cov.jl):
% P(T <= x) = int F((x - rb z)/rho) phi(z) dz, f_T(x) = -(1/rb) int F((x - rb z)/rho) z phi(z) dz, the trapezoidal rule
% with step h in z on |z| <= 8.5 (0.25/(1 + rb/rho), half that for a numerator of rank one); for c0 > 0 Gauss-Legendre
% panels of 12 nodes graded toward the kink of the distribution function at q_* = -2 sqrt(c0 a), no longer than 4h; F is
% taken as 0 or 1 outside a generous bound on the bulk of tau (tau_bulk)
zmax = 8.5; n = 12;
rb = sqrt(1 - rho^2);
Mo = dfdet_startmodel(D,c,start,g0,[]);
if isempty(h), if Mo.r == 1, h = 0.125/(1 + rb/rho); else, h = 0.25/(1 + rb/rho); end, end
if c0 > 0, qstar = -2*sqrt(c0*Mo.a); else, qstar = NaN; end
[gz,gw] = glrule(n); gz = gz(:); gw = gw(:);
nx = numel(x); Z = cell(1,nx); Wz = Z;
for i = 1:nx
    zstar = (x(i) - rho*qstar)/rb;
    if c0 > 0 && abs(zstar) < zmax
        pan = graded_panels(-zmax,zmax,zstar,1e-8);
        zz = cell(size(pan,1),1); ww = zz;
        for k = 1:size(pan,1)
            a0 = pan(k,1); b0 = pan(k,2); m = max(1,ceil((b0 - a0)/(4*h)));
            ab = a0 + (b0 - a0)*(0:m)/m; L = diff(ab);
            zz{k} = reshape(ab(1:end-1) + L.*gz,[],1); ww{k} = reshape(L.*gw,[],1);
        end
        Z{i} = vertcat(zz{:}); Wz{i} = vertcat(ww{:});
    else
        nz = ceil(zmax/h); Z{i} = (-nz:nz)'*h; Wz{i} = h*ones(2*nz + 1,1);
    end
end
len = cellfun(@numel,Z); q = cell(1,nx);
for i = 1:nx, q{i} = (x(i) - rb*Z{i})/rho; end
qa = vertcat(q{:});
[lo,hi] = tau_bulk(c,c0,start,g0);
Fq = double(qa >= lo);
ib = lo <= qa & qa <= hi;
if any(ib)
    Fq(ib) = dfdet_cdf(qa(ib),D,c,struct('c0',c0,'start',start,'g0',g0,'verified',verified));
end
F = zeros(1,nx); f = F; k0 = 0;
for i = 1:nx
    sel = k0 + (1:len(i)); k0 = k0 + len(i);
    w = Wz{i}.*exp(-Z{i}.^2/2)/sqrt(2*pi).*min(max(Fq(sel),0),1);
    F(i) = sum(w); f(i) = -sum(w.*Z{i})/rb;
end
end

function [lo,hi] = tau_bulk(c,c0,start,g0)
% m -+ 40 s with m = (c - c0) sqrt(E int J^2) and s = 1 + |c - c0| sqrt(E int J^2), E int J^2 of the process with its start
if abs(c) < 1e-6, ej = 0.5; else, ej = (expm1(2*c) - 2*c)/(4*c^2); end
if strcmp(start,'stationary'), ej = -1/(2*c); end
if strcmp(start,'fixed'), if abs(c) < 1e-6, ej = ej + g0^2; else, ej = ej + g0^2*expm1(2*c)/(2*c); end, end
r = sqrt(ej); m = (c - c0)*r; s = 1 + abs(c - c0)*r;
lo = m - 40*s; hi = m + 40*s;
end

function [F,f] = cov_verified(D,x,c,rho,c0,start,g0,m,vtol)
% the other route as well; where the two differ by more than vtol the mixture with half its step on the checked
% distribution function of tau, with a warning (cov_verified)
[F,f] = cov_both(D,x,c,rho,c0,start,g0,m,false,[]);
rb = sqrt(1 - rho^2); Mo = dfdet_startmodel(D,c,start,g0,[]);
if Mo.r == 1, hh = 0.0625/(1 + rb/rho); else, hh = 0.125/(1 + rb/rho); end
if strcmp(m,'lines'), [F2,f2] = cov_both(D,x,c,rho,c0,start,g0,'mixture',false,[]);
else, [F2,f2] = cov_both(D,x,c,rho,c0,start,g0,'mixture',true,hh); end
d = max([0 abs(F - F2) abs(f - f2)]);
if d <= vtol, return, end
if strcmp(m,'lines'), [F2,f2] = cov_both(D,x,c,rho,c0,start,g0,'mixture',true,hh); end
warning('dfdet:verified','dfdet, verified: the two computations of the normal component differ by %.2g (above %g); the normal average on the checked distribution function is returned (c = %g, c0 = %g, rho = %g).',d,vtol,c,c0,rho);
F = F2; f = f2;
end

function q = cov_icdf(p,D,c,opts)
% the quantiles with rho ~= 1: Newton's method with the density, the step limited to 1, then the safeguarded search, from
% opts.x0 or the normal approximation with the moments of the t-ratio (those of rho tau + rb Z); rho = 0 the normal
% quantile, rho < 0 by the reflection
fn = 'dfdet_icdf';
[rho,c0,start,g0,meth] = cov_args(fn,D,opts,{});
tol = getf(opts,'tol',1e-10); x0 = getf(opts,'x0',[]);
if rho == 0, q = sqrt(2)*erfinv(2*p - 1); return, end
if rho < 0
    o1 = opts; if ~isempty(x0), o1.x0 = -x0; end
    if rho == -1
        o1 = struct('c0',c0,'start',start,'g0',g0,'tol',tol); if ~isempty(x0), o1.x0 = -x0; end
        q = -dfdet_icdf(1 - p,D,c,o1);
    else
        o1.rho = -rho; q = -cov_icdf(1 - p,D,c,o1);
    end
    return
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','dfdet: c must be a finite real scalar.'); end
m = cov_method(fn,meth,rho,c,strcmp(start,'fixed') && g0 ~= 0);
mu = -2*rho; sd = 1;
if isempty(x0) && strcmp(D.detrend,'onestep')
    try
        mm = dfdet_mom(2,D,c,struct('c0',c0,'start',start,'g0',g0));
        mu = rho*mm(1); sd = sqrt(max(rho^2*(mm(2) - mm(1)^2) + 1 - rho^2,1e-6));
    catch ME
        if ~strncmp(ME.identifier,'dfdet:',6), rethrow(ME); end
    end
end
FF = @(y) ffun(D,y,c,rho,c0,start,g0,m);
q = zeros(size(p));
for k = 1:numel(p)
    if isempty(x0), xs = mu + sd*sqrt(2)*erfinv(2*p(k) - 1); else, xs = x0; end
    q(k) = dfdet_newton(FF,p(k),xs,1,tol,100,30);
end
end
function v = ffun(D,y,c,rho,c0,start,g0,m)
[F,f] = cov_both(D,y,c,rho,c0,start,g0,m,false,[]);
v = [F f];
end

function out = graded_panels(A,B,pos,h)
% the panels of [A, B] graded toward pos with the smallest panel h (graded_panels of the engine, one focus point)
pos = min(max(pos,A),B); tol = 1e-12*(B - A);
out = zeros(0,2);
if abs(pos - A) <= tol, out = seg(out,A,B,h); return, end
if abs(B - pos) <= tol, out = seg(out,B,A,h); return, end
out = seg(out,pos,A,h); out = seg(out,pos,B,h);
out = sortrows(out);
end
function out = seg(out,st,en,h0)
len = abs(en - st); if en >= st, dir = 1; else, dir = -1; end
if len == 0, return, end
a = 0; h = min(h0,len);
while a < len*(1 - 1e-14)
    b = min(a + h,len); if (len - b) < 0.5*h, b = len; end
    u1 = st + dir*a; u2 = st + dir*b;
    out(end+1,:) = [min(u1,u2) max(u1,u2)]; %#ok<AGROW>
    a = b; h = h*2;
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
