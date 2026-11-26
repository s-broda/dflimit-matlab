function varargout = dfdet_ext(op,varargin)
%DFDET_EXT  The extensions of the limit in the user-level routines of the general engine (opts.vol, opts.rho).
%
%   tf = DFDET_EXT('has',c,opts): whether the call has an extension (a
%   volatility path opts.vol, opts.rho ~= 1, or c a vector).
%   [F,f] = DFDET_EXT('dist',x,D,c,opts,mode): DFDET_CDF ('cdf'), DFDET_PDF
%   ('pdf') or both ('both') with the extensions (ext_eval of DFDet_cov.jl): a
%   volatility path (VOLPATH; c a number or a vector of one value per regime;
%   one regime with a scalar c is the engine itself) by the inversion on
%   paths of the volatility transform on its outer contour by the rule of one
%   evaluation of the inner model at the smallest c, or the checked
%   computation (opts.verified); or a normal component (DFDET_COV).
%   q = DFDET_EXT('icdf',p,D,c,opts): DFDET_ICDF with them, Newton's method
%   from opts.x0 or from the normal approximation with the moments of the
%   t-ratio at the mean of the c of the regimes, the step limited to 1, then
%   the safeguarded search.
%
%   See also VOLPATH, DFDET_VOL, DFDET_COV, DFDET_CDF.
switch op
    case 'has'
        c = varargin{1}; opts = varargin{2};
        varargout = {~isscalar(c) || (isfield(opts,'vol') && ~isempty(opts.vol)) || (isfield(opts,'rho') && ~isempty(opts.rho) && ~isequal(opts.rho,1))};
    case 'dist', [F,f] = ext_dist(varargin{:}); varargout = {F,f};
    case 'icdf', varargout = {ext_icdf(varargin{:})};
    otherwise, error('dfdet:args','dfdet_ext: unknown operation %s.',op);
end
end

function [F,f] = ext_dist(x,D,c,opts,mode)
fn = ['dfdet_' mode]; if strcmp(mode,'both'), fn = 'dfdet_cdf'; end
dcheck(D);
vol = getf(opts,'vol',[]);
if isempty(vol)
    if ~isscalar(c), error('dfdet:args','%s: c is a vector (one value per regime): that needs a volatility path, opts.vol = volpath(breaks, omega)',fn); end
    pathopts(fn,opts);
    [F,f] = dfdet_cov('dist',x,D,c,opts,mode);
    if strcmp(mode,'pdf'), F = f; end            % (the first output holds the values of the mode, as DFDET_DIST)
    return
end
if ~isequal(getf(opts,'rho',1),1), error('dfdet:args','%s: rho and vol together are not implemented (the limit with a normal component under a volatility path is not derived)',fn); end
[start,g0,gls] = startof(opts);
[VM,c1] = dfdet_vol('args',fn,D,c,opts,start,g0,gls);
if isempty(VM)
    o1 = rmfield(opts,'vol');
    [F,f] = dfdet_dist(x,D,c1,o1,mode);
    return
end
c0 = getf(opts,'c0',0);
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','%s: opts.c0 must be a finite real scalar.',fn); end
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','%s: x must be real.',fn); end
x = double(x);
F = zeros(size(x)); fin = isfinite(x);
F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
f = zeros(size(x)); f(isnan(x)) = NaN;
if strcmp(mode,'both'), F(x == Inf) = 1; end
if ~any(fin(:)), return, end
o = pathfields(opts); verified = logical(getf(o,'verified',false)); cm = min(VM.cr);
[hf,n,nt,rl] = dfdet_rules('rule',getf(o,'rule','auto'),verified,VM,cm,c0,o);
CP = dfdet_rules('path',VM,cm,getf(o,'phi',[]),hf,rl,c0);
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
Mo = CP.Mo;
if which >= 2 && c0 > 0 && Mo.r == 1 && Mo.rm == 0
    xs_ = -sqrt(4*c0*Mo.a); tg = abs(xf - xs_) <= 1e-12*max(1,abs(xf));
    if which == 2, v(tg) = Inf; else, w(tg) = Inf; end
end
F(fin) = v;
if which == 3, f(fin) = w; end
end

function q = ext_icdf(p,D,c,opts)
fn = 'dfdet_icdf';
dcheck(D);
if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('dfdet:args','dfdet_icdf: p must lie in (0,1).'); end
vol = getf(opts,'vol',[]); rho = getf(opts,'rho',1);
if isempty(vol)
    if ~isscalar(c), error('dfdet:args','%s: c is a vector (one value per regime): that needs a volatility path, opts.vol = volpath(breaks, omega)',fn); end
    pathopts(fn,opts);
    q = dfdet_cov('icdf',p,D,c,opts);
    return
end
if ~isequal(rho,1), error('dfdet:args','%s: rho and vol together are not implemented (the limit with a normal component under a volatility path is not derived)',fn); end
[start,g0,gls] = startof(opts);
[VM,c1] = dfdet_vol('args',fn,D,c,opts,start,g0,gls);
if isempty(VM), q = dfdet_icdf(p,D,c1,rmfield(opts,'vol')); return, end
c0 = getf(opts,'c0',0);
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','%s: opts.c0 must be a finite real scalar.',fn); end
tol = getf(opts,'tol',1e-10); x0 = getf(opts,'x0',[]);
mu = -2; sd = 1;
if isempty(x0)
    % (the normal approximation with the moments of the t-ratio at the mean of the c of the regimes, ext_icdf of merged.jl)
    try
        mo = struct('c0',c0); m = dfdet_mom(2,D,mean(VM.cr),mo);
        mu = m(1); sd = sqrt(max(m(2) - m(1)^2,1e-6));
    catch ME
        if ~strncmp(ME.identifier,'dfdet:',6), rethrow(ME); end
    end
end
o = pathfields(opts); cm = min(VM.cr);
[hf,n,nt,rl] = dfdet_rules('rule',getf(o,'rule','auto'),false,VM,cm,c0,o);
CP = dfdet_rules('path',VM,cm,getf(o,'phi',[]),hf,rl,c0);
co = struct('c0',c0,'n',n,'ntail',nt);
for fl = {'alpha','filon','warn'}, if isfield(o,fl{1}), co.(fl{1}) = o.(fl{1}); end, end
FF = @(y) ev(CP,y,co,3).';
q = zeros(size(p));
for k = 1:numel(p)
    if isempty(x0), xs = mu + sd*sqrt(2)*erfinv(2*p(k) - 1); else, xs = x0; end
    q(k) = dfdet_newton(FF,p(k),xs,1,tol,100,30);
end
end

function v = ev(CP,xf,co,which)
% (the paths of the volatility model are not refined: vol_eval of DFDet_vol.jl)
[F,f] = dfdet_contour('tratio',CP,xf,co);
if which == 1, v = F; elseif which == 2, v = f; else, v = [F(:).'; f(:).']; end
end

function dcheck(D)
if ~(isstruct(D) && isfield(D,'br') && isfield(D,'C') && isfield(D,'detrend')), error('dfdet:args','dfdet: the design must be a structure built by dfdet_determ, dfdet_trend, dfdet_perron, ...'); end
end
function pathopts(fn,opts)
for f = {'phi','hfac'}, if isfield(opts,f{1}) && ~isempty(opts.(f{1})), error('dfdet:args','%s: with rho ~= 1 the options phi, hfac and rule of the inversion on paths do not apply',fn); end, end
if isfield(opts,'rule') && ~isempty(opts.rule) && ~strcmpi(char(opts.rule),'auto'), error('dfdet:args','%s: with rho ~= 1 the options phi, hfac and rule of the inversion on paths do not apply',fn); end
end
function [start,g0,gls] = startof(opts)
g0 = getf(opts,'g0',0);
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('dfdet:args','dfdet: opts.g0 must be a finite real scalar.'); end
start = getf(opts,'start','');
if isempty(start), if g0 ~= 0, start = 'fixed'; else, start = 'zero'; end, end
start = lower(char(start));
if ~any(strcmp(start,{'zero','fixed','stationary'})), error('dfdet:args','dfdet: opts.start must be ''zero'', ''fixed'' or ''stationary''.'); end
gls = getf(opts,'gls',[]);
end
function o = pathfields(opts)
o = struct;
for f = {'verified','vtol','phi','hfac','n','ntail','alpha','filon','warn','tol','rule'}
    if isfield(opts,f{1}) && ~isempty(opts.(f{1})), o.(f{1}) = opts.(f{1}); end
end
end
function o = merge(o,p)
for f = fieldnames(p)', o.(f{1}) = p.(f{1}); end
end
function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
