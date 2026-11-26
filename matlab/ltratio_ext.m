function varargout = ltratio_ext(op,varargin)
%LTRATIO_EXT  The extensions of the limit in the classical user-level routines (opts.vol, opts.rho, c a vector).
%
%   tf = LTRATIO_EXT('has',c,opts): whether the call of LTRATIO_CDF,
%   LTRATIO_PDF or LTRATIO_ICDF has an extension of the limit: a volatility path
%   opts.vol (VOLPATH; c a number or a vector of one value per regime) or a
%   normal component opts.rho ~= 1 (the covariate-augmented t-ratio, the
%   predictive regression: CADF_C0, PREDICTIVE_C0).
%   [F,info] = LTRATIO_EXT('dist',x,cas,c,opts,mode) and [q,info] =
%   LTRATIO_EXT('icdf',p,cas,c,opts): the general engine computes them for the
%   case string or the design cas (ext_dist, ext_icdf of merged.jl), on the C
%   router where the MEX file is ahead on the path (LTRATIO_CROUTE), with the
%   MATLAB programs otherwise (DFDET_CDF, DFDET_PDF, DFDET_ICDF with opts.vol,
%   opts.rho); opts.engine = 'classical' is refused; opts.method 'auto',
%   'paths' (vol), 'lines' or 'mixture' (rho).  The quantiles start at the
%   normal approximation from the moments of the t-ratio (at the mean of the c
%   of the regimes; with rho those of rho tau + sqrt(1 - rho^2) Z).
%
%   See also LTRATIO_CDF, LTRATIO_ICDF, VOLPATH, CADF_C0.
switch op
    case 'has'
        c = varargin{1}; opts = varargin{2};
        varargout = {isstruct(opts) && (~isscalar(c) || (isfield(opts,'vol') && ~isempty(opts.vol)) || (isfield(opts,'rho') && ~isempty(opts.rho) && ~isequal(opts.rho,1)))};
    case 'dist'
        [x,cas,c,opts,mode] = varargin{:}; t0 = tic;
        [S,o,e] = spec_(cas,c,opts);
        if ltratio_croute('on')
            [spec,ev,cs] = mexargs_(cas,S,c,opts);
            F = reshape(ltratio_mex('route',mode,e,spec{:},double(x),cs,double(o.c0),method_(opts),ev),size(x));
        else
            F = dfdet_dist(x,S.D,c,o,mode);
            if strcmp(mode,'pdf'), F(F < 0) = 0; end
        end
        varargout = {F,struct('method',lower(char(getf(opts,'method','auto'))),'engine','general','ngeneral',nnz(isfinite(x)),'ndouble',0,'ncontour',0,'time',toc(t0))};
    case 'icdf'
        [p,cas,c,opts] = varargin{:};
        [S,o,e] = spec_(cas,c,opts);
        if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','ltratio_icdf: p must lie in (0,1).'); end
        tol = getf(opts,'tol',1e-10);
        if ltratio_croute('on')
            [spec,ev,cs] = mexargs_(cas,S,c,opts);
            q = reshape(ltratio_mex('route','icdf',e,spec{:},double(p),cs,double(o.c0),double(tol),method_(opts),ev),size(p));
        else
            o.tol = tol;
            q = dfdet_icdf(p,S.D,c,o);
        end
        varargout = {q,struct('engine','general')};
    otherwise
        error('ltratio:args','ltratio_ext: unknown operation %s.',op);
end
end

function [S,o,e] = spec_(cas,c,opts)
engine = lower(char(getf(opts,'engine','auto')));
if ~any(strcmp(engine,{'auto','classical','general'})), error('ltratio:args','ltratio: opts.engine must be ''auto'', ''classical'' or ''general''.'); end
if strcmp(engine,'classical'), error('ltratio:args','ltratio: engine = ''classical'': the extensions of the limit (opts.vol and opts.rho) are computed by the general engine'); end
e = 0; if strcmp(engine,'general'), e = 2; end
if ~(isnumeric(c) && isreal(c) && ~isempty(c)), error('ltratio:args','ltratio: c must be a real scalar, or with opts.vol a vector of one value per regime.'); end
S = ltratio_spec(cas,opts);
o = struct('c0',getf(opts,'c0',0),'start',S.start,'g0',S.g0,'gls',S.gls);
m = lower(char(getf(opts,'method','auto'))); if strcmp(m,'single'), m = 'line'; end
o.method = m;
for f = {'verified','vtol','vol','rho'}, if isfield(opts,f{1}) && ~isempty(opts.(f{1})), o.(f{1}) = opts.(f{1}); end, end
end

function [spec,ev,cs] = mexargs_(cas,S,c,opts)
% the computation for LTRATIO_MEX('route', ...): a case string as [model start g0], a design as 0, br, C, loc, first,
% start, g0, gls; the extensions as [verified rho nreg breaks omega c]
if ischar(cas) || isstring(cas)
    cs_ = upper(char(cas));
    if strcmp(cs_,'DG'), spec = {[4 0 0]};
    else, st = 0; if cs_(1) == 'S', st = 2; elseif S.g0 ~= 0, st = 1; end, spec = {[find('KUT' == cs_(2)) st S.g0]}; end
else
    sc = struct('zero',0,'fixed',1,'stationary',2); gls = NaN; if ~isempty(S.gls), gls = S.gls; end
    spec = [{0} ltratio_croute('design',S.D) {sc.(S.start), S.g0, gls}];
end
ev = [double(logical(getf(opts,'verified',false))) double(getf(opts,'rho',1))];
V = getf(opts,'vol',[]);
if ~isempty(V)
    if ~(isstruct(V) && isfield(V,'lam') && isfield(V,'om')), error('ltratio:args','ltratio: opts.vol must be a volatility path, volpath(breaks, omega)'); end
    k = numel(V.om);
    if ~isscalar(c) && numel(c) ~= k, error('ltratio:args','ltratio: c must be a number or a vector with one value per regime of the volatility path (%d), not %d',k,numel(c)); end
    ev = [ev k V.lam(2:end-1) V.om];
    if ~isscalar(c), ev = [ev double(c(:).')]; end
elseif ~isscalar(c)
    error('ltratio:args','ltratio: c is a vector (one value per regime): that needs a volatility path, opts.vol = volpath(breaks, omega)');
end
cs = double(c(1));
end

function m = method_(opts)
% the router's method: 0 'auto', 1 'line' ('single'), 2 'double', 3 'contour', 4 'paths', 8 'lines', 9 'mixture'
s = lower(char(getf(opts,'method','auto')));
switch s
    case 'auto', m = 0; case {'line','single'}, m = 1; case {'double','conditional'}, m = 2; case 'contour', m = 3;
    case 'paths', m = 4; case 'lines', m = 8; case 'mixture', m = 9;
    otherwise, error('ltratio:args','ltratio: opts.method must be ''auto'', ''paths'', ''lines'' or ''mixture'' with the extensions of the limit.');
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
