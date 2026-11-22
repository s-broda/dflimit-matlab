function varargout = ltratio_croute(op,varargin)
%LTRATIO_CROUTE  The C router of the two engines (the MEX file LTRATIO_MEX) under the user-level routines.
%
%   tf = LTRATIO_CROUTE('on'): whether the folder of the MEX file of the C
%   engine (c/mex in the repository, the folder of a downloaded zip: the MEX
%   file with its drop-ins DFX_CORE, TR_LIMDIST, ...) is ahead of these
%   programs on the path (LTRATIO_MEX exists and DFX_CORE resolves to its
%   folder).  Then LTRATIO_CDF, LTRATIO_PDF, LTRATIO_ICDF, LTRATIO_MOM,
%   COEF_CDF, COEF_ICDF, ADFMAX_CDF, ADFMAX_PDF, LTRATIO_LINCDF,
%   POWER_ENVELOPE and LTRATIO_REGION hand a call with a design, or with a
%   case string and opts.engine, without options of the MATLAB engines to
%   the C router (a case string without opts.engine keeps the classical
%   programs, on the C engine through the drop-ins, as LTRATIO of the MEX
%   folder and the Python package do) (dfl_engine_* of c/dfl_det_route.c): its
%   choice of engine (dfl_autoengine, set by the timings of the two C
%   engines) and both engines in C, with the warnings of the C engine
%   (ltratio:brittle where the router takes the classical engine in a regime
%   in which one evaluation is brittle, ltratio:bound, ltratio:det).  Such a
%   call has in opts at most the fields c0, g0, start, gls, engine ('auto'
%   or 'general'), method ('auto'), verified (false) and those of the
%   function (tol for the quantiles, cbar for the linear form, alpha for the
%   envelope).  With engine = 'classical', a representation named, the
%   checked computation or further options the routines compute as without
%   the router: the classical programs on the C engine through the drop-ins,
%   the general engine in MATLAB.  The checks of the arguments are those of
%   the MATLAB routines (the same errors), and the second output (info) has
%   the field engine = 'C router' (and the fields the function names below;
%   the counts of the representations of LTRATIO_CDF, ndouble, ncontour,
%   ngeneral, are NaN: the C router does not report them).
%
%   tf = LTRATIO_CROUTE('use',cas,opts,extra): whether the call goes to the
%   C router (on, a design or opts.engine, plain).
%   tf = LTRATIO_CROUTE('plain',opts,extra): whether opts has only the
%   fields of such a call (and those in the cell extra).
%   [F,info] = LTRATIO_CROUTE('dist',x,cas,c,opts,mode), mode 'cdf' or 'pdf';
%   [q,info] = LTRATIO_CROUTE('icdf',p,cas,c,opts) (info.F, info.conv: the
%   distribution function at q, |F - p| < 1e-6); [m,info] =
%   LTRATIO_CROUTE('mom',p,cas,c,opts); F = LTRATIO_CROUTE('coef',y,cas,c,
%   opts); y = LTRATIO_CROUTE('coeficdf',p,cas,c,opts); [F,info] =
%   LTRATIO_CROUTE('adfmax',x,cas,c,opts,mode); P = LTRATIO_CROUTE('lincdf',
%   b,a,cas,c,opts); P = LTRATIO_CROUTE('envelope',c,cas,opts); [P,info] =
%   LTRATIO_CROUTE('region',bnd,cas,c,opts) (tf = LTRATIO_CROUTE('bnd',bnd):
%   the boundaries of the C router, types line, parabola and wap without v0,
%   brks, g).
%   d = LTRATIO_CROUTE('design',D): {br, C, loc, first} of a design for the
%   MEX file (the coefficient matrices of the pieces side by side, padded to
%   a common degree).
%
%   See also LTRATIO_ENGINE, LTRATIO_ROUTE, LTRATIO_CDF.
switch op
    case 'on'
        m = which('ltratio_mex');
        tf = ~isempty(m) && exist('ltratio_mex','file') == 3;
        if tf, d = which('dfx_core'); tf = ~isempty(d) && strcmp(fileparts(d),fileparts(m)); end
        varargout = {tf};
    case 'use'
        % (a design, or a case string with opts.engine: a case string without it keeps the classical programs, as
        % LTRATIO of the MEX folder and the Python package do)
        cas = varargin{1}; opts = varargin{2}; extra = {}; if numel(varargin) > 2, extra = varargin{3}; end
        named = isstruct(opts) && isfield(opts,'engine') && ~isempty(opts.engine);
        varargout = {(~(ischar(cas) || isstring(cas)) || named) && ltratio_croute('plain',opts,extra) && ltratio_croute('on')};
    case 'plain'
        opts = varargin{1}; extra = {}; if numel(varargin) > 1, extra = varargin{2}; end
        if ~isstruct(opts), varargout = {false}; return, end
        ok = all(ismember(fieldnames(opts),[{'c0','g0','start','gls','engine','method','verified'} extra]));
        if ok && isfield(opts,'method') && ~isempty(opts.method), ok = (ischar(opts.method) || isstring(opts.method)) && strcmpi(char(opts.method),'auto'); end
        if ok && isfield(opts,'verified') && ~isempty(opts.verified), ok = ~logical(opts.verified); end
        if ok && isfield(opts,'engine') && ~isempty(opts.engine), ok = (ischar(opts.engine) || isstring(opts.engine)) && any(strcmpi(char(opts.engine),{'auto','general'})); end
        varargout = {ok};
    case 'design'
        varargout = {design_(varargin{1})};
    case 'bnd'
        b = varargin{1};
        varargout = {isstruct(b) && isfield(b,'type') && any(strcmpi(char(b.type),{'line','parabola','wap'})) && ~any(isfield(b,{'v0','brks','g'}))};
    case 'dist'
        [x,cas,c,opts,mode] = varargin{:}; t0 = tic;
        [spec,e,c0] = spec_(cas,c,opts);
        if ~(isnumeric(x) && isreal(x)), error('ltratio:args','ltratio_%s: x must be real.',mode); end
        F = reshape(ltratio_mex('route',mode,e,spec{:},double(x),double(c),double(c0)),size(x));
        varargout = {F,struct('method','auto','engine','C router','ndouble',NaN,'ncontour',NaN,'ngeneral',NaN,'time',toc(t0))};
    case 'icdf'
        [p,cas,c,opts] = varargin{:};
        [spec,e,c0] = spec_(cas,c,rmf_(opts,'tol'));
        if ~(isnumeric(p) && isreal(p) && all(p(:) > 0 & p(:) < 1)), error('ltratio:args','ltratio_icdf: p must lie in (0,1).'); end
        tol = getf_(opts,'tol',1e-10);
        q = reshape(ltratio_mex('route','icdf',e,spec{:},double(p),double(c),double(c0),double(tol)),size(p));
        info = struct('engine','C router');
        if nargout > 1
            ws = warning('off','ltratio:brittle');       % (the regime was named by the quantile search)
            Fq = reshape(ltratio_mex('route','cdf',e,spec{:},q,double(c),double(c0)),size(p));
            warning(ws);
            info = struct('F',Fq,'iter',NaN(size(p)),'conv',abs(Fq - p) < 1e-6,'engine','C router');
        end
        varargout = {q,info};
    case 'mom'
        [p,cas,c,opts] = varargin{:};
        [spec,e,c0] = spec_(cas,c,opts);
        if ~(isnumeric(p) && isscalar(p) && p == round(p) && p >= 1), error('ltratio:args','ltratio_mom: p must be a positive integer.'); end
        m = ltratio_mex('route','mom',e,spec{:},double(p),double(c),double(c0));
        varargout = {m,struct('err',NaN(size(m)),'engine','C router')};
    case {'coef','coeficdf'}
        [y,cas,c,opts] = varargin{:};
        [spec,e,c0] = spec_(cas,c,rmf_(opts,'tol'));
        if strcmp(op,'coef')
            if ~(isnumeric(y) && isreal(y)), error('ltratio:args','coef_cdf: y must be real.'); end
            v = ltratio_mex('route','coefcdf',e,spec{:},double(y),double(c),double(c0));
        else
            if ~(isnumeric(y) && isreal(y) && all(y(:) > 0 & y(:) < 1)), error('ltratio:args','coef_icdf: p must lie in (0,1).'); end
            v = ltratio_mex('route','coeficdf',e,spec{:},double(y),double(c),double(c0),getf_(opts,'tol',1e-9));
        end
        varargout = {reshape(v,size(y))};
    case 'adfmax'
        [x,cas,c,opts,mode] = varargin{:};
        [spec,e] = spec_(cas,c,opts,'adfmax');
        if ~(isnumeric(x) && isreal(x)), error('ltratio:args','adfmax: x must be real.'); end
        F = reshape(ltratio_mex('route',['adf' mode],e,spec{:},double(x),double(c)),size(x));
        varargout = {F,struct('engine','C router')};
    case 'lincdf'
        [b,a,cas,c,opts] = varargin{:};
        if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_lincdf: c must be a finite real scalar.'); end
        if c < -550, error('ltratio:args','ltratio_lincdf: c must be at least -550 (the range in which the programs are checked; the moments fail from about c = -600, the line integral from about c = -650).'); end
        if ~(isnumeric(a) && isscalar(a) && isreal(a)), error('ltratio:args','ltratio_lincdf: a must be a real scalar (NaN: the distribution of V).'); end
        if ~(isnumeric(b) && isreal(b)), error('ltratio:args','ltratio_lincdf: b must be real.'); end
        e = engine_(opts); cb = getf_(opts,'cbar',NaN);
        if ischar(cas) || isstring(cas)
            cs = upper(char(cas)); if strcmp(cs,'DU'), cs = 'DK'; end
            if ~any(strcmp(cs,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
        else
            S = ltratio_spec(cas,opts); cs = '';
            if strcmp(S.start,'zero') && any(strcmp(S.case,{'DK','DG'})), cs = S.case; cb = NaN; end
        end
        if ~isempty(cs)
            P = ltratio_mex('route','lincdf',e,[1 + 3*strcmp(cs,'DG') 0 0],double(b),double(c),double(a),double(cb));
        else
            % (a design that the router's linear form, which takes the cases DK and DG, does not have: the general engine in C)
            d = design_(S.D);
            P = ltratio_mex('det','lincdf',d{:},double(b),double(a),double(c),detopts_(S,0));
        end
        varargout = {reshape(P,size(b))};
    case 'envelope'
        [c,cas,opts] = varargin{:};
        e = engine_(opts); alpha = getf_(opts,'alpha',0.05);
        if ~(isnumeric(alpha) && isscalar(alpha) && alpha > 0 && alpha < 1), error('ltratio:args','power_envelope: the level opts.alpha must lie in (0, 1).'); end
        if ~(isnumeric(c) && isreal(c)), error('ltratio:args','power_envelope: c must be real.'); end
        if ischar(cas) || isstring(cas)
            cs = upper(char(cas)); if strcmp(cs,'DU'), cs = 'DK'; end
            if ~any(strcmp(cs,{'DK','DG'})), error('ltratio:args','power_envelope: the case must be DK (or DU) or DG.'); end
        else
            cs = '';
            for q = -1:1
                T = dfdet_trend(q); D = cas;
                if isequal(D.br,T.br) && numel(D.C) == 1 && isequal(size(D.C{1}),size(T.C{1})) && isequal(D.C{1},T.C{1}) && ~D.loc && strcmp(D.detrend,'onestep')
                    if q < 1, cs = 'DK'; else, cs = 'DG'; end
                end
            end
        end
        if ~isempty(cs)
            P = ltratio_mex('route','envelope',e,[1 + 3*strcmp(cs,'DG') 0 0],double(c),double(alpha));
        else
            % (another design: the envelope of the general construction, the point-optimal statistic against cbar = c, in C)
            d = design_(cas);
            P = ltratio_mex('det','ptenvelope',d{:},double(c),double(alpha));
        end
        varargout = {reshape(P,size(c))};
    case 'region'
        [bnd,cas,c,opts] = varargin{:};
        if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio_region: c must be a finite real scalar.'); end
        if c < -550, error('ltratio:args','ltratio_region: c must be at least -550 (the range in which the programs are checked).'); end
        e = engine_(opts);
        if ischar(cas) || isstring(cas)
            if ~strcmpi(cas,'DK'), error('ltratio:args','ltratio_region: the case must be DK (zero start).'); end
            spec = {[1 0 0]};
        else
            spec = spec_(cas,c,opts);
        end
        switch lower(char(bnd.type))
            case 'line', ty = 1; p1 = bnd.a; p2 = bnd.b;
            case 'parabola', ty = 2; p1 = bnd.x; p2 = 0; if isfield(bnd,'c0') && ~isempty(bnd.c0), p2 = bnd.c0; end
            otherwise, ty = 3; p1 = bnd.k; p2 = bnd.omega;
        end
        P = ltratio_mex('route','region',e,spec{:},ty,double(c),double(p1),double(p2));
        varargout = {P,struct('engine','C router')};
    otherwise
        error('ltratio:args','ltratio_croute: unknown operation %s.',op);
end
end

function [spec,e,c0] = spec_(cas,c,opts,what)
% the computation for LTRATIO_MEX('route', ...): a case string (after the checks of LTRATIO_ARGS), or a design that is a
% classical case the classical programs accept, as [model start g0]; another design as 0, br, C, loc, first, start, g0, gls
if nargin < 4, what = ''; end
e = engine_(opts);
oc = rmf_(rmf_(rmf_(opts,'engine'),'start'),'gls');
c0 = getf_(opts,'c0',0); g0 = getf_(opts,'g0',0);
if ischar(cas) || isstring(cas)
    cs = upper(char(cas));
    if isempty(what), ltratio_args(cs,c,oc); end        % (the errors of the MATLAB routines; ADFmax: adfmax_dist checks)
    if strcmp(cs,'DG'), spec = {[4 0 0]}; return, end
    if ~(numel(cs) == 2 && any(cs(1) == 'DS') && any(cs(2) == 'KUT')), error('ltratio:args','ltratio: the case must be one of DK, DU, DT, SK, SU, ST, DG.'); end
    st = 0; if cs(1) == 'S', st = 2; elseif g0 ~= 0, st = 1; end
    spec = {[find('KUT' == cs(2)) st g0]};
    return
end
S = ltratio_spec(cas,opts);
sc = struct('zero',0,'fixed',1,'stationary',2);
if ~isempty(S.case) && ltratio_route('valid',S,c,oc)
    spec = spec_(S.case,c,setfield(oc,'g0',S.g0),'checked'); %#ok<SFLD>
    return
end
d = design_(S.D); gls = NaN; if ~isempty(S.gls), gls = S.gls; end
spec = [{0} d {sc.(S.start), S.g0, gls}];
end

function d = design_(D)
% br, the coefficient matrices side by side (padded to a common degree), loc, first
np = numel(D.C); dg = max(cellfun(@(M) size(M,1),D.C)); K = size(D.C{1},2);
C = zeros(dg,K*np);
for j = 1:np, M = D.C{j}; C(1:size(M,1),(j-1)*K+1:j*K) = M; end
if K == 0, C = zeros(1,0); end
d = {double(D.br(:)'), C, double(logical(D.loc)), double(strcmp(D.detrend,'first'))};
end

function o = detopts_(S,c0)
% [c0 start g0 gls verified] of the general engine of the MEX file
sc = struct('zero',0,'fixed',1,'stationary',2); gls = NaN; if ~isempty(S.gls), gls = S.gls; end
o = [c0 sc.(S.start) S.g0 gls 0];
end

function e = engine_(opts)
e = 0;
if isfield(opts,'engine') && ~isempty(opts.engine) && strcmpi(char(opts.engine),'general'), e = 2; end
end

function o = rmf_(o,f)
if isfield(o,f), o = rmfield(o,f); end
end

function v = getf_(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = double(o.(f)); else, v = d; end
end
