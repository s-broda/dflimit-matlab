function varargout = ltratio(what,varargin)
%LTRATIO  The limiting Dickey-Fuller t-ratio distributions from the C engine (MEX file LTRATIO_MEX), self-contained.
%
%   F = LTRATIO('cdf',x,cas,c,opts)    P[t <= x]
%   f = LTRATIO('pdf',x,cas,c,opts)    density
%   q = LTRATIO('icdf',p,cas,c,opts)   quantiles
%   m = LTRATIO('mom',p,cas,c,opts)    moments 1..p
%   P = LTRATIO('lincdf',b,a,cas,c,opts)   P[U - aV <= b], cases 'DK' and 'DG' (a = NaN: P[V <= b], case DK), or a design
%   P = LTRATIO('envelope',c,cas,opts)     the power envelope of Elliott, Rothenberg and Stock, cases 'DK' ('DU') and 'DG'
%   F = LTRATIO('coefcdf',y,cas,c,opts)    the coefficient statistic n(phi^ - phi0): distribution function, also
%   f = LTRATIO('coefpdf',y,cas,c,opts)    the density, the quantiles (opts.tol, 1e-9) and the moments 1..p (cas a
%   y = LTRATIO('coeficdf',p,cas,c,opts)   case or a design)
%   m = LTRATIO('coefmom',p,cas,c,opts)
%   [v,err] = LTRATIO(...)              the error estimates of the single inversion of Bao and Kan (coefcdf, coefpdf,
%                                      coeficdf, coefmom, mom, lincdf, envelope; NaN where the engines computed)
%   n = LTRATIO('threads'[,n])         the number of threads (sets it with n)
%   D = LTRATIO('design',kind,...)     a design of the general engine (any deterministic terms): kind 'trend' p,
%                                      'perron' model lam ('A', 'B', 'C', breaks lam), 'ltrend' p, 'matrix' X [TB],
%                                      'first' D (the statistic that detrends first), 'fun' f K [breaks]
%   F = LTRATIO(what,x,D,c,opts)       with a design D in place of the case: cdf, pdf, icdf, mom by the router over
%                                      the two engines (opts.engine 'auto', the default, 'classical' or 'general';
%                                      opts.start 'zero', 'fixed' (with opts.g0), 'stationary'; opts.gls = cbar for
%                                      GLS detrending; opts.method 'auto', 'line' (the line integral with the
%                                      explicit kernel: for a design that is a classical case the classical one) or
%                                      'paths' (the general engine's inversion on paths around the roots))
%   v = LTRATIO('det',fn,D,args...,opts)  the general engine's functions (fn as in LTRATIO_MEX('det', ...): coefcdf,
%                                      coefpdf, coeficdf, lincdf, linpdf, linicdf, adfcdf, adfpdf, adficdf, power,
%                                      envelope, ptcdf, ptpdf, pticdf, ptpower, ptenvelope, coefgp (coef_cdf_gp),
%                                      lmcdf, lmicdf, mzacdf, mzaicdf, msbcdf, msbicdf, and cdf, pdf, icdf, mom without
%                                      the router; the point-optimal functions and coefgp take method (0 automatic: the
%                                      single inversion of Bao and Kan up to c = 4, 5 for coefgp; 5 the inversion, 7 the
%                                      fixed grid), tol, N, umax after their arguments, and return the error estimates
%                                      as a second output)
%
%   The arguments are those of LTRATIO_CDF, LTRATIO_PDF, LTRATIO_ICDF and LTRATIO_MOM of the MATLAB programs: cas is
%   'DK', 'DU', 'DT', 'SK', 'SU', 'ST' or 'DG' (the GLS-detrended trend case: zero start, the line integral only, no
%   moments), c defaults to 0 (at least -550), opts is a structure with the fields c0 (0), g0 (0), method
%   ('auto', the default: the faster representation at each x, the contour integral with a fixed start g0 ~= 0 in
%   cases DU, DT (a point it flags gets the conditional integral's value, and the line integral's where that is
%   flagged too), the line integral at c0 = 0, with a fixed start in case DK and with c0 > 0 at x >= 0, the
%   conditional integral otherwise with c0 ~= 0, and for the values below 1e-3 of the line and the contour integral;
%   'single' or 'line', the line integral representation; 'double' or 'conditional', the conditional integral
%   representation; 'contour', the contour integral representation, a double integral with an explicit integrand,
%   for a fixed start g0 ~= 0 in cases DU, DT only)
%   and, for 'icdf', tol (1e-10).  The coefficient statistic, 'mom', 'lincdf' and 'envelope' take by default the
%   devices of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in
%   First-Order Autoregressions: O(1) Algorithms" in the six classical cases where they are checked (one inversion of the
%   closed-form characteristic function along the ray of the half-plane; the moments from one jet per node; see the
%   README): opts.method 'auto' (the default), 'inversion' ('jets' for the moments) forces them, 'engine' the engines'
%   computation of before (and for 'coefcdf', 'coeficdf' the t-ratio's representations); opts.tol their tolerance (1e-12;
%   for the moments 1e-10 relative), opts.verified recomputes by the engines the values whose estimate exceeds it.  'lincdf' and
%   'envelope' take the arguments of LTRATIO_LINCDF and POWER_ENVELOPE, with the fields cbar (-13.5, case DG of
%   'lincdf'), alpha (0.05, 'envelope'), N (60000) and umax (4000) of opts.  Nothing but
%   this file and the MEX file is needed; with the MATLAB programs on the path too, DFX_CORE, TR_LIMDIST, TR_CONTOUR
%   and TR_LIMMOM next to this file replace their engines, so that LTRATIO_CDF and the others run on the C engine.
%
%   Examples
%      ltratio('icdf', [0.01 0.05 0.10], 'DU')              % asymptotic critical values of tau_mu
%      ltratio('cdf', -2.86157, 'DU', -5)                   % limiting power of the 5% test at c = -5
%      ltratio('cdf', -2, 'DU', -5, struct('g0', 1.5))      % fixed start (the contour integral)
%      ltratio('pdf', [-1 0 1], 'DK', 2, struct('c0', 2))   % t-ratio centered at the true value
%      ltratio('mom', 4, 'SU', -5)                          % moments
%      ltratio('lincdf', -0.2639, -3.5, 'DK', 0)            % 0.05: the size of the point-optimal test at cbar = -7
%      ltratio('envelope', [-5 -10], 'DG')                  % the power envelope of the trend case
if strcmpi(what,'threads'), varargout{1} = ltratio_mex('threads', varargin{:}); return, end
if strcmpi(what,'lincdf')                             % the linear form U - aV, zero start (LTRATIO_LINCDF)
    if nargin < 4, error('ltratio:args','usage: ltratio(''lincdf'', b, a, cas, c, opts)'); end
    b = varargin{1}; a = varargin{2}; cas = varargin{3}; if ~isstruct(cas), cas = upper(char(cas)); end
    c = 0; if nargin > 4 && ~isempty(varargin{4}), c = varargin{4}; end
    if isstruct(varargin{3})                          % a design: the router (the six classical cases by the single inversion)
        D = varargin{3}; opts = struct; if nargin > 5 && ~isempty(varargin{5}), opts = varargin{5}; end
        o = detvec(rmmethod(opts));
        spec = {D.br, D.C, double(D.loc), double(D.first), nz(o(2),0), nz(o(3),0), o(4)};
        [P,E] = ltratio_mex('route', 'lincdf', engcode(opts), 0, spec{:}, double(b), double(c), double(a), NaN, imethod(opts,'lincdf'), getf(opts,'tol',NaN), getf(opts,'verified',0));
        varargout{1} = reshape(P, size(b)); varargout{2} = reshape(E, size(b)); return
    end
    if ~any(strcmp(cas,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
    opts = struct; if nargin > 5 && ~isempty(varargin{5}), opts = varargin{5}; end
    model = 1; if strcmp(cas,'DG'), model = 4; end
    E = NaN(size(b));
    if ~isfield(opts,'N') && ~isfield(opts,'umax')    % (the C router, opts.engine 'auto' by default)
        [P,E] = ltratio_mex('route', 'lincdf', engcode(opts), [model 0 0], double(b), double(c), double(a), getf(opts,'cbar',NaN), imethod(opts,'lincdf'), getf(opts,'tol',NaN), getf(opts,'verified',0));
    else
        P = ltratio_mex('lincdf', double(b), double(a), model, double(c), getf(opts,'cbar',-13.5), [getf(opts,'N',NaN) getf(opts,'umax',NaN)]);
    end
    varargout{1} = reshape(P, size(b)); varargout{2} = reshape(E, size(b)); return
end
if strcmpi(what,'envelope')                           % the power envelope of Elliott, Rothenberg and Stock (POWER_ENVELOPE)
    if nargin < 2, error('ltratio:args','usage: ltratio(''envelope'', c, cas, opts)'); end
    c = varargin{1}; cas = 'DK'; if nargin > 2 && ~isempty(varargin{2}), cas = upper(char(varargin{2})); end
    if strcmp(cas,'DU'), cas = 'DK'; end
    if ~any(strcmp(cas,{'DK','DG'})), error('ltratio:args','power_envelope: the case must be DK (or DU) or DG.'); end
    opts = struct; if nargin > 3 && ~isempty(varargin{3}), opts = varargin{3}; end
    model = 1; if strcmp(cas,'DG'), model = 4; end
    E = NaN(size(c));
    if ~isfield(opts,'N') && ~isfield(opts,'umax')    % (the C router, opts.engine 'auto' by default)
        [P,E] = ltratio_mex('route', 'envelope', engcode(opts), [model 0 0], double(c), getf(opts,'alpha',0.05), imethod(opts,'envelope'), getf(opts,'tol',NaN), getf(opts,'verified',0));
    else
        P = ltratio_mex('envelope', double(c), model, getf(opts,'alpha',0.05), [getf(opts,'N',NaN) getf(opts,'umax',NaN)]);
    end
    varargout{1} = reshape(P, size(c)); varargout{2} = reshape(E, size(c)); return
end
if any(strcmpi(what,{'coefcdf','coefpdf','coeficdf','coefmom'}))   % the coefficient statistic by the router
    if nargin < 3, error('ltratio:args','usage: ltratio(''%s'', y, cas, c, opts)',lower(what)); end
    x = varargin{1}; cas = varargin{2};
    c = 0; if nargin > 3 && ~isempty(varargin{3}), c = varargin{3}; end
    opts = struct; if nargin > 4 && ~isempty(varargin{4}), opts = varargin{4}; end
    o = detvec(rmmethod(opts)); eng = engcode(opts); c0 = nz(o(1),0);
    if isstruct(cas)
        spec = {0, cas.br, cas.C, double(cas.loc), double(cas.first), nz(o(2),0), nz(o(3),0), o(4)};
    else
        [model,stat] = casecode(cas); g0 = nz(o(3),0); st = 0; if stat, st = 2; elseif g0 ~= 0, st = 1; end
        spec = {[model st g0]};
    end
    im = imethod(opts,lower(what)); vf = getf(opts,'verified',0);
    switch lower(what)
        case 'coefcdf', [v,E] = ltratio_mex('route', 'coefcdf', eng, spec{:}, double(x), double(c), double(c0), im, getf(opts,'tol',NaN), vf);
        case 'coefpdf', [v,E] = ltratio_mex('route', 'coefpdf', eng, spec{:}, double(x), double(c), double(c0), im, getf(opts,'tol',NaN), vf);
        case 'coeficdf', [v,E] = ltratio_mex('route', 'coeficdf', eng, spec{:}, double(x), double(c), double(c0), getf(opts,'tol',1e-9), im, vf);
        case 'coefmom', [v,E] = ltratio_mex('route', 'coefmom', eng, spec{:}, double(x), double(c), double(c0), im, getf(opts,'tol',NaN));
    end
    if ~strcmpi(what,'coefmom'), v = reshape(v, size(x)); E = reshape(E, size(x)); end
    varargout{1} = v; varargout{2} = E; return
end
if strcmpi(what,'design')                             % a design of the general engine (dfl_design)
    if nargin < 2, error('ltratio:args','usage: D = ltratio(''design'', kind, ...)'); end
    kind = lower(char(varargin{1})); a = varargin(2:end);
    switch kind
        case {'trend','ltrend'}, [br,C,loc,first] = ltratio_mex('design', kind, double(a{1}));
        case 'perron'
            m = find(upper(char(a{1})) == 'ABC'); if isempty(m), error('ltratio:args','perron: model A, B or C.'); end
            [br,C,loc,first] = ltratio_mex('design', 'perron', m, double(a{2}));
        case 'matrix', TB = []; if numel(a) > 1, TB = double(a{2}); end; [br,C,loc,first] = ltratio_mex('design', 'matrix', double(a{1}), TB);
        case 'first', D = a{1}; [br,C,loc,first] = ltratio_mex('design', 'first', D.br, D.C, D.loc, D.first);
        case 'fun', b = []; if numel(a) > 2, b = double(a{3}); end; [br,C,loc,first] = ltratio_mex('design', 'fun', a{1}, double(a{2}), b);
        otherwise, error('ltratio:args','the kind of design must be trend, perron, ltrend, matrix, first or fun.');
    end
    varargout{1} = struct('br',br,'C',C,'loc',loc,'first',first); return
end
if strcmpi(what,'det')                                % a function of the general engine with a design
    if nargin < 3, error('ltratio:args','usage: ltratio(''det'', fn, D, args..., opts)'); end
    fn = lower(char(varargin{1})); D = varargin{2}; args = varargin(3:end);
    o = NaN(1,5);
    if ~isempty(args) && isstruct(args{end}), o = detvec(args{end}); args = args(1:end-1); end
    nexp = struct('cdf',2,'pdf',2,'icdf',4,'coefcdf',2,'coefpdf',2,'coeficdf',4,'mom',2,'lincdf',3,'linpdf',3,'linicdf',5, ...
        'adfcdf',2,'adfpdf',2,'adficdf',4,'power',3,'envelope',2,'ptcdf',7,'ptpdf',7,'pticdf',6,'ptpower',7,'ptenvelope',6,'coefgp',6, ...
        'lmcdf',2,'lmicdf',2,'mzacdf',3,'mzaicdf',3,'msbcdf',3,'msbicdf',3);
    if ~isfield(nexp,fn), error('ltratio:args','unknown function %s of the general engine.',fn); end
    for k = 1:numel(args), args{k} = double(args{k}); end
    args(end+1:nexp.(fn)) = {[]};
    [v,E] = ltratio_mex('det', fn, D.br, D.C, double(D.loc), double(D.first), args{:}, o);
    if ~isempty(args) && ~any(strcmp(fn,{'mom','ptpower'})), v = reshape(v, size(args{1})); if numel(E) == numel(v), E = reshape(E, size(v)); end, end
    varargout{1} = v; varargout{2} = E; return
end
if nargin < 3, error('ltratio:args','usage: ltratio(what, x, cas, c, opts)'); end
if isstruct(varargin{2})                              % a design: the router over the two engines
    x = varargin{1}; D = varargin{2};
    c = 0; if nargin > 3 && ~isempty(varargin{3}), c = varargin{3}; end
    opts = struct; if nargin > 4 && ~isempty(varargin{4}), opts = varargin{4}; end
    if strcmpi(what,'mom'), o = detvec(rmmethod(opts)); else, o = detvec(opts); end
    eng = engcode(opts);
    c0 = o(1); if isnan(c0), c0 = 0; end
    st = o(2); if isnan(st), st = 0; end
    g0 = o(3); if isnan(g0), g0 = 0; end
    spec = {D.br, D.C, double(D.loc), double(D.first), st, g0, o(4)};
    rm = 0; if ~isnan(o(6)), rm = [0 1 4]; rm = rm(o(6) + 1); end   % (the router's codes: 1 the line integral, 4 the general engine's paths)
    switch lower(what)
        case {'cdf','pdf'}, v = ltratio_mex('route', lower(what), eng, 0, spec{:}, double(x), double(c), double(c0), rm);
        case 'icdf', tol = getf(opts,'tol',1e-10); v = ltratio_mex('route', 'icdf', eng, 0, spec{:}, double(x), double(c), double(c0), tol, rm);
        case 'mom', [v,E] = ltratio_mex('route', 'mom', eng, 0, spec{:}, double(x), double(c), double(c0), imethod(opts,'mom'), getf(opts,'tol',NaN), getf(opts,'verified',0)); varargout{2} = E;
        otherwise, error('ltratio:args','with a design, what must be cdf, pdf, icdf or mom (the other functions: ltratio(''det'', ...)).');
    end
    if ~strcmpi(what,'mom'), v = reshape(v, size(x)); end
    varargout{1} = v; return
end
x = varargin{1}; cas = upper(char(varargin{2}));
if strcmp(cas,'DG')                                   % the GLS-detrended trend case (model 4): zero start, the line integral only
    model = 4; stat = 0;
elseif numel(cas) == 2 && any(cas(1) == 'DS') && any(cas(2) == 'KUT')
    model = find('KUT' == cas(2)); stat = double(cas(1) == 'S');
else
    error('ltratio:args','the case must be DK, DU, DT, SK, SU, ST or DG.');
end
c = 0; if nargin > 3 && ~isempty(varargin{3}), c = varargin{3}; end
opts = struct; if nargin > 4 && ~isempty(varargin{4}), opts = varargin{4}; end
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
g0 = 0; if isfield(opts,'g0') && ~isempty(opts.g0), g0 = opts.g0; end
method = 0;
if isfield(opts,'method') && ~isempty(opts.method) && ~strcmpi(what,'mom')
    switch lower(char(opts.method))
        case 'auto', method = 0;
        case {'single','line'}, method = 1;                 % the line integral representation
        case {'double','conditional'}, method = 2;          % the conditional integral representation
        case 'contour', method = 3;                         % the contour integral representation (a fixed start in cases DU, DT)
        otherwise, error('ltratio:args','opts.method must be ''auto'', ''single'' (''line''), ''double'' (''conditional'') or ''contour''.');
    end
end
verified = isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified);
if strcmpi(what,'mom')                                % the moments by the router: the jets by default, the circles with 'engine'
    st = 0; if stat, st = 2; elseif g0 ~= 0, st = 1; end
    [v,E] = ltratio_mex('route', 'mom', engcode(opts), [model st g0], double(x), double(c), double(c0), imethod(opts,'mom'), getf(opts,'tol',NaN), double(verified));
    varargout{1} = v; varargout{2} = E; return
end
if any(strcmpi(what,{'cdf','pdf','icdf'}))           % the router over the two engines (opts.engine 'auto' by default)
    eng = engcode(opts); st = 0; if stat, st = 2; elseif g0 ~= 0, st = 1; end
    if eng ~= 1 && (method ~= 0 || verified), eng = 1; end   % (a representation named, or the checked computation: the classical engine)
    if eng == 1
        % (the classical engine: opts.engine = 'classical', a method named or the checked computation)
    else
        switch lower(what)
            case {'cdf','pdf'}, v = ltratio_mex('route', lower(what), eng, [model st g0], double(x), double(c), double(c0));
            case 'icdf', v = ltratio_mex('route', 'icdf', eng, [model st g0], double(x), double(c), double(c0), getf(opts,'tol',1e-10));
        end
        if ~strcmpi(what,'mom'), v = reshape(v, size(x)); end
        varargout{1} = v; return
    end
end
if ~verified && any(strcmpi(what,{'cdf','pdf','icdf'}))   % the warning of a brittle regime of the one-evaluation default
    bm = ltratio_mex('brittle', 2, model, double(c), double(c0), double(g0), stat, method);
    if ~isempty(bm), warning('ltratio:brittle','%s',bm); end
end
switch lower(what)
    case 'cdf'
        varargout{1} = ltratio_mex('dist', 2, double(x), model, double(c), double(c0), double(g0), stat, method, double(verified));
    case 'pdf'
        varargout{1} = ltratio_mex('dist', 1, double(x), model, double(c), double(c0), double(g0), stat, method, double(verified));
    case 'icdf'
        tol = 1e-10; if isfield(opts,'tol') && ~isempty(opts.tol), tol = opts.tol; end
        [q,Fq,it,conv] = ltratio_mex('icdf', double(x), model, double(c), double(c0), double(g0), stat, method, double(tol), double(verified));
        varargout{1} = q; if nargout > 1, varargout{2} = struct('F',Fq,'iter',it,'conv',logical(conv)); end
        if ~all(conv(:)), warning('ltratio:icdf','ltratio: Newton''s method did not converge at %d of %d points.',nnz(~conv),numel(conv)); end
    case 'mom'
        varargout{1} = ltratio_mex('mom', double(x), model, double(c), double(c0), double(g0), stat);
    otherwise
        error('ltratio:args','what must be cdf, pdf, icdf, mom, lincdf, envelope or threads.');
end
end
function v = getf(o,f,d), if isfield(o,f) && ~isempty(o.(f)), v = double(o.(f)); else, v = d; end, end
function o = detvec(opts)                             % [c0 start g0 gls verified method] of the general engine (NaN: the default)
o = NaN(1,6);
if isfield(opts,'c0') && ~isempty(opts.c0), o(1) = opts.c0; end
if isfield(opts,'g0') && ~isempty(opts.g0), o(3) = opts.g0; if opts.g0 ~= 0, o(2) = 1; end, end
if isfield(opts,'start') && ~isempty(opts.start)
    s = opts.start;
    if ischar(s) || isstring(s)
        k = find(strcmpi(char(s),{'zero','fixed','stationary'})); if isempty(k), error('ltratio:args','opts.start must be ''zero'', ''fixed'' or ''stationary''.'); end
        o(2) = k - 1;
    else, o(2) = s; end
end
if isfield(opts,'gls') && ~isempty(opts.gls), o(4) = opts.gls; end
if isfield(opts,'verified') && ~isempty(opts.verified), o(5) = double(logical(opts.verified)); end
if isfield(opts,'method') && ~isempty(opts.method), o(6) = gmethod(opts.method); end
end
function m = gmethod(s)                               % opts.method of the general engine: 0 'auto', 1 'line', 2 'paths'
k = find(strcmpi(char(s),{'auto','line','paths'}));
if isempty(k), error('ltratio:args','with a design opts.method must be ''auto'', ''line'' (the line integral with the explicit kernel) or ''paths'' (the inversion on paths around the roots).'); end
m = k - 1;
end
function e = engcode(opts)                            % opts.engine: 0 auto, 1 classical, 2 general
e = 0;
if isfield(opts,'engine') && ~isempty(opts.engine)
    k = find(strcmpi(char(opts.engine),{'auto','classical','general'})); if isempty(k), error('ltratio:args','opts.engine must be ''auto'', ''classical'' or ''general''.'); end
    e = k - 1;
end
end
function v = nz(x,d), if isnan(x), v = d; else, v = x; end, end
function o = rmmethod(o), if isfield(o,'method'), o = rmfield(o,'method'); end, end
function [model,stat] = casecode(cas)                 % the model and the start of a case string
cas = upper(char(cas)); stat = 0;
if strcmp(cas,'DG'), model = 4;
elseif numel(cas) == 2 && any(cas(1) == 'DS') && any(cas(2) == 'KUT'), model = find('KUT' == cas(2)); stat = double(cas(1) == 'S');
else, error('ltratio:args','the case must be DK, DU, DT, SK, SU, ST or DG.'); end
end
function m = imethod(opts,what)                       % opts.method of the functions with the devices of Bao and Kan (dfl_inv_opts)
m = 0;
if ~isfield(opts,'method') || isempty(opts.method), return, end
s = lower(char(opts.method)); rep = any(strcmp(what,{'coefcdf','coeficdf'}));
switch s
    case 'auto', m = 0;
    case {'inversion','jets'}, m = 5;
    case 'engine', m = 6;
    case {'single','line'}, if rep, m = 1; else, m = -1; end
    case {'double','conditional'}, if rep, m = 2; else, m = -1; end
    case 'contour', if rep, m = 3; else, m = -1; end
    case 'paths', if rep, m = 4; else, m = -1; end
    otherwise, m = -1;
end
if m < 0
    extra = ''; if rep, extra = ', or a representation of the t-ratio (''line'', ''conditional'', ''contour'', ''paths'')'; end
    error('ltratio:args','%s: opts.method must be ''auto'', ''inversion'' (the single inversion of Bao and Kan; ''jets'' for the moments) or ''engine'' (the engines'' computation)%s.',what,extra);
end
end
