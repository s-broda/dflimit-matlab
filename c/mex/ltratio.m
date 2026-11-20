function varargout = ltratio(what,varargin)
%LTRATIO  The limiting Dickey-Fuller t-ratio distributions from the C engine (MEX file LTRATIO_MEX), self-contained.
%
%   F = LTRATIO('cdf',x,cas,c,opts)    P[t <= x]
%   f = LTRATIO('pdf',x,cas,c,opts)    density
%   q = LTRATIO('icdf',p,cas,c,opts)   quantiles
%   m = LTRATIO('mom',p,cas,c,opts)    moments 1..p
%   P = LTRATIO('lincdf',b,a,cas,c,opts)   P[U - aV <= b], cases 'DK' and 'DG' (a = NaN: P[V <= b], case DK)
%   P = LTRATIO('envelope',c,cas,opts)     the power envelope of Elliott, Rothenberg and Stock, cases 'DK' ('DU') and 'DG'
%   n = LTRATIO('threads'[,n])         the number of threads (sets it with n)
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
%   and, for 'icdf', tol (1e-10).  'lincdf' and
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
    b = varargin{1}; a = varargin{2}; cas = upper(char(varargin{3}));
    if ~any(strcmp(cas,{'DK','DG'})), error('ltratio:args','ltratio_lincdf: the case must be DK or DG.'); end
    c = 0; if nargin > 4 && ~isempty(varargin{4}), c = varargin{4}; end
    opts = struct; if nargin > 5 && ~isempty(varargin{5}), opts = varargin{5}; end
    model = 1; if strcmp(cas,'DG'), model = 4; end
    P = ltratio_mex('lincdf', double(b), double(a), model, double(c), getf(opts,'cbar',-13.5), [getf(opts,'N',NaN) getf(opts,'umax',NaN)]);
    varargout{1} = reshape(P, size(b)); return
end
if strcmpi(what,'envelope')                           % the power envelope of Elliott, Rothenberg and Stock (POWER_ENVELOPE)
    if nargin < 2, error('ltratio:args','usage: ltratio(''envelope'', c, cas, opts)'); end
    c = varargin{1}; cas = 'DK'; if nargin > 2 && ~isempty(varargin{2}), cas = upper(char(varargin{2})); end
    if strcmp(cas,'DU'), cas = 'DK'; end
    if ~any(strcmp(cas,{'DK','DG'})), error('ltratio:args','power_envelope: the case must be DK (or DU) or DG.'); end
    opts = struct; if nargin > 3 && ~isempty(varargin{3}), opts = varargin{3}; end
    model = 1; if strcmp(cas,'DG'), model = 4; end
    P = ltratio_mex('envelope', double(c), model, getf(opts,'alpha',0.05), [getf(opts,'N',NaN) getf(opts,'umax',NaN)]);
    varargout{1} = reshape(P, size(c)); return
end
if nargin < 3, error('ltratio:args','usage: ltratio(what, x, cas, c, opts)'); end
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
if isfield(opts,'method') && ~isempty(opts.method)
    switch lower(char(opts.method))
        case 'auto', method = 0;
        case {'single','line'}, method = 1;                 % the line integral representation
        case {'double','conditional'}, method = 2;          % the conditional integral representation
        case 'contour', method = 3;                         % the contour integral representation (a fixed start in cases DU, DT)
        otherwise, error('ltratio:args','opts.method must be ''auto'', ''single'' (''line''), ''double'' (''conditional'') or ''contour''.');
    end
end
verified = isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified);
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
