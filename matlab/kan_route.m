function varargout = kan_route(op,varargin)
%KAN_ROUTE  The rules of the single inversion of Bao and Kan in the user-level routines (inv_take of DFLimit.jl).
%
%   The coefficient statistic (COEF_CDF, COEF_ICDF, COEF_PDF, COEF_MOM), the linear form (LTRATIO_LINCDF), the power
%   envelope (POWER_ENVELOPE) and the moments of the t-ratio (LTRATIO_MOM) take by default the devices of Yong Bao and
%   Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in First-Order
%   Autoregressions: O(1) Algorithms" (KAN_CORE) in the six classical cases where they were compared with the engines
%   (merged.jl of DFLimit.jl and the C router: the same rules and constants).
%
%   m = KAN_ROUTE('method',opts,fn,extra): opts.method of the function fn: 'auto' (the default), 'inversion' (forces the
%   single inversion; 'jets' the same for the moments), 'engine' (the engines' computation of before), and where extra
%   the representations of the t-ratio ('single' or 'line', 'double' or 'conditional', 'contour'), which take the
%   engines too.
%   tf = KAN_ROUTE('take',fn,m,engine,cas,c,plain): whether the single inversion computes: with m = 'inversion' always
%   (an error with engine 'general', for a case it does not have, or with options of the engines), with m = 'auto' under
%   engine 'auto' in the cases DK, DU, DT, SK, SU, ST with plain options for c up to 20 (INV_CMAX); otherwise the engines.
%   tf = KAN_ROUTE('plain',opts,extra): opts has no options of the engines (only c0, g0, method, verified, tol, engine,
%   start, gls and the fields in the cell extra).
%   tol = KAN_ROUTE('tol',opts,default): opts.tol (a positive finite number), or the default.
%   KAN_ROUTE('warn',fn,err,T): the warning ltratio:inversion where the estimates err exceed their tolerances T (a
%   number, or one per point).
%   [model,stat] = KAN_ROUTE('model',cas): the model (1, 2, 3) and the start of a case string.
%
%   See also KAN_CORE, COEF_CDF, LTRATIO_LINCDF, POWER_ENVELOPE, LTRATIO_MOM.
switch op
    case 'method'
        opts = varargin{1}; fn = varargin{2}; extra = numel(varargin) > 2 && varargin{3};
        m = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), m = lower(char(opts.method)); end
        if strcmp(m,'jets'), m = 'inversion'; end
        if strcmp(m,'line'), m = 'single'; elseif strcmp(m,'conditional'), m = 'double'; end
        ok = any(strcmp(m,{'auto','inversion','engine'})) || (extra && any(strcmp(m,{'single','double','contour'})));
        if ~ok
            x = ''; if extra, x = ', or a representation of the t-ratio (''single'' or ''line'', ''double'' or ''conditional'', ''contour'')'; end
            error('ltratio:args','%s: opts.method must be ''auto'', ''inversion'' (the single inversion of Bao and Kan) or ''engine'' (the engines'' computation)%s.',fn,x);
        end
        varargout = {m};
    case 'take'
        [fn,m,engine,cas,c,plain] = varargin{:};
        ic = (ischar(cas) || isstring(cas)) && any(strcmpi(char(cas),{'DK','DU','DT','SK','SU','ST'}));
        if strcmp(m,'inversion')
            if strcmp(engine,'general'), error('ltratio:args','%s: opts.method = ''inversion'' (the devices of Bao and Kan) does not apply with engine = ''general''.',fn); end
            if ~ic, error('ltratio:args','%s: opts.method = ''inversion'' applies to the cases DK, DU, DT, SK, SU, ST.',fn); end
            if ~plain, error('ltratio:args','%s: opts.method = ''inversion'' takes no options of the engines.',fn); end
            varargout = {true}; return
        end
        varargout = {strcmp(m,'auto') && strcmp(engine,'auto') && plain && ic && c <= 20};
    case 'plain'
        opts = varargin{1}; extra = {}; if numel(varargin) > 1, extra = varargin{2}; end
        varargout = {all(ismember(fieldnames(opts),[{'c0','g0','method','verified','tol','engine','start','gls'} extra]))};
    case 'tol'
        opts = varargin{1}; t = varargin{2};
        if isfield(opts,'tol') && ~isempty(opts.tol), t = opts.tol; end
        if ~(isnumeric(t) && isscalar(t) && isreal(t) && t > 0 && isfinite(t)), error('ltratio:args','ltratio: opts.tol must be a positive finite number.'); end
        varargout = {t};
    case 'warn'
        [fn,err,T] = varargin{:};
        bad = ~(err(:) <= T(:));
        if any(bad)
            warning('ltratio:inversion','%s: the error estimate of the single inversion exceeds its tolerance at %d of %d points (largest %g); opts.verified = true recomputes them by the engines'' checked computation.',fn,nnz(bad),numel(err),max(err(bad)));
        end
    case 'model'
        cas = upper(char(varargin{1}));
        varargout = {find('KUT' == cas(2)), cas(1) == 'S'};
    otherwise
        error('ltratio:args','kan_route: unknown operation %s.',op);
end
end
