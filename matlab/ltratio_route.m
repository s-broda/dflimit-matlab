function varargout = ltratio_route(op,varargin)
%LTRATIO_ROUTE  Helpers of the router between the classical programs and the general engine.
%
%   [engine,oc] = LTRATIO_ROUTE('engine',opts): opts.engine ('auto', the
%   default, 'classical' or 'general') and the options without the fields of
%   the router (engine, start, gls).
%   ok = LTRATIO_ROUTE('valid',S,c,oc): whether the classical programs accept
%   the case S.case of the spec S (LTRATIO_SPEC) at c with the options oc (the
%   checks of LTRATIO_ARGS).
%   ok = LTRATIO_ROUTE('plain',oc,allowed): no representation named, not the
%   checked computation, no option of the classical engines beyond the cell
%   allowed (c0, g0 and the method 'auto' always): the rule of LTRATIO_ENGINE
%   applies (otherwise the classical programs, as before).
%   LTRATIO_ROUTE('nomethod',opts): the error of a representation named with
%   engine = 'general'.
%   o = LTRATIO_ROUTE('gopts',S,opts): the options of the general engine for
%   the spec S (start, g0, gls, c0 and the rules).
%
%   See also LTRATIO_SPEC, LTRATIO_ENGINE.
switch op
    case 'engine'
        opts = varargin{1};
        if ~isstruct(opts), error('ltratio:args','ltratio: opts must be a structure.'); end
        engine = 'auto'; if isfield(opts,'engine') && ~isempty(opts.engine), engine = lower(char(opts.engine)); end
        if ~any(strcmp(engine,{'auto','classical','general'})), error('ltratio:args','ltratio: opts.engine must be ''auto'', ''classical'' or ''general''.'); end
        oc = opts;
        for f = {'engine','start','gls'}, if isfield(oc,f{1}), oc = rmfield(oc,f{1}); end, end
        varargout = {engine,oc};
    case 'valid'
        S = varargin{1}; c = varargin{2}; o = varargin{3}; o.g0 = S.g0;
        for f = {'engine','start','gls'}, if isfield(o,f{1}), o = rmfield(o,f{1}); end, end
        try
            ltratio_args(S.case,c,o); ok = true;
        catch ME
            if ~strcmp(ME.identifier,'ltratio:args'), rethrow(ME); end
            ok = false;
        end
        varargout = {ok};
    case 'plain'
        oc = varargin{1}; allowed = {}; if numel(varargin) > 1, allowed = varargin{2}; end
        fn = fieldnames(oc);
        ok = all(ismember(fn,[{'c0','g0','method','verified'} allowed]));
        if isfield(oc,'method') && ~isempty(oc.method), ok = ok && strcmpi(char(oc.method),'auto'); end
        if isfield(oc,'verified') && ~isempty(oc.verified), ok = ok && ~logical(oc.verified); end
        varargout = {ok};
    case 'nomethod'
        opts = varargin{1};
        if isfield(opts,'method') && ~isempty(opts.method) && ~strcmpi(char(opts.method),'auto')
            error('ltratio:args','ltratio: opts.method = ''%s'' names a representation of the classical programs; it does not apply with engine = ''general''.',char(opts.method));
        end
    case 'gopts'
        S = varargin{1}; opts = varargin{2};
        o = struct('start',S.start,'g0',S.g0,'gls',S.gls);
        for f = {'c0','verified','vtol','phi','hfac','n','ntail','alpha','filon','warn','tol'}
            if isfield(opts,f{1}), o.(f{1}) = opts.(f{1}); end
        end
        varargout = {o};
    otherwise
        error('ltratio:args','ltratio_route: unknown operation %s.',op);
end
end
