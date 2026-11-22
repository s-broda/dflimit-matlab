function S = ltratio_spec(cas,opts)
%LTRATIO_SPEC  The design, the start and the classical case of a computation (the router of the two engines).
%
%   S = LTRATIO_SPEC(cas,opts): cas a case string ('DK', 'DU', 'DT', 'SK',
%   'SU', 'ST', 'DG') or a design of the general engine (DFDET_TREND,
%   DFDET_PERRON, ...).  S.D the design, S.start ('zero', 'fixed',
%   'stationary'), S.g0, S.gls (GLS detrending at cbar, [] for least
%   squares) and S.case, the classical case of the design and start, or ''
%   where only the general engine has the computation: dfdet_trend(-1),
%   dfdet_trend(0), dfdet_trend(1) (the monomial basis, the one-step
%   regression) with a zero, fixed or stationary start, and dfdet_trend(1)
%   GLS-detrended at cbar = -13.5 with a zero start (DG).  For a design the
%   start is opts.start, or 'fixed' when opts.g0 ~= 0, and the detrending
%   opts.gls; for a case string the start and the detrending of the case.
%
%   See also LTRATIO_ENGINE, LTRATIO_CDF.
if nargin < 2 || isempty(opts), opts = struct; end
g0 = 0; if isfield(opts,'g0') && ~isempty(opts.g0), g0 = opts.g0; end
if ischar(cas) || isstring(cas)
    cas = upper(char(cas));
    if ~any(strcmp(cas,{'DK','DU','DT','SK','SU','ST','DG'})), error('ltratio:args','ltratio: the case must be one of DK, DU, DT, SK, SU, ST, DG.'); end
    if strcmp(cas,'DG'), S = struct('D',dfdet_trend(1),'start','zero','g0',0,'gls',-13.5,'case','DG'); return, end
    p = find('KUT' == cas(2)) - 2;
    if cas(1) == 'S', st = 'stationary'; elseif g0 ~= 0, st = 'fixed'; else, st = 'zero'; end
    S = struct('D',dfdet_trend(p),'start',st,'g0',g0,'gls',[],'case',cas);
    return
end
if ~(isstruct(cas) && isfield(cas,'br') && isfield(cas,'C')), error('ltratio:args','ltratio: the case must be one of DK, DU, DT, SK, SU, ST, DG, or a design of the general engine (dfdet_trend, dfdet_perron, ...).'); end
D = cas;
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('ltratio:args','ltratio: opts.g0 must be a finite real scalar.'); end
st = ''; if isfield(opts,'start') && ~isempty(opts.start), st = lower(char(opts.start)); end
if isempty(st), if g0 ~= 0, st = 'fixed'; else, st = 'zero'; end, end
if ~any(strcmp(st,{'zero','fixed','stationary'})), error('ltratio:args','ltratio: opts.start must be ''zero'', ''fixed'' or ''stationary''.'); end
if strcmp(st,'stationary') && g0 ~= 0, error('ltratio:args','ltratio: g0 (a fixed start) and start = ''stationary'' exclude each other.'); end
if strcmp(st,'fixed') && g0 == 0, st = 'zero'; end
gls = []; if isfield(opts,'gls') && ~isempty(opts.gls), gls = double(opts.gls); end
cs = '';
L = 'KUT';
for p = -1:1
    T = dfdet_trend(p);
    if isequal(size(D.C),[1 1]) && isequal(D.br,T.br) && isequal(size(D.C{1}),size(T.C{1})) && isequal(D.C{1},T.C{1}) && ~D.loc && strcmp(D.detrend,'onestep')
        if isempty(gls)
            if strcmp(st,'stationary'), cs = ['S' L(p + 2)]; else, cs = ['D' L(p + 2)]; end
        elseif p == 1 && strcmp(st,'zero') && gls == -13.5
            cs = 'DG';
        end
        break
    end
end
S = struct('D',D,'start',st,'g0',g0,'gls',gls,'case',cs);
end
