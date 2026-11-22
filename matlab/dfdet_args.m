function [c,c0,start,g0,gls,o] = dfdet_args(D,c,opts)
%DFDET_ARGS  Validation of the arguments of the user-level routines of the general engine.
%
%   [c,c0,start,g0,gls,o] = DFDET_ARGS(D,c,opts) checks the design D (a
%   structure from DFDET_DETERM and its constructors), the local-to-unity
%   parameter c (a finite real scalar) and the options: c0 (the local null
%   value, 0), start ('zero', 'fixed' or 'stationary'; by default 'fixed' when
%   opts.g0 ~= 0, else 'zero'), g0 (the fixed start J(0) = g0, y_0 = g0 sigma
%   sqrt(n); 0), gls (GLS detrending at cbar = gls < 0, [] for least squares);
%   a stationary start needs c < 0 and excludes g0; GLS detrending and the
%   statistic that detrends first need the zero start.  o: the options
%   passed on to the engine (verified, vtol, phi, hfac, n, ntail, alpha,
%   filon, warn, tol).
%
%   See also DFDET_CDF.
if ~(isstruct(D) && isfield(D,'br') && isfield(D,'C') && isfield(D,'detrend')), error('dfdet:args','dfdet: the design must be a structure built by dfdet_determ, dfdet_trend, dfdet_perron, ...'); end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','dfdet: c must be a finite real scalar.'); end
c = double(c);
if ~isstruct(opts), error('dfdet:args','dfdet: opts must be a structure.'); end
c0 = getf(opts,'c0',0);
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('dfdet:args','dfdet: opts.c0 must be a finite real scalar.'); end
g0 = getf(opts,'g0',0);
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('dfdet:args','dfdet: opts.g0 must be a finite real scalar.'); end
start = getf(opts,'start','');
if isempty(start), if g0 ~= 0, start = 'fixed'; else, start = 'zero'; end, end
start = lower(char(start));
if ~any(strcmp(start,{'zero','fixed','stationary'})), error('dfdet:args','dfdet: opts.start must be ''zero'', ''fixed'' or ''stationary''.'); end
if strcmp(start,'stationary') && g0 ~= 0, error('dfdet:args','dfdet: opts.g0 (a fixed start) and start = ''stationary'' exclude each other.'); end
if strcmp(start,'fixed') && g0 == 0, start = 'zero'; end
if strcmp(start,'stationary') && ~(c < 0), error('dfdet:args','dfdet: %s.',dfdet_msg('stat_c')); end
gls = getf(opts,'gls',[]);
if ~isempty(gls)
    if ~(isnumeric(gls) && isscalar(gls) && isreal(gls) && gls < 0), error('dfdet:args','dfdet: GLS detrending needs opts.gls = cbar < 0.'); end
    if ~strcmp(start,'zero'), error('dfdet:args','dfdet: %s.',dfdet_msg('gls_start')); end
    if strcmp(D.detrend,'first'), error('dfdet:args','dfdet: %s.',dfdet_msg('first_gls')); end
    gls = double(gls);
end
if strcmp(D.detrend,'first') && ~strcmp(start,'zero'), error('dfdet:args','dfdet: %s.',dfdet_msg('first_start')); end
o = struct;
for f = {'verified','vtol','phi','hfac','n','ntail','alpha','filon','warn','tol'}
    if isfield(opts,f{1}) && ~isempty(opts.(f{1})), o.(f{1}) = opts.(f{1}); end
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
