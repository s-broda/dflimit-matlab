function [cas,c,c0,g0,method] = ltratio_args(cas,c,opts,engines)
%LTRATIO_ARGS  Validation of the arguments of the LTRATIO functions.
%
%   [cas,c,c0,g0,method] = LTRATIO_ARGS(cas,c,opts) checks the case ('DK',
%   'DU', 'DT', 'SK', 'SU', 'ST', or 'DG', the GLS-detrended trend case of
%   Elliott, Rothenberg and Stock at cbar = -13.5, zero start, line integral
%   only; upper case returned), the local-to-unity
%   parameter c (a finite real scalar, at least -550; a stationary start
%   requires c <= -1e-4),
%   and the options c0 (the local null value, default 0), g0 (the fixed start
%   y_0 = g0 sigma sqrt(n), default 0; fixed-start cases only) and method
%   ('auto' by default: the faster representation point by point, see
%   LTRATIO_METHOD; 'single', the line integral representation, also
%   accepted as 'line'; 'double', the conditional integral representation,
%   also accepted as 'conditional'; 'contour', the contour integral representation of
%   TR_CONTOUR, for a fixed start g0 ~= 0 in cases DU, DT only, and an error
%   elsewhere), returned as 'auto', 'single', 'double' or 'contour'.  These
%   programs have no line integral with a fixed start: there 'single' is
%   treated as 'auto' in cases DU, DT (the contour integral) and a fixed start
%   forces 'double' in case DK.  A fixed start is checked against the region
%   of the C and Julia engines (case DK c <= 2, |g0| <= 2, |g0| <= 1 for
%   c > 0; DU -100 <= c <= 2, |g0| <= 1.5, (c g0)^2/2 <= 4050; DT
%   -100 <= c <= 2, |g0| <= 3, (c g0)^2/2 <= 7200), and where the conditional
%   integral is the method (case DK, and opts.method = 'double' in cases DU,
%   DT) also against c >= -20, or c >= -40 with (c g0)^2/2 <= 100.
%
%   See also LTRATIO_CDF, LTRATIO_PDF, LTRATIO_ICDF, LTRATIO_MOM.
if ~(ischar(cas) || isstring(cas)), error('ltratio:args','ltratio: the case must be one of DK, DU, DT, SK, SU, ST, DG.'); end
cas = upper(char(cas));
if ~any(strcmp(cas,{'DK','DU','DT','SK','SU','ST','DG'}))
    error('ltratio:args','ltratio: the case must be one of DK, DU, DT, SK, SU, ST, DG.');
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','ltratio: c must be a finite real scalar.'); end
% the range of c with a stationary start (the same check and message in the C and Julia engines: DFL_CSTAT, CSTAT, where
% the measurements are): the variance -1/(2c) of the start grows without bound as c -> 0-, and in case SK the conditional
% integral flags points off by up to 2e-8 from |c| ~ 1e-4 and is off by 1e-3 without a flag at c = -1e-6
if cas(1) == 'S' && ~(c <= -1e-4), error('ltratio:args','ltratio: a stationary start (cases SK, SU, ST) requires c <= -1e-4: closer to zero the variance -1/(2c) of the start is too large for the quadrature.'); end
% the range of c (the same check and message in the C and Julia engines: DFL_CMIN, CMIN)
if c < -550, error('ltratio:args','ltratio: c must be at least -550 (the range in which the programs are checked; the moments fail from about c = -600, the line integral from about c = -650).'); end
if ~isstruct(opts), error('ltratio:args','ltratio: opts must be a structure.'); end
c0 = 0; if isfield(opts,'c0') && ~isempty(opts.c0), c0 = opts.c0; end
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('ltratio:args','ltratio: opts.c0 must be a finite real scalar.'); end
g0 = 0; if isfield(opts,'g0') && ~isempty(opts.g0), g0 = opts.g0; end
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('ltratio:args','ltratio: opts.g0 must be a finite real scalar.'); end
if g0 ~= 0 && cas(1) == 'S', error('ltratio:args','ltratio: opts.g0 (a fixed start) applies to cases DK, DU, DT only.'); end
method = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), method = lower(char(opts.method)); end
% the companion paper's names, the line integral and the conditional integral, are synonyms
if strcmp(method,'line'), method = 'single'; elseif strcmp(method,'conditional'), method = 'double'; end
if ~any(strcmp(method,{'auto','single','double','contour'})), error('ltratio:args','ltratio: opts.method must be ''auto'', ''single'' (''line''), ''double'' (''conditional'') or ''contour''.'); end
% the contour integral (TR_CONTOUR): a fixed start in cases DU, DT only (the same check and message in the C and Julia
% engines)
dut = any(strcmp(cas,{'DU','DT'})) && g0 ~= 0;
if strcmp(method,'contour') && ~dut, error('ltratio:args','ltratio: method ''contour'' (the contour integral) applies to a fixed start g0 != 0 in cases DU, DT.'); end
% with a fixed start these programs have no line integral (it needs a central form of U + 1/2 here): 'single' is the
% automatic choice in cases DU, DT (the contour integral, LTRATIO_METHOD), and the conditional integral in case DK
if g0 ~= 0 && strcmp(method,'single'), if dut, method = 'auto'; else, method = 'double'; end, end
if g0 ~= 0 && ~dut, method = 'double'; end
% the region of a fixed start (the C and Julia engines: dfl_fixed_ok, fixed_ok), and where the conditional integral is
% the method (case DK, and an explicit 'double' in cases DU, DT) also c >= -20, or c >= -40 with (c g0)^2/2 <= 100
% (which covers a start of up to three stationary standard deviations, g0 = xi/sqrt(-2c), xi <= 3, as in the companion
% paper's table of power with the initial observation: within 5e-13 of the line integral of the Julia engine there); beyond, it flags most values, off by up to
% 3e-2 at c = -200 in case DK, and at c = -60 it is off by 1e-10 without a flag.  The contour integral covers the
% region of the C and Julia engines in cases DU, DT.
if g0 ~= 0 && any(strcmp(cas,{'DK','DU','DT'}))
    l = (c*g0)^2/2;
    switch cas
        case 'DK', ok = c <= 2 && abs(g0) <= 2 - (c > 0);
        case 'DU', ok = c >= -100 && c <= 2 && abs(g0) <= 1.5 && l <= 4050;
        otherwise, ok = c >= -100 && c <= 2 && abs(g0) <= 3 && l <= 7200;
    end
    if ~ok, error('ltratio:args','ltratio: with a fixed start g0 != 0 the programs are checked in case DK for c <= 2, |g0| <= 2 (|g0| <= 1 for c > 0), in case DU for -100 <= c <= 2, |g0| <= 1.5 and (c g0)^2/2 <= 4050, in case DT for -100 <= c <= 2, |g0| <= 3 and (c g0)^2/2 <= 7200.'); end
    if nargin < 4, engines = true; end
    % (engines false: the checks of the classical engines of DFLimit.jl and the C engine alone, for the single inversion
    % of Bao and Kan, which does not take the conditional integral)
    if engines && (~dut || strcmp(method,'double')) && c < -20 && ~(c >= -40 && l <= 100), error('ltratio:args','ltratio: with a fixed start g0 ~= 0 the conditional integral of the MATLAB programs (case DK, and opts.method = ''double'' in cases DU, DT) is checked for c >= -20, and for c >= -40 where (c g0)^2/2 <= 100; in case DK the MEX file (c/mex) and the Julia engine, which have the line integral, cover the other fixed-start bounds without it (all c <= 2), and in cases DU, DT the contour integral (opts.method ''auto'' or ''contour'') covers them down to c = -100.'); end
end
if strcmp(cas,'DG')                                  % the GLS-detrended trend case: zero start, the line integral only
    if g0 ~= 0, error('ltratio:args','ltratio: opts.g0 (a fixed start) applies to cases DK, DU, DT only.'); end
    if strcmp(method,'double'), error('ltratio:args','ltratio: case DG has the line integral representation only.'); end
    method = 'single';
end
end
