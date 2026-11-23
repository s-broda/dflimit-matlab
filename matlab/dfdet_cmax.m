function b = dfdet_cmax(D,start,g0)
%DFDET_CMAX  The largest c at which the general engine computes.
%
%   b = DFDET_CMAX(D,start,g0): 20 for a design without deterministic terms
%   (no regressors), with any start; with deterministic terms 16, and 10 with
%   a fixed start (start 'fixed' and g0 ~= 0) (GEN_CMAX0, GEN_CMAX,
%   GEN_CMAX_FIXED of DFDet.jl, DFL_DET_CMAX0, DFL_DET_CMAX, DFL_DET_CMAX_FIXED
%   of the C engine).  Above the unit root the
%   transform has one root of the order 2c e^{-2c}, and the reciprocal roots,
%   the eigenvalues of a small matrix, span e^{2c}: with deterministic terms
%   the others are lost in double precision from about c = 17 (the paths of
%   case DU fail there, of case DT at c = 18), and with a fixed start the
%   coefficient statistic was off by 3e-5 in case DU at c = 12.  The range in which the engine
%   is checked above the unit root against the classical engine and the
%   single inversion of DFKan (DFLimit.jl, November 2026).  Above it every
%   routine of the general engine refuses (DFDET_ARGS), and LTRATIO_ENGINE
%   takes the classical programs.
%
%   See also DFDET_ARGS, LTRATIO_ENGINE.
if nargin < 2 || isempty(start), start = 'zero'; end
if nargin < 3 || isempty(g0), g0 = 0; end
if size(D.C{1},2) == 0, b = 20; elseif strcmp(start,'fixed') && g0 ~= 0, b = 10; else, b = 16; end
end
