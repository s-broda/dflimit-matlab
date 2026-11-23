function [grid,tol] = dfdet_ptmethod(opts,c,fn)
%DFDET_PTMETHOD  The computation of the point-optimal functions: the single inversion of Bao and Kan or the fixed grid.
%
%   [grid,tol] = DFDET_PTMETHOD(opts,c,fn): opts.method 'auto' (the default) and 'inversion' the single inversion of
%   Yong Bao and Raymond Kan ("Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in
%   First-Order Autoregressions: O(1) Algorithms") with the branch of the transform from the pivots of the precision
%   matrix and the tail completed by integration by parts (DFDET_KAN), up to c = 4 with 'auto'; 'grid', or opts.N or
%   opts.umax given, the inversion of Gil-Pelaez on the fixed grid of N steps up to umax (DFDET_QGRID), the computation of
%   before, with its refusals; and the tolerance opts.tol of the single inversion (1e-12).  (PT_INV_CMAX of DFLimit.jl:
%   compared in November 2026 with the contour integral on the statistic's quadratic form for an intercept, a linear and
%   a quadratic trend and Perron's model B, c from -30 to 7: within 3.3e-12, where the grid is off by up to 6e-8 or
%   refuses above the unit root; beyond c = 4 one value can take seconds, and the grid stays the default there.)
%
%   See also PT_CDF, DFDET_KAN.
m = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), m = lower(char(opts.method)); end
if ~any(strcmp(m,{'auto','inversion','grid'})), error('dfdet:args','%s: opts.method must be ''auto'', ''inversion'' (the single inversion of Bao and Kan) or ''grid'' (the fixed grid of Gil-Pelaez).',fn); end
given = (isfield(opts,'N') && ~isempty(opts.N)) || (isfield(opts,'umax') && ~isempty(opts.umax));
if strcmp(m,'inversion') && given, error('dfdet:args','%s: opts.N and opts.umax belong to opts.method = ''grid''.',fn); end
grid = strcmp(m,'grid') || (strcmp(m,'auto') && (given || ~(c <= 4)));
tol = 1e-12; if isfield(opts,'tol') && ~isempty(opts.tol), tol = opts.tol; end
if ~(isnumeric(tol) && isscalar(tol) && tol > 0 && isfinite(tol)), error('dfdet:args','%s: opts.tol must be a positive finite number.',fn); end
end
