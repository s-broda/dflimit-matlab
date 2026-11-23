function [F,err] = pt_cdf(q,D,cbar,c,opts)
%PT_CDF  Limiting distribution function of the point-optimal statistic of Elliott, Rothenberg and Stock, any design.
%
%   F = PT_CDF(q,D,cbar,c,opts) is the limit of P(P_T <= q) at the points q (an
%   array), P_T the point-optimal statistic against cbar < 0,
%   {S(abar) - abar S(1)}/s^2 (S(a) the sum of squared residuals of the
%   regression of the data quasi-differenced at a = 1 + cbar/n on the trend
%   functions transformed in the same way) for the design D, under
%   phi = 1 + c/n (c defaults to 0), zero start: with delta the trend
%   functions that vanish at the origin and do not jump, psi = delta' -
%   cbar delta, Xi = int psi psi', eta = int psi (dJ - cbar J dr), and psi_0,
%   Xi_0, eta_0 the same at cbar = 0, the limit is Q = cbar^2 int J^2 -
%   cbar J(1)^2 - eta' Xi^{-1} eta + eta_0' Xi_0^{-1} eta_0, a quadratic form
%   whose characteristic function is a determinant: one inversion of its
%   characteristic function.  The test rejects for small values.  By default
%   (opts.method 'auto', up to c = 4) the single inversion of Yong Bao and
%   Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least
%   Squares Estimator in First-Order Autoregressions: O(1) Algorithms":
%   adaptive quadrature over decades with the branch of the square root of
%   the determinant from the pivots of the precision matrix at each point and
%   the tail completed by integration by parts (DFDET_KAN), for any q;
%   opts.tol (1e-12) its absolute tolerance, err the error estimates, and the
%   warning ltratio:inversion names the points at which they exceed it.
%   opts.method 'inversion' forces it; 'grid' (or opts.N, opts.umax given, and
%   above c = 4) the inversion of Gil-Pelaez on the fixed grid of opts.N =
%   60000 steps up to opts.umax = 4000 (DFDET_QGRID), the computation of
%   before: an error is raised where the grid does not resolve the transform
%   (far above the unit root) and at |q| > 0.15/h, h = asinh(umax)/N (1000
%   with the default grid), where the steps of u near the origin do not
%   resolve e^{-iuq} (F came out as 3.0 at q = 1e5 at c = 0), and err is NaN.
%   0, 1 and NaN at -Inf, Inf and NaN.
%
%   Example: pt_cdf(6.789, dfdet_perron('B', 0.5), -17.8)   % 0.0455
%
%   See also PT_PDF, PT_QUANTILE, PT_POWER, PT_ENVELOPE, DFDET_PTF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
[grid,tol] = dfdet_ptmethod(opts,c,'pt_cdf');
F = zeros(size(q)); err = NaN(size(q));
if grid
    G = dfdet_qgrid(D,double(cbar),double(c),getf(opts,'N',60000),getf(opts,'umax',4000));
    for k = 1:numel(q), F(k) = dfdet_ptf(G,q(k)); end
    return
end
K = pt_transform(D,cbar,c);
for k = 1:numel(q), [F(k),err(k)] = dfdet_kan('pt_cdf1',K,double(cbar),q(k),tol); end
pt_warn('pt_cdf',err,tol);
end

function K = pt_transform(D,cbar,c)
% the transform of the point-optimal statistic at c for the single inversion (the checks of DFDET_QGRID)
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','the point-optimal statistic: c must be a finite real scalar.'); end
K = dfdet_kan('cf',dfdet_model(dfdet_prep(D,'pt',double(cbar))),double(c));
end

function pt_warn(fn,err,tol)
bad = ~(err(:) <= tol);
if any(bad), warning('ltratio:inversion','%s: the error estimate of the single inversion exceeds the tolerance %g at %d of %d points (largest %g); opts.method = ''grid'' gives the inversion on the fixed grid.',fn,tol,nnz(bad),numel(err),max(err(bad))); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
