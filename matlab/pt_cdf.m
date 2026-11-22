function F = pt_cdf(q,D,cbar,c,opts)
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
%   whose characteristic function is a determinant (DFDET_QGRID); one
%   inversion of Gil-Pelaez on the grid of opts.N = 60000 steps up to
%   opts.umax = 4000.  The test rejects for small values.
%
%   Example: pt_cdf(6.789, dfdet_perron('B', 0.5), -17.8)   % 0.0455
%
%   See also PT_PDF, PT_QUANTILE, PT_POWER, PT_ENVELOPE.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
G = dfdet_qgrid(D,double(cbar),double(c),getf(opts,'N',60000),getf(opts,'umax',4000));
F = zeros(size(q));
for k = 1:numel(q), [F(k),~] = ptF(G,q(k)); end
end

function [F,f] = ptF(G,q)
z = exp(G.lM - 1i*G.u*q);
F = 0.5 - sum(imag(z)./G.u.*G.w)/pi; f = sum(real(z).*G.w)/pi;
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
