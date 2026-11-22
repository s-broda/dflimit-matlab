function G = dfdet_qgrid(D,cbar,c,N,umax)
%DFDET_QGRID  The characteristic function of the limit of the point-optimal statistic on the grid of Gil-Pelaez.
%
%   G = DFDET_QGRID(D,cbar,c,N,umax): Q = cbar^2 V_g + ell' A ell (DFDET_PREP,
%   kind 'pt'), so that E exp(iuQ) = e^{iu/2} M(iu, iu cbar^2) with the
%   numerator of the point-optimal statistic in place of U + 1/2, the square
%   root of the determinant continued along u from the origin; the grid
%   u = sinh(v), v = k h, k = 0..N (60000), h = asinh(umax)/N (umax 4000),
%   u_0 = 1e-7, trapezoidal weights in v.  G.u, G.w, G.lM (the logarithm of
%   the characteristic function); cached.  Shared by PT_CDF, PT_PDF,
%   PT_QUANTILE, PT_POWER.
%
%   See also PT_CDF.
if nargin < 4 || isempty(N), N = 60000; end
if nargin < 5 || isempty(umax), umax = 4000; end
key = sprintf('%s|pt%.17g|%.17g|%d|%.17g',dfdet_key(D),cbar,c,N,umax);
G = dfdet_cache('get','qgrid',key);
if ~isempty(G), return, end
P = dfdet_prep(D,'pt',cbar); Mo = dfdet_model(P);
h = asinh(umax)/N; v = (0:N)'*h; u = sinh(v); u(1) = 1e-7; w = cosh(v)*h; w(1) = w(1)/2;
L0 = complex(zeros(N + 1,1)); LD = L0;
for k = 1:N+1
    [L0(k),LD(k)] = dfdet_logmgf(Mo,c,complex(0,u(k)),complex(0,cbar^2*u(k)));
end
ph = imag(LD(1)) + [0; cumsum(angle(exp(1i*diff(imag(LD)))))];
lM = L0 - complex(real(LD),ph)/2 + 1i*u/2;
G = struct('u',u,'w',w,'lM',lM);
dfdet_cache('put','qgrid',key,G);
end
