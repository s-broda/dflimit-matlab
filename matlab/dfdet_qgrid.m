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
%   PT_QUANTILE, PT_POWER.  The transform of Q has a singularity at the
%   distance d from the real line (the first zeros of its determinant on the
%   real line, QEDGE below), which shrinks like e^{-2c} above the unit root,
%   and the rule errs by about exp(-2 pi d/h): below d/h = 3.3 an error is
%   raised that names the N that would serve (with a linear trend at cbar =
%   -13.5 from about c = 3.7; the distribution function came out between
%   -0.1 and 3.4 at c = 5).  G.h is the step h.
%
%   See also PT_CDF.
if nargin < 4 || isempty(N), N = 60000; end
if nargin < 5 || isempty(umax), umax = 4000; end
key = sprintf('%s|pt%.17g|%.17g|%d|%.17g',dfdet_key(D),cbar,c,N,umax);
G = dfdet_cache('get','qgrid',key);
if ~isempty(G), return, end
P = dfdet_prep(D,'pt',cbar); Mo = dfdet_model(P);
h = asinh(umax)/N;
d = qedge(Mo,cbar,c);
if ~(d/h >= 3.3)
    error('dfdet:args','the point-optimal statistic at c = %g, cbar = %g: the transform of Q has a singularity at %.3g from the real line, which the grid of the inversion (step %.3g near the origin) does not resolve (far above the unit root); it needs N >= %d (opts.N).',c,cbar,d,h,ceil(3.3*asinh(umax)/d));
end
v = (0:N)'*h; u = sinh(v); u(1) = 1e-7; w = cosh(v)*h; w(1) = w(1)/2;
L0 = complex(zeros(N + 1,1)); LD = L0;
for k = 1:N+1
    [L0(k),LD(k)] = dfdet_logmgf(Mo,c,complex(0,u(k)),complex(0,cbar^2*u(k)));
end
ph = imag(LD(1)) + [0; cumsum(angle(exp(1i*diff(imag(LD)))))];
lM = L0 - complex(real(LD),ph)/2 + 1i*u/2;
G = struct('u',u,'w',w,'lM',lM,'h',h);
dfdet_cache('put','qgrid',key,G);
end

function e = qedge(Mo,cbar,c)
% the distance of the nearest singularity of phi(u) = E e^{iuQ} from the real line: Q is a quadratic form, whose transform
% M(s, cbar^2 s) is singular on the real line only, at the first zeros s_+ > 0 and s_- < 0 of its determinant (scanned by
% factors of 1.1 from 1e-10, then bisected), the smaller of s_+ and |s_-|
e = Inf;
for sg = [1 -1]
    s = 1e-10;
    while s < 1e4 && ~bad(Mo,c,cbar,sg*1.1*s), s = s*1.1; end
    if ~(s < 1e4), continue, end
    lo = s; hi = 1.1*s;
    for it = 1:40, m = sqrt(lo*hi); if bad(Mo,c,cbar,sg*m), hi = m; else, lo = m; end, end
    e = min(e,lo);
end
end

function b = bad(Mo,c,cbar,s)
[~,ld] = dfdet_logmgf(Mo,c,complex(s,0),complex(cbar^2*s,0));
b = ~(isfinite(real(ld)) && abs(imag(ld)) < 1);
end
