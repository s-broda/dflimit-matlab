function [F,f] = dfdet_ptf(G,q)
%DFDET_PTF  P(Q <= q) and the density of the limit of the point-optimal statistic at one point, on the grid of DFDET_QGRID.
%
%   [F,f] = DFDET_PTF(G,q): the inversion formula of Gil-Pelaez on the grid G.
%   |q| at most 0.15/h, h = asinh(umax)/N the step of the grid in v (1000
%   with the default grid): beyond, the steps of u near the origin do not
%   resolve e^{-iuq} (with a linear trend at c = 0, F came out as 1 - 8e-8 at
%   q = 1e4, 3.0 at 1e5 and 24 at 1e6), and an error is raised.  0, 1 and
%   NaN at -Inf, Inf and NaN.  Shared by PT_CDF, PT_PDF, PT_QUANTILE.
%
%   See also DFDET_QGRID, PT_CDF.
if isnan(q), F = NaN; f = NaN; return, end
if isinf(q), F = double(q > 0); f = 0; return, end
if ~(abs(q)*G.h <= 0.15)
    error('dfdet:args','the point-optimal statistic: q = %g lies beyond the range of the inversion (|q| at most %.3g = 0.15 N/asinh(umax); a larger opts.N extends it).',q,0.15/G.h);
end
z = exp(G.lM - 1i*G.u*q);
F = 0.5 - sum(imag(z)./G.u.*G.w)/pi; f = sum(real(z).*G.w)/pi;
end
