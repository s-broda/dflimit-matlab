function [f,err] = pt_pdf(q,D,cbar,c,opts)
%PT_PDF  Limiting density of the point-optimal statistic against cbar for the design D (arguments as in PT_CDF).
%
%   [f,err] = PT_PDF(q,D,cbar,c,opts): by default (up to c = 4) the single
%   inversion of Bao and Kan, the density of a quadratic form by one integral
%   of the real part of its transform (DFDET_KAN); opts.method 'grid' (or
%   opts.N, opts.umax) the fixed grid of before (err NaN).
%
%   Example: pt_pdf(6.789, dfdet_perron('B', 0.5), -17.8)
%
%   See also PT_CDF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
[grid,tol] = dfdet_ptmethod(opts,c,'pt_pdf');
f = zeros(size(q)); err = NaN(size(q));
if grid
    G = dfdet_qgrid(D,double(cbar),double(c),getf(opts,'N',60000),getf(opts,'umax',4000));
    for k = 1:numel(q), [~,f(k)] = dfdet_ptf(G,q(k)); end   % (the range of q and the non-finite points as PT_CDF)
    return
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','the point-optimal statistic: c must be a finite real scalar.'); end
K = dfdet_kan('cf',dfdet_model(dfdet_prep(D,'pt',double(cbar))),double(c));
for k = 1:numel(q), [f(k),err(k)] = dfdet_kan('pt_pdf1',K,double(cbar),q(k),tol); end
bad = ~(err(:) <= tol);
if any(bad), warning('ltratio:inversion','pt_pdf: the error estimate of the single inversion exceeds the tolerance %g at %d of %d points (largest %g); opts.method = ''grid'' gives the inversion on the fixed grid.',tol,nnz(bad),numel(err),max(err(bad))); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
