function [F,info] = tr_contour(x,cas,c,opts,mode)
%TR_CONTOUR  Limiting distribution and density of the t-ratio with a fixed start by the contour integral (C engine, MEX file LTRATIO_MEX).
%
%   [F,info] = TR_CONTOUR(x,cas,c,opts,mode): a drop-in replacement for the MATLAB engine of the same name, with the
%   same arguments: mode 'cdf' (default) or 'pdf'; cases 'DU' and 'DT' only, with a fixed start opts.g0 ~= 0
%   (required); opts.c0 the local null value (default 0, as in the MATLAB engine; LTRATIO_CDF passes its own c0).  The
%   contour integral representation (dfl_contour in c/dfl_gen.c, a port of tr_contour of julia/DFExact.jl): the
%   inversion formula on paths around the singularities of the moment generating function, with the integral over v
%   in closed form, for each t of an outer hyperbola two U-shaped paths in s around the roots
%   of D(., t), a double integral with an explicit integrand, computed twice on different paths with different rules (the
%   error estimate is their difference).  info: time; err, the error estimate per point; status, 0
%   where the estimate is at most 1e-10 and the value is finite, otherwise 2 (the warning tr_contour:status); nodes, the
%   number of outer nodes (of the two evaluations together); level, the level of the value (1 to 4).  Non-finite x: F = 0 at -Inf, 1 (cdf)
%   or 0 (pdf) at +Inf, NaN at NaN.
%
%   See also LTRATIO_MEX, TR_LIMDIST, LTRATIO_CDF, LTRATIO_PDF.
if nargin < 4 || isempty(opts), opts = struct; end
if nargin < 5 || isempty(mode), mode = 'cdf'; end
msg = 'tr_contour: the contour integral applies to a fixed start g0 ~= 0 in cases DU, DT.';
if ~(ischar(cas) || isstring(cas)), error('tr_contour:args',msg); end
switch upper(char(cas))
    case 'DU', model = 2;
    case 'DT', model = 3;
    otherwise, error('tr_contour:args',msg);
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('tr_contour:args','tr_contour: c must be a finite real scalar.'); end
if ~(isnumeric(x) && isreal(x)), error('tr_contour:args','tr_contour: x must be real.'); end
if ~isstruct(opts), error('tr_contour:args','tr_contour: opts must be a structure.'); end
g0 = getf(opts,'g0',0); c0 = getf(opts,'c0',0);
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('tr_contour:args','tr_contour: opts.g0 must be a finite real scalar.'); end
if g0 == 0, error('tr_contour:args',msg); end
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('tr_contour:args','tr_contour: opts.c0 must be a finite real scalar.'); end
switch mode
    case 'cdf', kind = 2;
    case 'pdf', kind = 1;
    otherwise, error('tr_contour:args','tr_contour: mode must be ''cdf'' or ''pdf''.');
end
t0 = tic;
verified = isfield(opts,'verified') && ~isempty(opts.verified) && logical(opts.verified);
[F,err,status,nodes,level] = ltratio_mex('contour', kind, double(x(:)), model, double(c), double(c0), double(g0), double(verified));
sz = size(x); F = reshape(F,sz);
info = struct('time',toc(t0),'err',reshape(err,sz),'status',reshape(status,sz),'nodes',reshape(nodes,sz),'level',reshape(level,sz));
if verified && any(status(:))
    iw = find(status,1);
    warning('tr_contour:status',['tr_contour: %d of %d points with nonzero status (first: x = %.6g, error estimate ' ...
        '%.3g); the result may be inaccurate there.'],nnz(status),numel(x),x(iw),err(iw));
end
end
function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
