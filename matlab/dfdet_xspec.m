function S = dfdet_xspec(c0,q,A,opts)
%DFDET_XSPEC  One parabola term of an event of the general engine.
%
%   S = DFDET_XSPEC(c0,q,A,opts): the term w P(R <= X(sqrt V), sqrt V <= vmax)
%   with R = U + a (a the shift of Mt = e^{a s} M, DFDET_MODEL) and the parabola
%   X(v) = c0 v^2 + q v + A; opts.vmax (Inf), opts.w (1), and the derivatives
%   (dc0, dq, dA) of the coefficients in the parameter whose density is
%   wanted (0); opts.bnd (true): with a finite vmax the boundary term
%   -e^{-t vmax^2} P(R <= X(vmax) under the tilt).  The contour integral of
%   such a term is that of the t-ratio with the constant A in place of a and
%   the pieces of v ending at vmax (DFDET_CONTOUR, mode 'specs'): the t-ratio
%   and the coefficient statistic have A = a; the line U - aV <= b has
%   X(v) = a v^2 + a + b.  An event is a row of such terms.
%
%   See also DFDET_CONTOUR, DFDET_LIN_CDF, DFDET_ADFMAX_CDF.
if nargin < 4 || isempty(opts), opts = struct; end
S = struct('c0',double(c0),'q',double(q),'A',double(A),'vmax',getf(opts,'vmax',Inf),'w',getf(opts,'w',1), ...
    'dc0',getf(opts,'dc0',0),'dq',getf(opts,'dq',0),'dA',getf(opts,'dA',0),'bnd',logical(getf(opts,'bnd',true)));
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = double(o.(f)); else, v = d; end
end
