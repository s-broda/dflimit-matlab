function F = ltratio_rangecheck(fn,args,xname,x,F)
%LTRATIO_RANGECHECK  The final check of a distribution function.
%
%   F = LTRATIO_RANGECHECK(fn,args,xname,x,F) returns F, or raises the error
%   ltratio:range that names the call fn, its arguments args, the point
%   (xname = x(k)) and the value where a value of F lies below -1e-8 or above
%   1 + 1e-8.  Where the representations hold their errors are far below 1e-8,
%   so such a value is a computation that failed, never rounding; the check
%   costs nothing and recomputes nothing (cdf_range_check of DFExact.jl and
%   DFDet.jl, dfl_cdf_range_check of the C engine).  LTRATIO_CDF, COEF_CDF,
%   DFDET_CDF and DFDET_COEF_CDF apply it to the values they return.
%
%   See also LTRATIO_CDF, COEF_CDF, DFDET_CDF, DFDET_COEF_CDF.
tol = 1e-8;
bad = find(F < -tol | F > 1 + tol, 1);
if isempty(bad), return, end
error('ltratio:range',['%s: %s, %s = %g: the computed distribution function is %.12g, outside [0, 1] by more than %g; ' ...
      'the computation fails here and no value is returned'],fn,args,xname,x(bad),F(bad),tol);
end
