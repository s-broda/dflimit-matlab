function dfdet_method(opts,fn)
%DFDET_METHOD  The option method of DFDET_CDF, DFDET_PDF and DFDET_ICDF.
%
%   DFDET_METHOD(opts,fn) checks opts.method: 'auto' (the default) and
%   'paths' take the inversion on paths around the roots, the only
%   representation of the general engine in these programs.  'line', the line
%   integral with the explicit kernel for the designs whose numerator has
%   rank one or two (polynomial trends, an intercept with sines and cosines,
%   the statistic that detrends first; the C engine, the Julia, R and Python
%   packages, where opts.method = 'auto' takes it at c0 = 0), is not
%   available in the MATLAB programs: their line integral (DFX_CORE, without
%   the fixed node layouts of the C and Julia engines) is 2.8 to 9 times
%   slower than their inversion on paths in the rank-two cases DU and DT
%   at c0 = 0 even with its closed forms (validation/engine_table_matlab.txt,
%   c = 0, -10, -50; November 27, 2026),
%   so that it would not be faster for any design.  With the folder of the
%   MEX file (c/mex) ahead on the path, LTRATIO_CDF, LTRATIO_PDF and
%   LTRATIO_ICDF with a design and opts.method = 'line' take the C engine's
%   line integral.
%
%   See also DFDET_CDF, LTRATIO_CROUTE.
m = 'auto';
if isfield(opts,'method') && ~isempty(opts.method)
    if ~(ischar(opts.method) || isstring(opts.method)), error('dfdet:args','%s: opts.method must be ''auto'', ''line'' or ''paths''.',fn); end
    m = lower(char(opts.method));
end
switch m
    case {'auto','paths'}
    case 'line'
        error('dfdet:line',['%s: opts.method = ''line'', the line integral with the explicit kernel, is not available in the MATLAB programs ' ...
            '(their line integral is slower than the inversion on paths for these designs); use the MEX file: with the folder c/mex ahead on the path, ' ...
            'ltratio_cdf(x,D,c,struct(''method'',''line'')) or ltratio(''cdf'',x,D,c,struct(''method'',''line'')).'],fn);
    otherwise
        % (also 'lines' and 'mixture' without a normal component, rho = 1: the refusal and the message of DFLimit.jl)
        error('dfdet:args','%s: opts.method must be ''auto'', ''line'' (the line integral with the explicit kernel, in the MEX file) or ''paths'' (the inversion on paths around the roots).',fn);
end
end
