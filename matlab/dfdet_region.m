function P = dfdet_region(bnd,D,c,opts)
%DFDET_REGION  The limiting probability of a region {U <= g(V)} of the plane by the general engine.
%
%   P = DFDET_REGION(bnd,D,c,opts) for the boundaries that the contour integral
%   delivers: bnd.type 'line' (fields a, b: g(v) = b + a v, DFDET_LIN_CDF) and
%   'parabola' (fields x, c0: g(v) = x sqrt(v) + c0 v, the t-ratio, DFDET_CDF).
%   The weighted-average-power boundary ('wap') and a general function ('fun')
%   are refused: the kernels of the contour integral are in closed form for
%   boundaries quadratic in sqrt(V) only (LTRATIO_REGION has them in case DK).
%   Options start, g0, gls as in DFDET_CDF.
%
%   Example: dfdet_region(struct('type','line','a',-3.5,'b',-0.26), dfdet_trend(-1), -7)
%
%   See also LTRATIO_REGION, DFDET_LIN_CDF, DFDET_CDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
if ~isstruct(bnd) || ~isfield(bnd,'type'), error('dfdet:args','dfdet_region: bnd must be a structure with the field type.'); end
switch lower(char(bnd.type))
    case 'line'
        if ~(isnumeric(bnd.a) && isscalar(bnd.a) && isreal(bnd.a) && isfinite(bnd.a)), error('dfdet:args','dfdet_region: the slope a of the line must be a finite real scalar.'); end
        P = dfdet_lin_cdf(bnd.b,bnd.a,D,c,opts);
    case 'parabola'
        o = opts; o.c0 = 0; if isfield(bnd,'c0') && ~isempty(bnd.c0), o.c0 = bnd.c0; end
        P = dfdet_cdf(bnd.x,D,c,o);
    case {'wap','fun'}
        error('dfdet:args','dfdet_region: the boundary type %s is not available by the contour integral (its kernels are in closed form for boundaries quadratic in sqrt(V) only).',lower(char(bnd.type)));
    otherwise
        error('dfdet:args','dfdet_region: unknown boundary type %s.',char(bnd.type));
end
end
