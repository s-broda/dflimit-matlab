function [P,info] = tr_region(bnd,c,opts)
%TR_REGION  P{U <= g(V)}, case DK, zero start, by the C engine (MEX file LTRATIO_MEX): a drop-in replacement for the
%   MATLAB engine of the same name for the boundary types 'line' (fields a, b), 'parabola' (x, c0) and 'wap' (k, omega)
%   of LTRATIO_REGION, with bnd.v0 (P(V < v0) included); the type 'fun' (a function handle) and bnd.brks need the
%   MATLAB engine (matlab/tr_region.m: remove this folder from the path).
%
%   See also LTRATIO_MEX, LTRATIO_REGION.
if nargin < 3 || isempty(opts), opts = struct; end
v0 = 0; if isfield(bnd,'v0') && ~isempty(bnd.v0), v0 = bnd.v0; end
if isfield(bnd,'brks') && ~isempty(bnd.brks), error('tr_region:args','tr_region: the C engine takes no breakpoints (bnd.brks); use the MATLAB engine.'); end
switch lower(bnd.type)
    case 'line', type = 1; par = [bnd.a bnd.b];
    case 'parabola', type = 2; c0 = 0; if isfield(bnd,'c0') && ~isempty(bnd.c0), c0 = bnd.c0; end, par = [bnd.x c0];
    case 'wap', type = 3; par = [bnd.k bnd.omega];
    otherwise, error('tr_region:args','tr_region: the C engine takes the boundary types line, parabola and wap; for %s use the MATLAB engine (matlab/).',bnd.type);
end
t0 = tic;
nm = {'M','Nv','Nl','maxlev','plen','otol','stepw','errtol','blo','bhi','bmin','mumin','tdec'}; dopt = NaN(1,numel(nm));
for k = 1:numel(nm), if isfield(opts,nm{k}) && ~isempty(opts.(nm{k})), dopt(k) = double(opts.(nm{k})); end, end
[P,err,status] = ltratio_mex('region', type, double(par), double(v0), double(c), dopt);
info = struct('err',err,'status',status,'time',toc(t0));
if status, warning('tr_limdist:status','tr_region: the error estimate %.1e is above the tolerance; the result may be inaccurate.',err); end
end
