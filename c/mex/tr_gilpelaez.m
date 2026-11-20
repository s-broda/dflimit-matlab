function [P,phi] = tr_gilpelaez(b,a,cas,c,opts,phi)
%TR_GILPELAEZ  P(U - aV < b) with a zero start by the C engine (MEX file LTRATIO_MEX): a drop-in replacement for the
%   MATLAB engine of the same name (cases 'DK' and 'DG', a = NaN: P(V < b) in case DK; the options cbar (-13.5),
%   N (60000) and umax (4000)), so that LTRATIO_LINCDF and POWER_ENVELOPE run on the C engine (dfl_lincdf in
%   c/dfl_gilpelaez.c).  The characteristic function is not returned (phi = []), and a sixth argument is ignored:
%   the engine computes it anew at every call.
%
%   See also LTRATIO_MEX, LTRATIO_LINCDF, POWER_ENVELOPE.
if nargin < 5 || isempty(opts), opts = struct; end
switch upper(char(cas))
    case 'DK', model = 1;
    case 'DG', model = 4; if isnan(a), error('tr_gilpelaez:args','tr_gilpelaez: P(V < b) is available in case DK only.'); end
    otherwise, error('tr_gilpelaez:args','tr_gilpelaez: the case must be DK or DG.');
end
P = ltratio_mex('lincdf', double(b), double(a), model, double(c), getf(opts,'cbar',-13.5), [getf(opts,'N',NaN) getf(opts,'umax',NaN)]);
P = reshape(P, size(b)); phi = [];
end
function v = getf(o,f,d), if isfield(o,f) && ~isempty(o.(f)), v = double(o.(f)); else, v = d; end, end
