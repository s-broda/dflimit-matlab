function [F,info] = tr_adfmax(x,cas,c,opts,mode)
%TR_ADFMAX  The ADFmax statistic by the C engine (MEX file LTRATIO_MEX): a drop-in replacement for the MATLAB engine of
%   the same name (the larger of the forward and reverse Dickey-Fuller t-ratios, cases 'DU', 'DT', 'SU', 'ST'; mode
%   'cdf' or 'pdf'; the options M, Nv, Nl, maxlev, plen, otol, stepw, errtol, blo, bhi, bmin, mumin, tdec), so that
%   ADFMAX_CDF, ADFMAX_PDF and ADFMAX_ICDF run on the C engine.  x finite.
%
%   See also LTRATIO_MEX, TR_LIMDIST, TR_REGION.
if nargin < 4 || isempty(opts), opts = struct; end
if nargin < 5 || isempty(mode), mode = 'cdf'; end
cas = upper(char(cas));
if ~any(strcmp(cas,{'DU','DT','SU','ST'})), error('tr_adfmax:args','tr_adfmax: the case must be DU, DT, SU or ST.'); end
kind = 2; if strcmpi(mode,'pdf'), kind = 1; end
model = find('KUT' == cas(2)); stat = double(cas(1) == 'S');
t0 = tic;
[F,err,status] = ltratio_mex('adfmax', kind, double(x), model, double(c), stat, dopts(opts));
info = struct('err',err,'status',status,'time',toc(t0)*ones(size(x))/max(1,numel(x)));
if any(status(:)), warning('tr_limdist:status','tr_adfmax: %d of %d points with nonzero status; the result may be inaccurate there.',nnz(status),numel(x)); end
end
function v = dopts(o)
nm = {'M','Nv','Nl','maxlev','plen','otol','stepw','errtol','blo','bhi','bmin','mumin','tdec'}; v = NaN(1,numel(nm));
for k = 1:numel(nm), if isfield(o,nm{k}) && ~isempty(o.(nm{k})), v(k) = double(o.(nm{k})); end, end
end
