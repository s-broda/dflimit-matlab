function f = pt_pdf(q,D,cbar,c,opts)
%PT_PDF  Limiting density of the point-optimal statistic against cbar for the design D (arguments as in PT_CDF).
%
%   Example: pt_pdf(6.789, dfdet_perron('B', 0.5), -17.8)
%
%   See also PT_CDF.
if nargin < 4 || isempty(c), c = 0; end
if nargin < 5 || isempty(opts), opts = struct; end
G = dfdet_qgrid(D,double(cbar),double(c),getf(opts,'N',60000),getf(opts,'umax',4000));
f = zeros(size(q));
for k = 1:numel(q), f(k) = sum(real(exp(G.lM - 1i*G.u*q(k))).*G.w)/pi; end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
