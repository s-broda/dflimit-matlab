function [f,info] = adfmax_pdf(x,cas,c,opts)
%ADFMAX_PDF  Limiting density of the ADFmax test of Leybourne (1995); arguments as in ADFMAX_CDF.
%
%   See also ADFMAX_CDF, ADFMAX_ICDF.
if nargin < 3 || isempty(c), c = 0; end
if nargin < 4 || isempty(opts), opts = struct; end
[f,info] = adfmax_dist(x,cas,c,opts,'pdf');
end
