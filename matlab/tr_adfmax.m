function [F,info] = tr_adfmax(x,cas,c,opts,mode)
%TR_ADFMAX  The engine of ADFMAX_CDF and ADFMAX_PDF: the limiting distribution of the ADFmax statistic.
%
%   [F,info] = TR_ADFMAX(x,cas,c,opts,mode), mode 'cdf' (default) or 'pdf', x an array of finite values: the limit
%   (|R| - 1/2)/sqrt(V), R = U + 1/2, of the larger of the Dickey-Fuller t-ratios of the forward and the reverse
%   regression (Leybourne 1995), so that {ADFmax <= x} = {|R| <= a(V)}, a(v) = 1/2 + x sqrt(v), empty where a < 0,
%   i.e. beyond v0 = 1/(4x^2) for x < 0.  By the conditional integral of TR_LIMDIST with that boundary (opts.bnd,
%   bnd.two = true): cases 'DU', 'DT' with a zero start, the Girsanov forms of the inner inversion (G(a, b) + G(-a, b)),
%   cases 'SU', 'ST' with a stationary start (c <= -1e-4), the inner inversion of the general route, with the integration by
%   parts of TR_LIMDIST and the outer rule cosine-mapped twice at v0.  info: per point the error estimate, the status
%   and the time.  The C engine has a drop-in of this function (c/mex).
%
%   See also ADFMAX_CDF, ADFMAX_PDF, ADFMAX_ICDF, TR_LIMDIST.
if nargin < 4 || isempty(opts), opts = struct; end
if nargin < 5 || isempty(mode), mode = 'cdf'; end
F = zeros(size(x)); info = struct('err',zeros(size(x)),'status',zeros(size(x)),'time',zeros(size(x)));
o = opts; if isfield(o,'method'), o = rmfield(o,'method'); end
o.c0 = 0; o.g0 = 0; o.verified = true;                % (the applications keep the agreement rule)
for k = 1:numel(x)
    xk = x(k);
    bnd = struct('a',@(b) 1/2 + xk*sqrt(b),'da',@(b) xk./(2*sqrt(b)),'two',true,'v0',Inf);
    if xk < 0, bnd.v0 = 1/(4*xk^2); end
    o.bnd = bnd;
    [F(k),inf1] = tr_limdist(0,cas,c,o,mode);
    info.err(k) = inf1.err; info.status(k) = inf1.status; info.time(k) = inf1.time;
end
end
