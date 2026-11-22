function D = dfdet_design_from_matrix(X,opts)
%DFDET_DESIGN_FROM_MATRIX  The design of the limit from the matrix of the deterministic regressors of a regression.
%
%   D = DFDET_DESIGN_FROM_MATRIX(X,opts): X the n x K matrix of the regressors
%   (column k sampled at t = 1..n, any scaling).  opts.breaks, if given, are
%   the indices T_B after which the regressors change (a regime of t <= T_B
%   and one of t > T_B, as 1{t > T_B}), at lambda = T_B/n.  Without them the
%   breaks are detected: on each column scaled to unit maximum the sixth
%   differences of a function smooth on the scale of the sample are of order
%   n^-6, while a jump or a kink at T leaves a cluster of differences of order
%   one or 1/n around T; each cluster gives the largest T (after which the
%   regime changes) at which polynomials of degree at most 8 fit 30
%   observations on either side exactly.  On each piece every column is then
%   fitted by least squares in Chebyshev polynomials of increasing degree up to
%   opts.maxdeg (40) until the largest residual is below opts.tol (1e-8)
%   (exact for polynomial and dummy columns, a smooth interpolant for others),
%   and the design is the function design of these fits (DFDET_FDESIGN).
%   Columns that are not piecewise smooth at the sample's resolution (a piece
%   shorter than opts.minlen = 10 observations, as a pulse dummy makes, or no
%   fit within maxdeg, as seasonal dummies give) are refused with an error:
%   they do not affect the limit and should be left out.  Columns linearly
%   dependent on the others are dropped, with a warning.  opts.name names it.
%
%   Example: n = 200; t = (1:n)';
%      D = dfdet_design_from_matrix([ones(n,1), t, t > 60]);   % a level shift after t = 60
%      dfdet_icdf(0.05, D)
%
%   See also DFDET_FDESIGN, DFDET_PERRON.
if nargin < 2 || isempty(opts), opts = struct; end
breaks = getf(opts,'breaks',[]); name = getf(opts,'name','design from a matrix');
tol = getf(opts,'tol',1e-8); maxdeg = getf(opts,'maxdeg',40); minlen = getf(opts,'minlen',10);
if ~(isnumeric(X) || islogical(X)) || ~ismatrix(X), error('dfdet:args','dfdet_design_from_matrix: X must be a numeric matrix.'); end
[n,K] = size(X);
Xs = double(X);
for k = 1:K
    m = max(abs(Xs(:,k))); if m > 0, Xs(:,k) = Xs(:,k)/m; end
end
% dependent columns
[~,R,E] = qr(Xs,0);
rk = nnz(abs(diag(R)) > 1e-9*abs(R(1,1)));
keep = sort(E(1:rk));
if rk < K
    warning('dfdet:design','dfdet_design_from_matrix: columns %s are linearly dependent on the others and are dropped.',mat2str(setdiff(1:K,keep)));
end
Xs = Xs(:,keep); K = rk;
if isempty(breaks)
    q = 6; T = zeros(1,0);
    flagged = false(n,1);
    for k = 1:K
        d = Xs(:,k);
        for it = 1:q, d = diff(d); end
        flagged(find(abs(d) > 1e-7)) = true; %#ok<FNDSB>
    end
    idx = find(flagged); cl = {};
    for t = idx(:).'
        if isempty(cl) || t - cl{end}(end) > q
            cl{end+1} = t; %#ok<AGROW>
        else
            cl{end}(end+1) = t;
        end
    end
    if numel(cl) > max(10,floor(n/50))
        error('dfdet:args','dfdet_design_from_matrix: %d clusters of jumps or kinks: some column is not piecewise smooth at the sample''s resolution (seasonal dummies?); such columns do not affect the limit and should be left out.',numel(cl));
    end
    w = 30;
    for ic = 1:numel(cl)
        c = cl{ic}; best = 0;
        for Tc = max(c(1) - 1,1):min(c(end) + q + 1,n - 1)
            lo = max(1,Tc - w + 1); hi = min(n,Tc + w);
            if Tc - lo + 1 < 3 || hi - Tc < 3, continue, end
            ok = true;
            for k = 1:K
                [~,r1] = pfit(Xs(lo:Tc,k),min(8,Tc - lo)); [~,r2] = pfit(Xs(Tc+1:hi,k),min(8,hi - Tc - 1));
                if ~(min(r1,1) <= 1e-9 && min(r2,1) <= 1e-9), ok = false; break; end
            end
            if ok, best = Tc; end
        end
        if best == 0
            error('dfdet:args','dfdet_design_from_matrix: no break position fits the cluster of changes at t = %d..%d (a pulse, or a column that is not piecewise smooth?); such columns do not affect the limit and should be left out.',c(1),c(end));
        end
        T(end+1) = best; %#ok<AGROW>
    end
    breaks = T;
end
TB = sort(round(double(breaks(:).')));
bk = [0 TB n];
if any(diff(bk) < minlen)
    error('dfdet:args','dfdet_design_from_matrix: a piece shorter than %d observations (breaks at %s): a pulse or a feature at the sample''s resolution; such columns do not affect the limit and should be left out.',minlen,mat2str(TB));
end
coefs = cell(numel(bk) - 1,K);
for j = 1:numel(bk)-1
    t1 = bk(j) + 1; t2 = bk(j+1);
    for k = 1:K
        y = Xs(t1:t2,k); got = false;
        for deg = 0:min(maxdeg,t2 - t1)
            [cf,res] = pfit(y,deg);
            if res <= tol, coefs{j,k} = cf; got = true; break; end
        end
        if ~got
            error('dfdet:args','dfdet_design_from_matrix: column %d is not smooth on observations %d..%d at the sample''s resolution (no Chebyshev fit of degree <= %d within %g); such columns do not affect the limit and should be left out.',keep(k),t1,t2,maxdeg,tol);
        end
    end
end
lams = TB/n;
D = dfdet_fdesign(@(r) fval(r,lams,TB,bk,n,coefs,K),lams,struct('K',K,'name',name));
end

function v = fval(r,lams,TB,bk,n,coefs,K)
% the fitted polynomials in the variable of piece j: t in [t1, t2] maps to x in [-1, 1]; t = n r
j = 1; while j <= numel(TB) && r > lams(j), j = j + 1; end
t1 = bk(j) + 1; t2 = bk(j+1);
if t2 == t1, x = 0; else, x = (2*(n*r) - t1 - t2)/(t2 - t1); end
v = zeros(K,1);
for k = 1:K, v(k) = chebval(coefs{j,k},x); end
end

function [c,res] = pfit(y,deg)
% least squares in the Chebyshev polynomials of degree <= deg on equispaced points of [-1, 1]
x = linspace(-1,1,numel(y))';
V = cos(acos(min(max(x,-1),1))*(0:deg));
c = V\y;
res = max(abs(V*c - y));
end

function v = chebval(c,x)
b1 = 0; b2 = 0;
for k = numel(c):-1:2
    t = 2*x*b1 - b2 + c(k); b2 = b1; b1 = t;
end
v = x*b1 - b2 + c(1);
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
