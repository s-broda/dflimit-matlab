function D = dfdet_fdesign(f,breaks,opts)
%DFDET_FDESIGN  Deterministic terms given by a function, as a design of the general engine.
%
%   D = DFDET_FDESIGN(f,breaks,opts): f(r) returns the vector of the K
%   regressors at r, smooth on each piece between the breaks (0 < breaks < 1,
%   default none; jumps and kinks only there; left-continuous at a break, as in
%   DFDET_DETERM).  On each piece every regressor is replaced by its Chebyshev
%   interpolant (first kind) of the smallest degree 2^m - 1 >= 15 at which the
%   coefficients beyond the last that exceeds tol times the largest have
%   fallen below it (an error beyond maxdeg).  For the covariance the
%   transfer-function formula needs Taylor coefficients, whose sum over the
%   whole piece is ill conditioned at high degree (the monomials of
%   sin(4 pi r) on [0, 1] reach 3e4); the piece is therefore cut into sub
%   computational sub-pieces (default ceil(deg/4)) and the interpolant is
%   re-expanded about the left end of each, where the terms g^(k) h^k / k!
%   stay of order one.  The design carries these local coefficients (the field
%   loc).  Options: K (the number of regressors, default numel(f(0.5))), name,
%   tol (1e-14), maxdeg (96), sub ([]: the default above).
%
%   Example: the Fourier terms of Enders and Lee, one frequency
%      D = dfdet_fdesign(@(r) [1, sin(2*pi*r), cos(2*pi*r)]);
%      dfdet_icdf(0.05, D)
%
%   See also DFDET_DETERM, DFDET_LTREND, DFDET_DESIGN_FROM_MATRIX.
if nargin < 2, breaks = []; end
if nargin < 3 || isempty(opts), opts = struct; end
if ~isa(f,'function_handle'), error('dfdet:args','dfdet_fdesign: f must be a function handle.'); end
K = getf(opts,'K',numel(f(0.5))); name = getf(opts,'name','function design');
tol = getf(opts,'tol',1e-14); maxdeg = getf(opts,'maxdeg',96); sub = getf(opts,'sub',[]);
bk = [0, sort(double(breaks(:).')), 1];
if ~all(diff(bk) > 0), error('dfdet:args','dfdet_fdesign: breaks must lie strictly inside (0, 1) and be distinct.'); end
np0 = numel(bk) - 1;
pa = zeros(1,np0); pb = pa; pdeg = pa; pcs = cell(1,np0);
for j = 1:np0
    a = bk(j); b = bk(j+1);
    N = 16;
    while true
        Cc = chebcoef(f,a,b,N,K);
        mx = max(abs(Cc(:)));
        if mx == 0
            last = 0;
        else
            last = find(any(abs(Cc) > tol*mx,2),1,'last');
        end
        if last <= N - 6 || N > maxdeg
            deg = last - 1;
            break
        end
        N = N*2;
    end
    if N > maxdeg && deg > maxdeg - 6
        error('dfdet:args','dfdet_fdesign: on piece [%g, %g] the Chebyshev coefficients do not fall below %g by degree %d; the regressors must be smooth between the breaks.',a,b,tol,maxdeg);
    end
    deg = max(deg,0);
    pa(j) = a; pb(j) = b; pdeg(j) = deg; pcs{j} = Cc(1:deg+1,:);
end
Dg = max(pdeg);
br = zeros(1,0); Cs = {};
for j = 1:np0
    a = pa(j); b = pb(j); deg = pdeg(j);
    if isempty(sub), m = max(1,ceil(deg/4)); else, m = sub; end
    for q = 0:m-1
        r0 = a + (b - a)*q/m;
        br(end+1) = r0; %#ok<AGROW>
        Cm = zeros(Dg + 1,K);
        for i = 1:K
            Cm(1:deg+1,i) = cheb_taylor(pcs{j}(:,i),a,b,r0,deg);
        end
        Cs{end+1} = Cm; %#ok<AGROW>
    end
end
br(end+1) = 1;
D = dfdet_determ(br,Cs,name,true);
end

function Cc = chebcoef(f,a,b,N,K)
% Chebyshev coefficients (degree N - 1) of the interpolant of the columns of f on [a, b] at the N points of the first
% kind (interior points only, so that f is never evaluated at a break)
x = cos(pi*((0:N-1) + 0.5)/N);
V = zeros(K,N);
for j = 1:N
    v = f((a + b)/2 + (b - a)/2*x(j));
    if numel(v) ~= K, error('dfdet:args','dfdet_fdesign: the function returns %d values, expected %d.',numel(v),K); end
    V(:,j) = double(v(:));
end
Cc = zeros(N,K);
jj = (1:N) - 0.5;
for k = 0:N-1
    Cc(k+1,:) = (2/N)*(V*cos(k*pi*jj/N).').';
end
Cc(1,:) = Cc(1,:)/2;
end

function d = chebder(c)
% the derivative of a Chebyshev series (coefficients c, variable x in [-1, 1])
n = numel(c) - 1;
if n == 0, d = 0; return, end
d = zeros(n + 1,1);
for k = n:-1:1
    if k + 1 <= n, dk2 = d(k+2); else, dk2 = 0; end
    d(k) = dk2 + 2*k*c(k+1);
end
d(1) = d(1)/2;
d = d(1:n);
end

function v = chebval(c,x)
% Clenshaw evaluation of a Chebyshev series at x
b1 = 0; b2 = 0;
for k = numel(c):-1:2
    t = 2*x*b1 - b2 + c(k); b2 = b1; b1 = t;
end
v = x*b1 - b2 + c(1);
end

function out = cheb_taylor(c,a,b,r0,n)
% the Taylor coefficients g^(k)(r0)/k!, k = 0..n, of the Chebyshev series c on [a, b] at r0 (the derivative series and
% Clenshaw, scaled by (2/(b - a))^k / k!)
x0 = (2*r0 - a - b)/(b - a); sc = 2/(b - a);
out = zeros(n + 1,1); cd = c(:); fct = 1;
for k = 0:n
    out(k+1) = chebval(cd,x0)*fct;
    cd = chebder(cd); fct = fct*sc/(k + 1);
end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
