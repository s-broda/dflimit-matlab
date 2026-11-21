function [P,L] = tr_lincond(x,a,cas,c,cbar)
%TR_LINCOND  P(U - aV <= x - 1/2) with a zero start by the conditional integral given R = U + 1/2.
%
%   P = TR_LINCOND(x,a,cas,c,cbar) (elementwise in x, a finite) in case 'DK' or 'DG' (detrended at cbar, default
%   -13.5), clamped to [0, 1].  With a zero start the joint moment generating function is E[e^(sR + tV)] =
%   B(t) (1 - s/s1(t))^(-1/2): R is a gamma(1/2) variable of rate kappa = s1(0), and the conditional moment
%   generating function of V given R = r is M(t | r) = B(t) (s1(t)/kappa)^(1/2) e^(-(s1(t) - kappa) r), so that
%     a > 0:  P(U - aV <= x - 1/2) = P(R <= x) + int_{r > x} f_R(r) P(V >= (r - x)/a | R = r) dr,
%     a < 0:  P(U - aV <= x - 1/2) = P(R <= x) - int_{r < x} f_R(r) P(V > (x - r)/|a| | R = r) dr  (0 for x <= 0),
%   P(R <= x) = erf(sqrt(kappa x)), exactly so at a = 0.  P(V > y | r) by Gil-Pelaez, the midpoint rule in w with the
%   step 2 pi/(2.2 ymax) beyond the range ymax = E V + 42/t_* of y, up to where the transform has fallen below
%   1e-17 w; the outer integral on panels graded from its lower end in y (the first E V/16 long, doubling) and at the
%   scale of f_R (rho = k/sqrt(kappa)), 24-point Gauss-Legendre in rho = sqrt(r) on each.  At most 2e6 nodes in w
%   (case DG for c up to about 3.5 at cbar = -13.5; an error beyond).  The engine of LTRATIO_LINCDF in case DG
%   where the inversion of Gil-Pelaez is not trusted (TR_GILPELAEZ, TR_LINEXACT).  [~,L] = TR_LINCOND([],a,cas,c,cbar)
%   returns the setup (L.kap the rate kappa, L.ev = E V).  As lincond of julia/DFExact.jl.
%
%   See also TR_LINEXACT, TR_GILPELAEZ, LTRATIO_LINCDF.
if nargin < 5 || isempty(cbar), cbar = -13.5; end
L = lc_setup(upper(char(cas)),c,cbar);
P = zeros(size(x));
for k = 1:numel(x), P(k) = lc_x(L,x(k),a); end
end

function L = lc_setup(cas,c,cbar)
% the setups of the last few (case, c, cbar), cached
persistent C
if isempty(C), C = struct('key',{},'L',{}); end
key = sprintf('%s %.17g %.17g',cas,c,cbar);
for j = 1:numel(C), if strcmp(C(j).key,key), L = C(j).L; return, end, end
L = lc_build(cas,c,cbar);
C = [struct('key',key,'L',L) C(1:min(end,7))];
end

function L = lc_build(cas,c,cbar)
LCNMAX = 2e6;
[~,s1] = lc_parts(cas,c,0,cbar); kap = real(s1);
[lb,s1] = lc_parts(cas,c,1e-6,cbar); ev = imag(lb - log(s1)/2)/1e-6;   % E V: log B(iw) ~ i w E V
if strcmp(cas,'DK'), ts = dfx_tstar(1,c); else, ts = tstar_gls(c,cbar); end
if ~(ts > 0), error('ltratio:args','ltratio_lincdf: case %s, c = %g: no zero of the moment generating function of V was found (t_*).',cas,c); end
ymax = ev + 42/ts; dw = 2*pi/(2.2*ymax);
W = 1;
while true
    lb = lc_parts(cas,c,W,cbar);
    if real(lb) - log(kap)/2 - log(W) < log(1e-17), break, end
    W = W*1.3;
    if W/dw > LCNMAX
        error('ltratio:args',['ltratio_lincdf: case %s, c = %g: the conditional integral of the linear form needs more than %d nodes in its ' ...
            'inner inversion (the law of V spreads like 1/t_* = %.3g as c grows); case DG is available for c up to about 3.5 at cbar = -13.5, ' ...
            'less for cbar nearer 0.'],cas,c,LCNMAX,1/ts);
    end
end
nw = ceil(W/dw); w = ((1:nw)' - 0.5)*dw;
[lB,s1] = lc_parts(cas,c,w,cbar);
L = struct('kap',kap,'ev',ev,'ymax',ymax,'dw',dw,'w',w,'lB',lB - log(kap)/2,'s1k',s1 - kap);
end

function [lB,s1] = lc_parts(cas,c,w,cbar)
% log(B s1^(1/2)) and s1 at t = i w (a column), the phases continued along w
w = w(:);
if strcmp(cas,'DK')
    mu = sqrt(complex(c^2,-2*w));
    q = -expm1(-2*mu)./mu; sm = abs(mu) < 1e-6; q(sm) = 2 - 2*mu(sm);
    s1 = ((1 + exp(-2*mu)) - c*q)./q;
    lq = real(mu) - log(2) + log(abs(q)); ph = unwrapc(imag(mu) + angle(q));
    lB = -c/2 - complex(lq,ph)/2;
else
    n = numel(w); Dv = complex(zeros(n,1)); s1 = Dv;
    for j0 = 1:1e5:n                                   % (in blocks: gls_parts holds some 30 arrays of its argument's size)
        j = (j0:min(n,j0 + 1e5 - 1))';
        [lch,D,s] = gls_parts(1i*w(j),c,cbar); Dv(j) = exp(lch).*D; s1(j) = s;
    end
    ph = unwrapc(angle(Dv)); ps = unwrapc(angle(s1));
    lB = -c/2 - complex(log(abs(Dv)),ph)/2 + complex(log(abs(s1)),ps)/2;
end
end

function ph = unwrapc(ph)
% the phase continued along the grid (ph_j - 2 pi round((ph_j - ph_{j-1})/(2 pi)), as unwrap! of DFExact.jl)
for j = 2:numel(ph), ph(j) = ph(j) - 2*pi*round((ph(j) - ph(j-1))/(2*pi)); end
end

function ts = tstar_gls(c,cbar)
% the abscissa t_* of the moment generating function of V in case DG detrended at cbar: the first zero of
% cosh th |I - Q Sig| at real t > 0 (tstar_gls of DFExact.jl)
g = @(t) gfun(t,c,cbar);
a = 1e-6;
while g(a) <= 0 && a > 1e-300, a = a/1e3; end
b = a;
while g(b) > 0
    a = b; b = b*1.05;
    if b > max(1e4,2*c^2), ts = NaN; return, end
end
for it = 1:100
    m = (a + b)/2;
    if g(m) > 0, a = m; else, b = m; end
    if b - a < 1e-15*b, break, end
end
ts = (a + b)/2;
end
function v = gfun(t,c,cbar)
[lch,D] = gls_parts(complex(t),c,cbar); v = real(exp(lch).*D);
end

function sv = lc_sv(L,y,r)
% P(V > y | R = r), elementwise in (y, r) (columns), in blocks of the grid
sv = ones(size(y)); k = find(y > 0); nb = max(1,floor(4e6/numel(L.w)));
for j0 = 1:nb:numel(k)
    j = k(j0:min(end,j0 + nb - 1));
    E = exp(L.lB - L.s1k*r(j).' - 1i*(L.w*y(j).'));
    sv(j) = 0.5 + (sum(imag(E)./L.w,1)).'*L.dw/pi;
end
end

function acc = lc_outer(L,x,a,y0,y1,zlo,zhi)
% the outer integral of |a| f_R(x + a y) P(V > y | R = x + a y) over [y0, y1]; zlo, zhi: r = 0 at y0, at y1 (taken
% exactly: the root of x + a y0 rounded off leaves r ~ 1e-18 there, and sqrt(r) ~ 1e-9 would cut off f_R(r) dr ~ 1e-9)
acc = 0;
if y1 <= y0, return, end
[u,gw] = glrule(24); u = u(:); gw = gw(:);
br = y0; h = L.ev/16;
while br(end) < y1, br(end+1) = min(y1,br(end) + h); h = h*2; end %#ok<AGROW>
rh = sqrt(max(x + a*br,0));
if zlo, rh(1) = 0; end
if zhi, rh(end) = 0; end
lo = min(rh); hi = max(rh); sk = 1/sqrt(L.kap);
for m = [0.25 0.5 1 1.5 2 3 4 6]
    if lo < m*sk && m*sk < hi, rh(end+1) = m*sk; end %#ok<AGROW>
end
rh = unique(rh);
% on each panel Gauss-Legendre with 24 points, checked against 12; a panel where the two differ by more than 1e-13 is
% halved (level by level, at most 30 times, at most 4096 panels a level): the conditional law of V given R = r can turn
% within a fraction of the scale 1/kappa of r (case DG for c > 0 detrended at cbar near 0), which fixed panels do not
% resolve
[u2,gw2] = glrule(12); u2 = u2(:); gw2 = gw2(:);
pa = rh(1:end-1); pb = rh(2:end);
for depth = 0:30
    if isempty(pa), break, end
    rho = [pa + (pb - pa).*u; pa + (pb - pa).*u2];      % 36 x (panels)
    r = rho(:).^2;
    f = reshape(2*sqrt(L.kap/pi)*exp(-L.kap*r).*lc_sv(L,(r - x)/a,r),36,[]);
    s24 = (pb - pa).*sum(gw.*f(1:24,:),1); s12 = (pb - pa).*sum(gw2.*f(25:36,:),1);
    ok = abs(s24 - s12) <= 1e-13 | depth == 30;
    bad = find(~ok); if numel(bad) > 2048, ok(bad(2049:end)) = true; bad = bad(1:2048); end
    acc = acc + sum(s24(ok));
    m = (pa(bad) + pb(bad))/2;
    pa = reshape([pa(bad); m],1,[]); pb = reshape([m; pb(bad)],1,[]);
end
end

function F = lc_x(L,x,a)
if x > 0, F0 = erf(sqrt(L.kap*x)); else, F0 = 0; end
if a == 0 || isnan(x), F = F0; if isnan(x), F = NaN; end, return, end
Y = L.ymax;
if a > 0
    y0 = max(0,-x/a);
    if y0 >= Y, F = F0; else, F = F0 + lc_outer(L,x,a,y0,Y,x <= 0,false); end
else
    y1 = x/(-a);
    if x <= 0, F = 0; else, F = F0 - lc_outer(L,x,a,0,min(Y,y1),false,y1 <= Y); end
end
F = min(1,max(0,F));
end
