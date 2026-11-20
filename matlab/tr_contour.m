function [F,info] = tr_contour(x,cas,c,opts,mode)
%TR_CONTOUR  Limiting distribution and density of the t-ratio with a fixed start: the contour integral.
%
%   [F,info] = TR_CONTOUR(x,cas,c,opts,mode), mode 'cdf' (default) or 'pdf': the
%   limit as n -> Inf of P[t <= x] (or its density) under phi = 1 + c/n, no drift,
%   with a fixed start y0 = g0 sigma sqrt(n), g0 ~= 0 (opts.g0, required), in
%   cases 'DU' and 'DT' only; t is the t-ratio for the null value phi0 = 1 + c0/n
%   (opts.c0, default 0, as in TR_LIMDIST and LTRATIO_CDF).  x may be an array.
%
%   The method is the contour integral representation: the inversion
%   formula on paths around the singularities of the moment generating
%   function, with the integral over v in closed form and the paths in s and
%   t deformed,
%   a double integral with an explicit integrand (a port, step by step, of
%   tr_contour of julia/DFExact.jl, November 2026):
%      f(x) = (1/2 pi i) int_C I_f(t) dt,   F(x) = (1/2 pi i) int_C I_F(t) dt/t,
%   I(t) = (1/2 pi i) [int_{G+} Mt W+ ds + int_{G-} Mt W- ds], with
%   Mt(s,t) = P(t) [(s1 - s)(s - s2)]^(-1/2) exp{g0^2 Q(s,t)}, s1, s2 the roots
%   of D(., t), Q the initial-value term of the joint moment generating
%   function with a fixed start, and W+-, the integrals
%   over the parts of v > 0 where X(v) = c0 v^2 + x v + 1/2 is positive and
%   negative of 2 v^2 e^(-t v^2 - s X(v)) (density) or X'(v) e^(-t v^2 - s X(v))
%   (distribution function), incomplete Gaussian integrals in closed form (the
%   functions h_k from the complex erfcx, FADDEEVA, for |z| < 6 and their
%   asymptotic series beyond); the distribution function adds the kernels
%   (1 - e^(-s/2))/s on G+ and -1/s on G-.  The outer contour C is the
%   hyperbola t(u) = tau + a (cosh u - 1) + i b sinh u (tau the line of the
%   line integral, its angle 45 degrees for c >= -10, rising to 70 degrees at
%   c = -40), integrated by the trapezoidal rule in u up to |t - tau| = 1e5 or
%   until a block of 32 nodes whose terms all stay below 1e-18; the branches
%   of P(t) and the labels of the roots are continued along C on eight points
%   per step.  At each node t the inner paths are U-shapes around the roots,
%   G+ around s1 opening to the right and G- around s2 opening to the left,
%   whose vertical segment is placed by minimizing a bound of the largest
%   term of the quadrature (a scan of 24 distances, then 12 steps of golden
%   section), graded by Gauss-Legendre panels toward the origin, the other
%   root and the branch point -t/c0, and continued by two horizontal rays.
%   The error estimate comes from two structurally different evaluations: A
%   as above, with 64 nodes per panel and ray, and B, with the asymptotic
%   angle of the hyperbola larger by 7.5 degrees (so that every t differs),
%   the grading cost weighted by 1 instead of 0.5 in the placement (so that
%   the U-shapes and their panels differ) and 48 nodes; the estimate is
%   |A - B| plus 10 eps times the sums of the largest terms of both, and while
%   it exceeds 1e-10 the next level (A: 96, 128, then 192 nodes; B: 64, 96,
%   128) is added to both at all their nodes.  The value is A at the last
%   level.
%
%   opts.verified (false by default): false is one evaluation, A at its first
%   level (info.err NaN, info.status 0, nodes those of A); true the pair A, B
%   and the levels described above, at 1.5 to 3 times the cost, with
%   info: time; err, the error estimate per point; status, 0 where the estimate
%   is at most 1e-10 and the value is finite, otherwise 2 (the warning
%   tr_contour:status); nodes, the number of outer nodes of A and B together;
%   level, the level of the value (1 to 4).  Non-finite x: F = 0 at -Inf, 1 (cdf) or 0
%   (pdf) at +Inf, NaN at NaN.  The engine pieces of the line integral (the
%   roots, P(t), t_* and the line tau) come from DFX_CORE.
%
%   See also LTRATIO_CDF, LTRATIO_PDF, TR_LIMDIST, DFX_CORE.
if nargin < 4 || isempty(opts), opts = struct; end
if nargin < 5 || isempty(mode), mode = 'cdf'; end
msg = 'tr_contour: the contour integral applies to a fixed start g0 ~= 0 in cases DU, DT.';
if ~(ischar(cas) || isstring(cas)), error('tr_contour:args',msg); end
switch upper(char(cas))
    case 'DU', model = 2;
    case 'DT', model = 3;
    otherwise, error('tr_contour:args',msg);
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('tr_contour:args','tr_contour: c must be a finite real scalar.'); end
if ~(isnumeric(x) && isreal(x)), error('tr_contour:args','tr_contour: x must be real.'); end
if ~isstruct(opts), error('tr_contour:args','tr_contour: opts must be a structure.'); end
g0 = getf(opts,'g0',0); c0 = getf(opts,'c0',0);
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('tr_contour:args','tr_contour: opts.g0 must be a finite real scalar.'); end
if g0 == 0, error('tr_contour:args',msg); end
if ~(isnumeric(c0) && isscalar(c0) && isreal(c0) && isfinite(c0)), error('tr_contour:args','tr_contour: opts.c0 must be a finite real scalar.'); end
switch mode
    case 'cdf', kind = 2;
    case 'pdf', kind = 1;
    otherwise, error('tr_contour:args','tr_contour: mode must be ''cdf'' or ''pdf''.');
end
c = double(c); c0 = double(c0); g0 = double(g0);
verified = logical(getf(opts,'verified',false));
t0 = tic;
F = zeros(size(x)); err = F; status = F; nodes = F; level = F;
for k = 1:numel(x)
    xk = double(x(k));
    if ~isfinite(xk)
        if isnan(xk), F(k) = NaN; elseif xk > 0 && kind == 2, F(k) = 1; end
        continue
    end
    [F(k),err(k),status(k),nodes(k),level(k)] = ctr_value(kind,xk,model,c,c0,g0,verified);
end
info = struct('time',toc(t0),'err',err,'status',status,'nodes',nodes,'level',level);
if verified && any(status(:))
    iw = find(status,1);
    warning('tr_contour:status',['tr_contour: %d of %d points with nonzero status (first: x = %.6g, error estimate ' ...
        '%.3g); the result may be inaccurate there.'],nnz(status),numel(x),x(iw),err(iw));
end
end

% =====================================================================================================================
% The constants of julia/DFExact.jl: CTR_LEVELS (64, 96, 128, 192; B: CTR_LEVELS_B = (48, 64, 96, 128)), CTR_NGRID (24),
% CTR_MARG (1e-3), CTR_PEN (0.5; B: CTR_PEN_B = 1), CTR_DPHI (pi/24), CTR_RMAX (1e5), CTR_SUB (8), CTR_ERRTOL (1e-10);
% the outer nodes in chunks of 32.
% The functions below follow ctr_* there one by one; where Julia loops over the outer nodes of a chunk (on several
% threads) or over the nodes of a path, these evaluate them together, with the same arithmetic element by element.

function [v,err,st,nodes,l] = ctr_value(kind,x,model,c,c0,g0,verified)
% f(x) (kind 1) or F(x) (kind 2): the value, the error estimate, the status, the number of outer nodes, the level.  Two
% structurally different evaluations, A (the hyperbola at the angle ctr_phi(c), placement weight 0.5, 64 nodes per
% panel and ray) and B (the angle larger by pi/24, placement weight 1, one level fewer nodes, 48); the error estimate
% is |A - B| plus 10 eps times the sums of the largest terms of both; while it exceeds 1e-10 the next level is added
% to both; the value is A at the last level
NL = 4; ERRTOL = 1e-10;
A = ctr_eval(kind,x,model,c,c0,g0,1);
if ~verified                                    % one evaluation: A at the first level, no estimate
    l = 1; v = ctr_sum(A,1); err = NaN; st = 0; if ~isfinite(v), st = 2; end; nodes = A.done; return
end
B = ctr_eval(kind,x,model,c,c0,g0,2);
canc = 10*eps*(sum(A.w.*A.gm(1:A.done)) + sum(B.w.*B.gm(1:B.done)));
l = 1; err = abs(ctr_sum(A,1) - ctr_sum(B,1)) + canc;
while err > ERRTOL && l < NL
    l = l + 1;
    A = ctr_level(A,kind,x,c0,c,g0,l); B = ctr_level(B,kind,x,c0,c,g0,l);
    err = abs(ctr_sum(A,l) - ctr_sum(B,l)) + canc;
end
v = ctr_sum(A,l);
st = 0; if err > ERRTOL || ~isfinite(v), st = 2; end
nodes = A.done + B.done;
end

function v = ctr_sum(E,l)
v = sum(E.w.*E.G(1:E.done,l));
end

function E = ctr_eval(kind,x,model,c,c0,g0,variant)
% one evaluation of the contour integral (variant 1: A, 2: B) at the first level: the outer nodes processed in chunks
% of 32 until a chunk whose terms all stay below 1e-18.  At each outer node only the pieces that are not negligible
% there take part (ctr_active); the nodes of a chunk with the same active pieces are evaluated together
if variant == 2, LEVELS = [48 64 96 128]; pen = 1; else, LEVELS = [64 96 128 192]; pen = 0.5; end
NL = numel(LEVELS); CH = 32;
C = ctr_path(model,c,variant);
[plus,minus] = ctr_pieces(x,c0);
K1 = numel(C.t);
G = zeros(K1,NL); gm = zeros(K1,1);
lay = struct('sp',zeros(K1,1),'sm',NaN(K1,1),'fp',complex(NaN(K1,3)),'fm',complex(NaN(K1,3)),'nfp',zeros(K1,1), ...
    'nfm',zeros(K1,1),'rp',zeros(K1,1),'rm',zeros(K1,1));
act = false(K1,size(plus,1) + size(minus,1));          % the active pieces of each node (plus, then minus)
done = 0;
while done < K1
    idx = (done+1:min(done+CH,K1))';
    act(idx,:) = ctr_active(plus,minus,real(C.t(idx)));
    [~,~,grp] = unique(act(idx,:),'rows','stable');
    for g = 1:max(grp)
        ii = idx(grp == g);
        [pl,mi] = ctr_sel(plus,minus,act(ii(1),:));
        L = ctr_layout(kind,C,ii,x,c0,c,g0,pl,mi,pen);
        lay.sp(ii) = L.sp; lay.sm(ii) = L.sm; lay.rp(ii) = L.rp; lay.rm(ii) = L.rm;
        lay.nfp(ii) = size(L.fp,2); lay.nfm(ii) = size(L.fm,2);
        lay.fp(ii,1:size(L.fp,2)) = L.fp; lay.fm(ii,1:size(L.fm,2)) = L.fm;
        wt = C.dt(ii); if kind == 2, wt = wt./C.t(ii); end
        [I,mx] = ctr_inner_n(kind,C,ii,x,c0,c,g0,pl,mi,LEVELS(1),L);
        G(ii,1) = imag(I.*wt); gm(ii) = mx.*abs(wt);
    end
    done = idx(end);
    if all(abs(G(idx,1)) < 1e-18) && all(gm(idx) < 1e-18), break, end
end
w = repmat(C.h/pi,done,1); w(1) = w(1)/2;
E = struct('C',C,'G',G,'gm',gm,'lay',lay,'act',act,'plus',plus,'minus',minus,'done',done,'w',w,'levels',LEVELS);
end

function E = ctr_level(E,kind,x,c0,c,g0,l)
% the next level l of an evaluation (all its nodes)
CH = 32; C = E.C; lay = E.lay;
for k1 = 1:CH:E.done
    idx = (k1:min(k1+CH-1,E.done))';
    [~,~,grp] = unique(E.act(idx,:),'rows','stable');
    for g = 1:max(grp)
        ii = idx(grp == g);
        [pl,mi] = ctr_sel(E.plus,E.minus,E.act(ii(1),:));
        L = struct('sp',lay.sp(ii),'sm',lay.sm(ii),'fp',lay.fp(ii,1:lay.nfp(ii(1))),'fm',lay.fm(ii,1:lay.nfm(ii(1))), ...
            'rp',lay.rp(ii),'rm',lay.rm(ii));
        wt = C.dt(ii); if kind == 2, wt = wt./C.t(ii); end
        I = ctr_inner_n(kind,C,ii,x,c0,c,g0,pl,mi,E.levels(l),L);
        E.G(ii,l) = imag(I.*wt);
    end
end
end

function A = ctr_active(plus,minus,tr)
% the pieces that are not negligible at the outer nodes of real part tr: Re(t) a^2 <= 50 at the start a of the piece
% (a piece whose Gaussian factor is below e^-50 from its start is dropped); rows the nodes, columns the pieces (plus,
% then minus)
A = [tr.*plus(:,1).'.*plus(:,1).' <= 50, tr.*minus(:,1).'.*minus(:,1).' <= 50];
end

function [pl,mi] = ctr_sel(plus,minus,a)
% the active pieces of a row a of ctr_active
np = size(plus,1);
pl = plus(a(1:np),:); mi = minus(a(np+1:end),:);
end

function C = ctr_path(model,c,variant)
% the outer contour of model, c (variant 1: evaluation A, the angle ctr_phi(c); variant 2: evaluation B, the angle
% ctr_phi(c) + pi/24): the nodes u_k = k h, the per-point data along it (branches continued on 8 points per step) and
% the derivative t'(u_k); cached per (model, c, variant).  s1 is the root with Re s1 > 0 at tau (s2 of TPARTS)
persistent CACHE
if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','any'); end
key = sprintf('%d|%.17g|%d',model,c,variant);
if isKey(CACHE,key), C = CACHE(key); return, end
SUB = 8; RMAX = 1e5; DPHI = pi/24;
tau = dfx_core('tau',model,c,false); ts = dfx_core('tstar',model,c,false);
if variant == 2, dp = DPHI; else, dp = 0; end
phi = ctr_phi(c) + dp; y = min(0.5,0.9*(pi/2 - phi));
b = 0.9*min(tau,ts - tau)/sin(y); a = b/tan(phi);
h = 2*pi*y/(37 + 37*tan(phi)*y);
K = ceil(asinh(RMAX/hypot(a,b))/h);
ud = (0:SUB*K)'*h/SUB;
tt = complex(tau + a*(cosh(ud) - 1),b*sinh(ud));
[P,s1e,s2e] = dfx_core('tparts',model,c,tt,false);   % (the engine's labels: s1e the root with Re < 0 at tau)
j = SUB*(0:K)' + 1; t = tt(j); Pj = P(j);
if model == 2, cd = c*c; else, cd = c^4; end
KQ = cd*exp(c)*dfx_core('gstat',model,c*c - 2*t).*(Pj.*Pj);
uk = (0:K)'*h;
C = struct('t',t,'s1',s2e(j),'s2',s1e(j),'P',Pj,'KQ',KQ,'dt',complex(a*sinh(uk),b*cosh(uk)), ...
    'h',h,'tau',tau,'a',a,'b',b,'phi',phi);
CACHE(key) = C;
end

function p = ctr_phi(c)
% the asymptotic angle of the outer hyperbola: 45 degrees for c >= -10, rising linearly to 70 degrees at c = -40
p = (45 + 25*min(max((-c - 10)/30,0),1))*pi/180;
end

function [plus,minus] = ctr_pieces(q,c0)
% the pieces [a b] of {v > 0} with X(v) = c0 v^2 + q v + 1/2 > 0 (plus) and < 0 (minus), split at the positive roots
% of X (a double root is no split)
if c0 == 0
    if q < 0, r = -0.5/q; else, r = zeros(0,1); end
else
    dsc = q*q - 2*c0;
    if dsc <= 0
        r = zeros(0,1);
    else
        r = [(-q - sqrt(dsc))/(2*c0); (-q + sqrt(dsc))/(2*c0)];
        r = sort(r(r > 0));
    end
end
bk = [0; r(:); Inf];
plus = zeros(0,2); minus = zeros(0,2);
for j = 1:numel(bk)-1
    if mod(j,2) == 1, plus(end+1,:) = [bk(j) bk(j+1)]; else, minus(end+1,:) = [bk(j) bk(j+1)]; end %#ok<AGROW>
end
end

% ---------- the closed-form v-kernels ----------
function [h0,h1,h2] = ctr_h(z)
% h_k(z) = int_0^Inf u^k exp(-u^2 - 2 z u) du, k = 0, 1, 2, Re z >= 0: erfcx and the forward recurrence for |z| < 6,
% the asymptotic series to its smallest term beyond (ctr_h of DFExact.jl)
h0 = complex(zeros(size(z))); h1 = h0; h2 = h0;
sm = abs(z) < 6;
if any(sm(:))
    zs = z(sm);
    a0 = (sqrt(pi)/2)*faddeeva(1i*zs);                % erfcx(z) = w(iz)
    a1 = (1 - 2*zs.*a0)/2;
    h0(sm) = a0; h1(sm) = a1; h2(sm) = (a0 - 2*zs.*a1)/2;
end
lg = ~sm;
if any(lg(:))
    % (the elements still summing are kept together in compact arrays; pos maps them to the output)
    zl = z(lg); pos = find(lg);
    w = 1./(2*zl); w2 = w.*w;
    a0 = w; a1 = w2; a2 = 2*w2.*w;
    s0 = a0; s1 = a1; s2 = a2;
    j = 0;
    while ~isempty(pos)
        r2 = -(2*j + 3)*(2*j + 4)/(j + 1)*w2;
        stop = abs(r2) >= 1;                         % the smallest term is reached: no update
        if any(stop)
            h0(pos(stop)) = s0(stop); h1(pos(stop)) = s1(stop); h2(pos(stop)) = s2(stop);
            k = ~stop; pos = pos(k); w2 = w2(k); r2 = r2(k);
            a0 = a0(k); a1 = a1(k); a2 = a2(k); s0 = s0(k); s1 = s1(k); s2 = s2(k);
            if isempty(pos), break, end
        end
        r0 = -(2*j + 1)*(2*j + 2)/(j + 1)*w2;
        r1 = -(2*j + 2)*(2*j + 3)/(j + 1)*w2;
        a0 = a0.*r0; a1 = a1.*r1; a2 = a2.*r2;
        s0 = s0 + a0; s1 = s1 + a1; s2 = s2 + a2;
        j = j + 1;
        cv = abs(a0) < 1e-17*abs(s0) & abs(a2) < 1e-17*abs(s2);
        if any(cv)
            h0(pos(cv)) = s0(cv); h1(pos(cv)) = s1(cv); h2(pos(cv)) = s2(cv);
            k = ~cv; pos = pos(k); w2 = w2(k);
            a0 = a0(k); a1 = a1(k); a2 = a2(k); s0 = s0(k); s1 = s1(k); s2 = s2(k);
        end
    end
end
end

function [n,R0,R1,R2] = ctr_rend(a,sa,beta,Ea)
% the remainders R_m(a) = int_a^Inf v^m e^(-t v^2 - s X(v)) dv less the full-line Gaussian when it belongs to a, and
% the flag n (1: the saddle lies beyond a in the Stokes sense, Re z_a < 0); a is a row of endpoints, Ea a matrix (rows
% the nodes, columns the endpoints), all endpoints in one evaluation of ctr_h
z = sa.*a + beta./(2*sa);
n = double(~(real(z) >= 0));
zz = z; ng = n == 1; zz(ng) = -z(ng);
[g0,g1,g2] = ctr_h(zz);
g0(ng) = -g0(ng); g2(ng) = -g2(ng);
i1 = 1./sa; i2 = i1.*i1; i3 = i2.*i1;
R0 = Ea.*i1.*g0;
R1 = Ea.*(a.*i1.*g0 + i2.*g1);
R2 = Ea.*(a.*a.*i1.*g0 + 2*a.*i2.*g1 + i3.*g2);
end

function [M0,M1,M2] = ctr_moments(t,s,c0,q,pcs)
% int over the pieces (a, b) of v^m e^(-t v^2 - s X(v)) dv, m = 0, 1, 2: the endpoint remainders, and the full-line
% Gaussian as often as the flags at the lower ends exceed those at the upper ends
alpha = t + c0*s; beta = q*s; sa = sqrt(alpha);
ends = pcs(:,1).'; fin = isfinite(pcs(:,2)).';
ends = [ends, pcs(fin,2).'];                           % the lower ends, then the finite upper ends
Ea = complex(zeros(numel(s),numel(ends)));
for e = 1:numel(ends)
    a = ends(e); Ea(:,e) = exp(-t*a*a - s*(c0*a*a + q*a + 0.5));
end
[nn,R0,R1,R2] = ctr_rend(ends,sa,beta,Ea);
M0 = complex(zeros(size(s))); M1 = M0; M2 = M0; nf = zeros(size(s));
np = size(pcs,1); iu = np;
for j = 1:np
    M0 = M0 + R0(:,j); M1 = M1 + R1(:,j); M2 = M2 + R2(:,j); nf = nf + nn(:,j);
    if fin(j)
        iu = iu + 1;
        M0 = M0 - R0(:,iu); M1 = M1 - R1(:,iu); M2 = M2 - R2(:,iu); nf = nf - nn(:,iu);
    end
end
m = nf ~= 0;
if any(m(:))
    al = alpha(m); be = beta(m); nm = nf(m);
    G0 = sqrt(pi)./sa(m).*exp(be.*be./(4*al) - s(m)/2); b2 = be./(2*al);
    M0(m) = M0(m) + nm.*G0; M1(m) = M1(m) - nm.*b2.*G0; M2(m) = M2(m) + nm.*(b2.*b2 + 1./(2*al)).*G0;
end
end

function K = ctr_kernel(kind,t,s,c0,q,pcs)
% the v-kernel of the density (kind 1, 2 v^2) or the distribution function (kind 2, X'(v) = q + 2 c0 v)
if isempty(pcs), K = complex(zeros(size(s))); return, end
[M0,M1,M2] = ctr_moments(t,s,c0,q,pcs);
if kind == 1, K = 2*M2; else, K = q*M0 + 2*c0*M1; end
end

function lm = ctr_logM(t,s1,s2,KQ,s,g0,c)
% log Mt(s,t) - log P: -(log(s1 - s) + log(s - s2))/2 + g0^2 Q(s,t) (principal logarithms)
a = s1 - s; b = s - s2;
lm = -(log(a) + log(b))/2 + g0*g0*KQ.*(t + s.*s/2 + c*s)./(a.*b);
end

% ---------- the placement of the U-shapes ----------
function m = ctr_kbound(tr,sig,c0,q,pcs)
% max over v in the pieces of Re(-t v^2) - sig X(v): +Inf when a semi-infinite piece has a growing Gaussian
sz = size(tr + sig);
if isempty(pcs), m = zeros(sz); return, end
m = -Inf(sz); hit = false(sz);
for j = 1:size(pcs,1)
    a = pcs(j,1); b = pcs(j,2);
    A = tr + sig*c0;
    pos = A > 0;
    vv = min(max(-sig*q./(2*A),a),b);
    hp = -A.*vv.*vv - sig*q.*vv - sig/2;
    m(pos) = jmax(m(pos),hp(pos));
    np = ~pos;
    if ~isfinite(b)
        hit = hit | np;
    else
        ha = -A*a*a - sig*q*a - sig/2; hb = -A*b*b - sig*q*b - sig/2;
        m(np) = jmax(jmax(m(np),ha(np)),hb(np));
    end
end
m(hit) = Inf;
end

function sc = ctr_xscale(tr,c0,q,pcs)
% the slowest decay scale along the rays of the v-kernel of the pieces: for each piece 1/min(max|X|, |X(v0)|), v0 the
% vertex of the parabola inside the piece, and |c0|/X'(v_j)^2 at an end that is a root v_j of X; semi-infinite pieces
% cut where e^(-Re(t) v^2) = e^-18; 0 without pieces
sc = zeros(size(tr));
for j = 1:size(pcs,1)
    a = pcs(j,1); b = pcs(j,2);
    if isfinite(b), bb = b + zeros(size(tr)); else, bb = a + 6./sqrt(tr); end
    m = zeros(size(tr));
    for k = 0:32
        v = a + (bb - a)*k/32;
        m = jmax(m,abs(c0*v.*v + q*v + 0.5));
    end
    sp = 0;
    if c0 ~= 0
        v0 = -q/(2*c0);
        if a < v0 && v0 < b, m = min(m,abs(c0*v0*v0 + q*v0 + 0.5)); end
        % at an end that is a root v_j of X the kernel decays like s^(-1/2) (the Gaussian width of c0 (v - v_j)^2) up to
        % |s| ~ |c0|/X'(v_j)^2 before the slope takes over: near a tangency of X (two close roots) on the pieces beside it
        % as well as between them, where the vertex is not inside the piece
        for e = [a b]
            if e > 0 && isfinite(e), sp = max(sp,abs(c0)/max((2*c0*e + q)^2,1e-12)); end
        end
    end
    sc = jmax(sc,jmax(1./max(m,1e-12),sp));
end
end

function m = ctr_bound(D,g0,c,sk,so,sig,side,c0,q,pcs,kfl)
% the largest of log|Mt| + kb(Re s) over the 8 probe points of the U-shape around sk at sig (so the other root); the
% rows are the outer nodes, the columns the abscissas.  kb is the exponent bound of the v-kernel of the pieces
% (ctr_kbound), for the distribution function (kfl 1 on G+, 2 on G-) at least that of the kernels of Psi_t(1/2),
% (1 - e^(-s/2))/s (max(0, -sig/2)) and -1/s (0)
rk = real(sk); y0 = imag(sk); E = zeros(size(sig));
dd = abs(rk - sig);
yo = min(max(imag(so),y0 - dd),y0 + dd);
% the 8 probes side by side (columns k M + 1 .. (k + 1) M), their maximum with NaN propagation (as Julia's max)
s = complex([sig, sig, sig, rk + E, rk + E, rk + side*dd, rk + side*dd, sig], ...
            [y0 + E, y0 + dd, y0 - dd, y0 + dd, y0 - dd, y0 + dd, y0 - dd, yo]);
kb = ctr_kbound(D.tr,real(s),c0,q,pcs);
if kfl == 1, kb = jmax(kb,jmax(0,-real(s)/2)); elseif kfl == 2, kb = jmax(kb,0); end
v = real(ctr_logM(D.t,D.s1,D.s2,D.KQ,s,g0,c)) + kb;
v = reshape(v,[size(sig) 8]);
m = max(v,[],3); m(any(isnan(v),3)) = NaN;
end

function V = ctr_objective(D,g0,c,sk,so,side,c0,q,pcs,kfl,foc,pen,sig)
% the objective B(sig) of the placement: ctr_bound plus the grading cost toward the focus points, weighted by pen
V = ctr_bound(D,g0,c,sk,so,sig,side,c0,q,pcs,kfl);
dd = abs(real(sk) - sig);
for p = 1:size(foc,2)
    fp = foc(:,p);
    dist = abs(complex(sig - real(fp),max(0,abs(imag(fp) - imag(sk)) - dd)));
    V = V + pen*log(1 + dd./max(dist,1e-300));
end
end

function x = ctr_golden(f,a,b,n)
% golden-section minimum of f on [a, b], n steps (elementwise over the outer nodes)
g = (sqrt(5) - 1)/2;
x1 = b - g*(b - a); x2 = a + g*(b - a); f1 = f(x1); f2 = f(x2);
for it = 1:n
    m = f1 < f2; nm = ~m;
    b(m) = x2(m); x2(m) = x1(m); f2(m) = f1(m); x1(m) = b(m) - g*(b(m) - a(m));
    a(nm) = x1(nm); x1(nm) = x2(nm); f1(nm) = f2(nm); x2(nm) = a(nm) + g*(b(nm) - a(nm));
    xn = x2; xn(m) = x1(m);
    fn = f(xn);
    f1(m) = fn(m); f2(nm) = fn(nm);
end
x = x2; s = f1 < f2; x(s) = x1(s);
end

function sig = ctr_place(D,g0,c,sk,so,side,lo,hi,c0,q,pcs,kfl,foc,pen)
% the abscissa of the U-shape around sk (side +1: opening right, -1: left) in the admissible interval [lo, hi]: the
% minimum of the objective, scanned on a geometric grid of 24 distances d = |Re sk - sig| and refined by 12
% golden-section steps between the neighbours of the best
NG = 24;
rk = real(sk);
if side == 1, dmin = rk - hi; dmax = rk - lo; else, dmin = lo - rk; dmax = hi - rk; end
dmin = max(dmin,1e-3*dmax);
% an empty admissible interval (the branch point -t/c0 or the condition Re s < 0 of G- beyond the root; not met in the
% supported region, met for c > 2): NaN, so that the value is not finite and flagged (status 2)
bad = ~(dmax > 0 & dmin < dmax); dmin(bad) = 1; dmax(bad) = 2;
N = numel(rk);
sgs = zeros(N,NG);
for k = 1:NG
    dk = dmin.*(dmax./dmin).^((k - 1)/(NG - 1));
    if side == 1, sgs(:,k) = rk - dk; else, sgs(:,k) = rk + dk; end
end
f = @(s) ctr_objective(D,g0,c,sk,so,side,c0,q,pcs,kfl,foc,pen,s);
V = f(sgs);
best = Inf(N,1); jb = ones(N,1);
for k = 1:NG
    u = V(:,k) < best; best(u) = V(u,k); jb(u) = k;
end
r = (1:N)';
a = sgs(sub2ind([N NG],r,max(jb - 1,1))); b = sgs(sub2ind([N NG],r,min(jb + 1,NG)));
sig = ctr_golden(f,min(a,b),max(a,b),12);
sig(bad) = NaN;
end

function L = ctr_layout(kind,C,idx,q,c0,c,g0,plus,minus,pen)
% the placement at the outer nodes idx: the abscissas sp, sm, the focus points fp, fm and the ray scales rp, rm
MARG = 1e-3;
t = C.t(idx); s1 = C.s1(idx); s2 = C.s2(idx); tr = real(t);
D = struct('t',t,'s1',s1,'s2',s2,'KQ',C.KQ(idx),'tr',tr);
lo = real(s2); hi = real(s1); Lw = hi - lo;
if c0 == 0, sb = complex(zeros(size(t))); else, sb = -t/c0; end
z0 = complex(zeros(size(t)));
semip = c0 > 0 && any(~isfinite(plus(:,2)));
if semip, fp = [z0 sb s2]; else, fp = [z0 s2]; end
lop = lo + MARG*Lw; if semip, lop = max(lop,real(sb) + MARG*Lw); end
% (the exponent bound of the kernel on G+ and on G- below: for the distribution function also that of the kernels of
% Psi_t(1/2), so that a negligible v-kernel does not leave the placement blind to Mt)
sp = ctr_place(D,g0,c,s1,s2,1,lop,hi - MARG*Lw,c0,q,plus,(kind == 2)*1,fp,pen);
semim = c0 < 0 && ~isempty(minus) && any(~isfinite(minus(:,2)));
if semim, fm = [z0 sb s1]; else, fm = [z0 s1]; end
sm = NaN(size(t));
if ~isempty(minus) || kind == 2
    him = hi - MARG*Lw; if semim, him = min(him,real(sb) - MARG*Lw); end
    if kind == 2, him = min(him,-MARG*Lw); end                     % Re s < 0 on G- (the kernel -1/s)
    sm = ctr_place(D,g0,c,s2,s1,-1,lo + MARG*Lw,him,c0,q,minus,(kind == 2)*2,fm,pen);
end
% (the kernel (1 - e^(-s/2))/s of the distribution function on G+ decays on the scale 2 = 1/X(0))
rp = ctr_xscale(tr,c0,q,plus); if kind == 2, rp = max(rp,2); end
L = struct('sp',sp,'sm',sm,'fp',fp,'fm',fm,'rp',rp,'rm',ctr_xscale(tr,c0,q,minus));
end

% ---------- the inner integrals ----------
function [s,wc,wr,ray] = ctr_unodes(sk,sig,side,n,rscale,foc,z,w)
% the nodes and weights of the U-shape around sk with its segment at sig (ctr_upath and ctr_seg of DFExact.jl): the
% segment split and graded toward the focus points, the rays r = rho u/(1 - u) (two scales when rho > 2d); the
% segment weights are i wc (complex), the ray weights wr (real), ray flags the latter
dd = abs(real(sk) - sig); y0 = imag(sk); ylo = y0 - dd; yhi = y0 + dd;
rho = max(dd,rscale);
lev = min(max(imag(foc),ylo),yhi); hs = max(abs(sig - real(foc)),1e-6*dd);
% the focus points in the order of their scales: one whose level lies within its own scale of a level already kept is
% dropped
kl = zeros(1,0); kh = zeros(1,0);
[~,o] = sort(hs);
for k = o(:).'
    if any(abs(lev(k) - kl) <= hs(k)), continue, end
    kl(end+1) = lev(k); kh(end+1) = hs(k); %#ok<AGROW>
end
pts = [ylo Inf];
[~,o] = sort(kl);
tol = 1e-12*max(dd,1);
for k = o(:).'
    if abs(kl(k) - pts(end,1)) <= tol
        pts(end,2) = min(pts(end,2),kh(k));
    else
        pts(end+1,:) = [kl(k) kh(k)]; %#ok<AGROW>
    end
end
if ~(abs(yhi - pts(end,1)) <= tol), pts(end+1,:) = [yhi Inf]; end
S = {}; W = {};
for k = 1:size(pts,1)-1
    A = pts(k,1); hA = pts(k,2); B = pts(k+1,1); hB = pts(k+1,2);
    if B - A <= 0, continue, end
    if ~isfinite(hA) && ~isfinite(hB)
        [S{end+1},W{end+1}] = ctr_seg(sig,A,B,Inf,z,w); %#ok<AGROW>
    elseif ~isfinite(hA)
        [S{end+1},W{end+1}] = ctr_seg(sig,B,A,hB,z,w); %#ok<AGROW>
    elseif ~isfinite(hB)
        [S{end+1},W{end+1}] = ctr_seg(sig,A,B,hA,z,w); %#ok<AGROW>
    else
        M = (A + B)/2;
        [S{end+1},W{end+1}] = ctr_seg(sig,A,M,hA,z,w); %#ok<AGROW>
        [S{end+1},W{end+1}] = ctr_seg(sig,B,M,hB,z,w); %#ok<AGROW>
    end
end
nseg = sum(cellfun(@numel,S));
R = {}; WR = {};
for e = [yhi 1; ylo -1].'
    yc = e(1); dir = e(2);
    cs = complex(sig,yc); if side == 1, sgn = dir; else, sgn = -dir; end
    r0 = 0;
    if rho > 2*dd
        % two scales: graded panels from the corner (first length d, growing by 3) up to rho, the map beyond
        a = 0; h = dd;
        while a < rho*(1 - 1e-14)
            b = min(a + h,rho); if (rho - b) < 0.5*h, b = rho; end
            R{end+1} = cs + side*(a + (b - a)*z); WR{end+1} = sgn*(b - a)*w; %#ok<AGROW>
            a = b; h = h*3;
        end
        r0 = rho;
    end
    r = r0 + rho*z./(1 - z);
    R{end+1} = cs + side*r; WR{end+1} = sgn*rho./((1 - z).*(1 - z)).*w; %#ok<AGROW>
end
s = [vertcat(S{:}); vertcat(R{:})];
wc = [vertcat(W{:}); zeros(numel(s) - nseg,1)];
wr = [zeros(nseg,1); vertcat(WR{:})];
ray = [false(nseg,1); true(numel(s) - nseg,1)];
end

function [s,wc] = ctr_seg(sig,st,en,h0,z,w)
% the nodes of int f(sig + i y) i dy between st and en by Gauss-Legendre panels graded from st (first length h0, Inf:
% one panel, growing by the factor 3); the weights wc of i dy
len = abs(en - st); if en >= st, dir = 1; else, dir = -1; end
a = 0; if isfinite(h0), h = min(h0,len); else, h = len; end
S = {}; W = {};
while a < len*(1 - 1e-14)
    b = min(a + h,len); if (len - b) < 0.5*h, b = len; end
    S{end+1} = complex(sig,st + dir*(a + (b - a)*z)); W{end+1} = (b - a)*w; %#ok<AGROW>
    a = b; h = h*3;
end
s = vertcat(S{:}); wc = vertcat(W{:});
if isempty(s), s = complex(zeros(0,1)); wc = zeros(0,1); end
end

function [I,mx] = ctr_inner_n(kind,C,idx,q,c0,c,g0,plus,minus,n,L)
% the inner function I_f(t) (kind 1) or I_F(t) (kind 2) at the outer nodes idx with n nodes per panel and ray, the
% placement L given: the value and the largest term
[z,w] = glrule(n); z = z(:); w = w(:);
% no admissible abscissa (ctr_place: NaN): the value is NaN, hence flagged
usem = ~isempty(minus) || kind == 2;
nan = isnan(L.sp(:)); if usem, nan = nan | isnan(L.sm(:)); end
I = complex(NaN(numel(idx),1)); mx = zeros(numel(idx),1);
ok = find(~nan);
if isempty(ok), I = reshape(I,size(idx)); mx = reshape(mx,size(idx)); return; end
Lk = L; Lk.sp = L.sp(ok); Lk.rp = L.rp(ok); Lk.fp = L.fp(ok,:);
if usem, Lk.sm = L.sm(ok); Lk.rm = L.rm(ok); Lk.fm = L.fm(ok,:); end
[vp,mp] = ctr_upaths(kind,C,idx(ok),q,c0,c,g0,plus,1,n,Lk.sp,Lk.rp,Lk.fp,z,w);
vm = complex(zeros(size(vp))); mm = zeros(size(vp));
if usem
    [vm,mm] = ctr_upaths(kind,C,idx(ok),q,c0,c,g0,minus,-1,n,Lk.sm,Lk.rm,Lk.fm,z,w);
end
P = C.P(idx(ok));
I(ok) = P(:).*(vp(:) + vm(:)); mx(ok) = abs(P(:)).*jmax(mp(:),mm(:));
I = reshape(I,size(idx)); mx = reshape(mx,size(idx));
end

function [v,mx] = ctr_upaths(kind,C,idx,q,c0,c,g0,pcs,side,n,sig,rs,foc,z,w)
% (1/2 pi i) int over the U-shapes of side (+1: G+ around s1 with the kernel of the pieces plus, -1: G- around s2)
% at the outer nodes idx, all nodes evaluated together; returns the values and the largest terms / (2 pi)
N = numel(idx);
Sc = cell(N,1); Wc = Sc; Wr = Sc; Ry = Sc; cnt = zeros(N,1);
for j = 1:N
    k = idx(j);
    if side == 1, sk = C.s1(k); else, sk = C.s2(k); end
    [Sc{j},Wc{j},Wr{j},Ry{j}] = ctr_unodes(sk,sig(j),side,n,rs(j),foc(j,:),z,w);
    cnt(j) = numel(Sc{j});
end
s = vertcat(Sc{:}); wc = vertcat(Wc{:}); wr = vertcat(Wr{:}); ray = vertcat(Ry{:});
g = repelem((1:N)',cnt);
t = C.t(idx(g)); s1 = C.s1(idx(g)); s2 = C.s2(idx(g)); KQ = C.KQ(idx(g));
f = exp(ctr_logM(t,s1,s2,KQ,s,g0,c)).*ctr_fkern(kind,side,t,s,c0,q,pcs);
vals = complex(zeros(size(s)));
vals(~ray) = f(~ray).*complex(0,wc(~ray));
vals(ray) = f(ray).*wr(ray);
v = complex(zeros(N,1)); mx = zeros(N,1);
e = cumsum(cnt); b = e - cnt + 1;
av = abs(vals);
for j = 1:N
    v(j) = sum(vals(b(j):e(j)));
    a = av(b(j):e(j)); m = max(a); if any(isnan(a)), m = NaN; end
    mx(j) = m;
end
v = v/(2*pi*1i); mx = mx/(2*pi);
end

function K = ctr_fkern(kind,side,t,s,c0,q,pcs)
% the kernel factor of the integrand: the v-kernel of the pieces, plus (1 - e^(-s/2))/s on G+ and -1/s on G- for the
% distribution function
K = ctr_kernel(kind,t,s,c0,q,pcs);
if kind == 2
    if side == 1
        ps = (1 - exp(-s/2))./s;
        sm = abs(s) < 1e-8; ps(sm) = 0.5 - s(sm)/8;
        K = K + ps;
    else
        K = K - 1./s;
    end
end
end

% ---------- helpers ----------
function m = jmax(a,b)
% max with NaN propagation (as Julia's max)
m = max(a,b); m(isnan(a) | isnan(b)) = NaN;
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
