function varargout = dfdet_contour(mode,CP,arg,opts,c0s)
%DFDET_CONTOUR  The inversion on paths around the roots: the inner integrals over the U-shapes and the outer sum.
%
%   [F,f,nodes,eF,ef] = DFDET_CONTOUR('tratio',CP,xs,opts): the distribution
%   function F(x) = P(tau <= x) and the density f(x) of tau = (U - c0 V)/sqrt(V)
%   at the points xs on the outer contour CP (DFDET_CPATH, which carries the
%   design and its start): the trapezoidal rule in u >= 0 on the hyperbola
%   (on the vertical line its Gauss-Legendre panels), F = sum' w_k Im{I_F(t_k)
%   t'(u_k)/t_k}, f likewise with Im{I_f(t_k) t'(u_k)}, the outer nodes in
%   chunks of 32 until a chunk whose terms all stay below opts.tol (1e-17).
%   At a node t
%      I_F(t) = (1/2 pi i) [ int_{G+} Mt (W+_F + (1 - e^{-a s})/s) ds
%                          + int_{G-} Mt (W-_F - 1/s) ds ],
%      I_f(t) = (1/2 pi i) [ int_{G+} Mt W+_f ds + int_{G-} Mt W-_f ds ],
%   Mt = e^{a s} M(s,t) (a = 1/2 except with a fixed start), W+-_F =
%   int_{V+-} X'(v) e^{-t v^2 - s X(v)} dv, W+-_f = int_{V+-} 2 v^2 e^{-t v^2 -
%   s X(v)} dv, X(v) = c0 v^2 + x v + a, V+ and V- the sets of v > 0 where
%   X > 0 and X < 0, incomplete Gaussian integrals in closed form (the
%   functions h_k from the complex erfcx, FADDEEVA, for |z| < 6 and their
%   asymptotic series beyond).  G+ is the U-shape around the positive group of
%   roots, its segment at alpha min Re s_+, and G- the one around the negative
%   group at alpha max Re s_- (none without a negative root): a vertical
%   segment graded by Gauss-Legendre panels toward the levels of the roots and
%   of the origin, and two horizontal rays graded toward the roots they pass,
%   continued by the map r = R + R' (u/(1 - u))^2.  With a fixed start the
%   abscissas are placed by minimizing the largest Re log Mt at probe points
%   of the U-shape plus a grading cost, the rays pass each root at its own
%   horizontal distance, and every panel is cut to the rate at which log Mt
%   varies.  opts: c0 (0), n and ntail (nodes per panel and on the tails of
%   the rays: 24 and 40, with a fixed start 48 and 64), alpha (0.5), tol, warn
%   (true), vpow (2; 3 replaces the density kernel 2 v^2 by 2 v^3: at x = 0 the
%   derivative of F in c0, the density of the coefficient statistic), filon
%   (true: on the vertical line the parts of the kernels at the root of the
%   parabola, which carry e^{-t v^2}, by Filon's rule on each panel).  eF, ef:
%   estimates of the rounding error (eps times the sums of (1 + |log Mt|)
%   |term| over all nodes), with the warning dfdet:rounding above 1e-4 or at
%   a value that is not finite, and dfdet:contour when the outer contour ends
%   with terms above 1e-10.
%
%   [F,f,nodes,eF,ef] = DFDET_CONTOUR('specs',CP,specs,opts): events made of
%   parabola terms: specs{i} a structure array with the fields c0, q, A, vmax,
%   w, dc0, dq, dA, bnd (DFDET_XSPEC); F(i) = sum over its terms of
%   w P(R <= X(sqrt V), sqrt V <= vmax), R = U + a, X(v) = c0 v^2 + q v + A, and
%   f(i) its derivative in the parameter whose derivatives of the
%   coefficients are (dc0, dq, dA).  The constant term P(R <= A under the tilt)
%   is split at R = 0: (1 - e^{-As})/s on G+ and -1/s on G- for A > 0, -e^{-As}/s
%   on G- for A < 0; with a finite vmax the boundary term
%   -e^{-t vmax^2} P(R <= X(vmax) under the tilt) is added (bnd).
%
%   CP = DFDET_CONTOUR('refine',CP,xs,amax,c0s): the vertical line with its
%   panels in w cut to the oscillation of the outer integrand at |x| <= max|xs|
%   (above the unit root the inner integral oscillates in w at the rate
%   kappa^2 x^2/4 + kappa |A|, kappa = |Im s_k|/w, and with a local null value
%   c0 ~= 0 (c0s) at the rate of the stationary value of its exponent, faster by
%   about 1/|1 - c0/c| in case DK): each panel cut into parts of length at most
%   30/Omega; cached; CP itself when nothing needs cutting.
%
%   See also DFDET_CPATH, DFDET_CDF, DFDET_XSPEC.
if nargin < 4 || isempty(opts), opts = struct; end
switch mode
    case 'tratio'
        [varargout{1:nargout}] = eval_tratio(CP,arg,opts);
    case 'specs'
        [varargout{1:nargout}] = eval_specs(CP,arg,opts);
    case 'refine'
        if nargin < 5, c0s = []; end
        varargout{1} = line_refined(CP,arg,opts,c0s);
    otherwise
        error('dfdet:args','dfdet_contour: unknown mode %s.',mode);
end
end

% =====================================================================================================================
% the two front ends
% =====================================================================================================================
function [F,f,done,eF,ef] = eval_tratio(CP,xs,opts)
Mo = CP.Mo; fixed = strcmp(Mo.start,'fixed'); A = Mo.a;
c0 = getf(opts,'c0',0); vpow = getf(opts,'vpow',2);
[n,ntail] = nodes(opts,fixed);
xs = xs(:).'; nx = numel(xs);
ev = cell(1,nx);
for ix = 1:nx
    [pl,mi] = xpieces(xs(ix),c0,A);
    ev{ix} = struct('c0',c0,'q',xs(ix),'A',A,'vmax',Inf,'w',1,'dc0',double(vpow == 3),'dq',double(vpow == 2),'dA',0,'bnd',true, ...
        'plus',pl,'minus',mi);
end
vr = [];
if getf(opts,'filon',true) && isnan(CP.h), vr = filon_roots_x(CP,xs,c0); end
if isempty(vr) && getf(opts,'warn',true) && isnan(CP.h)
    b2 = 0;
    for ix = 1:nx
        pc = [ev{ix}.plus; ev{ix}.minus];
        fin = pc(isfinite(pc(:,2)),2); if ~isempty(fin), b2 = max(b2,max(fin.^2)); end
    end
    line_edge_warn(CP,b2,c0 < 0);
end
U = struct('mode','tratio','c0s',c0,'A',A,'xq',max([0 xs(isfinite(xs))]),'n',n,'ntail',ntail,'alpha',getf(opts,'alpha',0.5), ...
    'nk',3 + (vpow == 3),'cache',c0 == 0 && ~(fixed && CP.c > 0));
[F,f,done,eF,ef] = csum(CP,ev,U,vr,getf(opts,'tol',1e-17),getf(opts,'warn',true));
end

function [F,f,done,eF,ef] = eval_specs(CP,specs,opts)
Mo = CP.Mo; fixed = strcmp(Mo.start,'fixed');
[n,ntail] = nodes(opts,fixed);
ne = numel(specs); ev = cell(1,ne);
c0s = zeros(1,0); ak = Mo.a; d3 = false; xq = 0;
for i = 1:ne
    S = specs{i};
    for j = 1:numel(S)
        [pl,mi] = xpieces_gen(S(j).q,S(j).c0,S(j).A,S(j).vmax);
        S(j).plus = pl; S(j).minus = mi; %#ok<AGROW>
        ak = max(ak,abs(S(j).A)); d3 = d3 || S(j).dc0 ~= 0; xq = max(xq,S(j).q);
        if S(j).c0 ~= 0 && ~any(c0s == S(j).c0), c0s(end+1) = S(j).c0; end %#ok<AGROW>
    end
    ev{i} = S;
end
vr = [];
if getf(opts,'filon',true) && isnan(CP.h), vr = filon_roots(CP,ev); end
if isempty(vr) && getf(opts,'warn',true) && isnan(CP.h)
    b2 = 0; neg = false;
    for i = 1:ne
        for j = 1:numel(ev{i})
            pc = [ev{i}(j).plus; ev{i}(j).minus];
            for q = 1:size(pc,1)
                b = pc(q,2);
                if isfinite(b) && b < ev{i}(j).vmax, b2 = max(b2,b*b); neg = neg || ev{i}(j).c0 < 0; end
            end
        end
    end
    line_edge_warn(CP,b2,neg);
end
U = struct('mode','specs','c0s',c0s,'A',ak,'xq',xq,'n',n,'ntail',ntail,'alpha',getf(opts,'alpha',0.5),'nk',3 + d3,'cache',false);
[F,f,done,eF,ef] = csum(CP,ev,U,vr,getf(opts,'tol',1e-17),getf(opts,'warn',true));
end

function [n,ntail] = nodes(opts,fixed)
if fixed, n = 48; ntail = 64; else, n = 24; ntail = 40; end
n = getf(opts,'n',n); ntail = getf(opts,'ntail',ntail);
end

% =====================================================================================================================
% the outer sum
% =====================================================================================================================
function [F,f,done,eF,ef] = csum(CP,ev,U,vr,tol,warn)
% the outer trapezoidal rule (or the line's Gauss-Legendre rule) over the nodes of CP, in chunks of 32 until a chunk
% whose terms all stay below tol; with vr (Filon's rule on the line) the parts at the roots v1 = vr(i) separately
K1 = numel(CP.t); nx = numel(ev); CH = 32;
ushape_guard(CP,U,tol);                                % (the guard on the nodes of the U-shapes, before the outer sum)
fil = ~isempty(vr);
GF = zeros(K1,nx); Gf = GF; AF = GF; Af = GF;
if fil, JF = complex(zeros(K1,nx)); Jf = JF; eb = exp(-CP.tau*vr.^2); else, JF = []; Jf = []; end
done = 0; last = Inf;
Uc = [];
if U.cache
    ukey = sprintf('%s|%.17g|%.17g|%d|%d|%.17g',CP.key,U.c0s,U.alpha,U.n,U.ntail,U.A);
    Uc = dfdet_cache('get','ushape',ukey);
    if isempty(Uc), Uc = struct('ch',{cell(1,ceil(K1/CH))},'n',0); end
end
while done < K1
    k1 = done + 1; k2 = min(done + CH,K1); ich = ceil(k2/CH);
    if ~isempty(Uc) && ~isempty(Uc.ch{ich})
        B = Uc.ch{ich};
        [GF(k1:k2,:),Gf(k1:k2,:),AF(k1:k2,:),Af(k1:k2,:),jF,jf] = chunk_terms(CP,k1:k2,B,ev,U,vr);
        if fil, JF(k1:k2,:) = jF; Jf(k1:k2,:) = jf; end
    else
        % the U-shapes of the outer nodes one by one, their terms in batches of at most 2e5 nodes of the U-shapes (where
        % the roots run off, as for designs without a proper intercept on the vertical line, a U-shape can carry 1e5
        % nodes and more)
        kb = k1; Bs = {};
        for k = k1:k2
            Bk = chunk_nodes(CP,k,ev,U);
            if isfield(Bk,'lm')
                % a large U-shape (the roots far out): the terms of the outer node from parts of 2e5 of its nodes, summed
                if ~isempty(Bs)
                    B = catnodes(Bs);
                    [GF(kb:k-1,:),Gf(kb:k-1,:),AF(kb:k-1,:),Af(kb:k-1,:),jF,jf] = chunk_terms(CP,kb:k-1,B,ev,U,vr);
                    if fil, JF(kb:k-1,:) = jF; Jf(kb:k-1,:) = jf; end
                end
                GF(k,:) = 0; Gf(k,:) = 0; AF(k,:) = 0; Af(k,:) = 0; if fil, JF(k,:) = 0; Jf(k,:) = 0; end
                n = numel(Bk.s);
                for a = 1:2e5:n
                    r = (a:min(a + 2e5 - 1,n))'; L = Bk.lm(Bk.s(r));
                    Bp = struct('s',Bk.s(r),'g',ones(numel(r),1),'side',Bk.side(r),'m',exp(L).*Bk.W(r),'E',1 + abs(L));
                    [g1,g2,g3,g4,jF,jf] = chunk_terms(CP,k,Bp,ev,U,vr);
                    GF(k,:) = GF(k,:) + g1; Gf(k,:) = Gf(k,:) + g2; AF(k,:) = AF(k,:) + g3; Af(k,:) = Af(k,:) + g4;
                    if fil, JF(k,:) = JF(k,:) + jF; Jf(k,:) = Jf(k,:) + jf; end
                end
                kb = k + 1; Bs = {}; Uc = []; continue
            end
            Bs{end+1} = Bk; %#ok<AGROW>
            nb = sum(cellfun(@(b) numel(b.s),Bs));
            if nb > 2e5 || k == k2
                B = catnodes(Bs);
                [GF(kb:k,:),Gf(kb:k,:),AF(kb:k,:),Af(kb:k,:),jF,jf] = chunk_terms(CP,kb:k,B,ev,U,vr);
                if fil, JF(kb:k,:) = jF; Jf(kb:k,:) = jf; end
                % (the nodes of a path are kept up to 2e6 of them, some 100 MB, a chunk at a time)
                if ~isempty(Uc) && kb == k1 && k == k2 && Uc.n + numel(B.s) <= 2e6, Uc.ch{ich} = B; Uc.n = Uc.n + numel(B.s); end
                kb = k + 1; Bs = {};
            end
        end
    end
    done = k2;
    mx = max([max(max(abs(GF(k1:k2,:)))) max(max(abs(Gf(k1:k2,:))))]);
    if fil
        mx = max([mx, max(eb.*max(abs(JF(k1:k2,:)),[],1)), max(eb.*max(abs(Jf(k1:k2,:)),[],1))]);
    end
    last = mx*max(CP.wt(k1:k2));
    if mx < tol, break, end
end
if ~isempty(Uc), dfdet_cache('put','ushape',ukey,Uc); end          % (not where a U-shape was large)
w = CP.wt(1:done);
F = (w.'*GF(1:done,:)); f = (w.'*Gf(1:done,:));
if fil
    [nq,wa,Ls] = linepanels(CP);
    Fd = filondata(nq);
    for ix = 1:nx
        b = vr(ix)^2; if ~(b > 0), continue, end         % (the part carries e^{-t b}, b = v1^2)
        for p = 1:numel(wa)
            j1 = (p - 1)*nq + 1; if j1 > done, break, end
            L = Ls(p);
            if b*L < 30
                % (the slower panels keep the Gauss-Legendre rule for the whole integrand)
                kk = j1:min(j1 + nq - 1,done);
                e = exp(-CP.t(kk)*b);
                F(ix) = F(ix) + sum(CP.wt(kk).*imag(e.*JF(kk,ix))); f(ix) = f(ix) + sum(CP.wt(kk).*imag(e.*Jf(kk,ix)));
                continue
            end
            ph = L*exp(-1i*b*wa(p))*eb(ix)/pi;
            F(ix) = F(ix) + imag(ph*filon_int(Fd,JF(j1:j1+nq-1,ix),b*L));
            f(ix) = f(ix) + imag(ph*filon_int(Fd,Jf(j1:j1+nq-1,ix),b*L));
        end
    end
end
eF = eps*(w.'*AF(1:done,:)); ef = eps*(w.'*Af(1:done,:));
if warn && (max([eF ef 0]) > 1e-4 || ~all(isfinite(F)) || ~all(isfinite(f)))
    warning('dfdet:rounding','dfdet: the terms of the sums are large: the rounding error may reach %.2g; the result may be inaccurate (start %s, g0 = %g, c = %g).',max([eF ef]),CP.Mo.start,CP.Mo.g0,CP.c);
end
if warn && last > 1e-10
    warning('dfdet:contour','dfdet: the outer contour ends with terms of size %.2g; the result may be inaccurate (start %s, g0 = %g, c = %g).',last,CP.Mo.start,CP.Mo.g0,CP.c);
end
end

function [GF,Gf,AF,Af,JF,Jf] = chunk_terms(CP,kk,B,ev,U,vr)
% the terms of the outer nodes kk for every event, all nodes of the U-shapes of the chunk together
nk = numel(kk); nx = numel(ev);
GF = zeros(nk,nx); Gf = GF; AF = GF; Af = GF; JF = complex(zeros(nk,nx)); Jf = JF;
t = CP.t(kk); dt = CP.dt(kk); tr = real(t);
s = B.s; g = B.g; sd = B.side; m = B.m; E = B.E;
tn = t(g); trn = tr(g);
sac = containers.Map('KeyType','double','ValueType','any');
for ix = 1:nx
    S = ev{ix};
    tF = complex(zeros(size(s))); tf = tF; jF = tF; jf = tF;
    v1 = 0; if ~isempty(vr), v1 = vr(ix); end
    for j = 1:numel(S)
        T = S(j);
        if isKey(sac,T.c0), sa = sac(T.c0); else, sa = sqrt(tn + T.c0*s); sac(T.c0) = sa; end
        i1 = 1./sa;
        for side = [1 -1]
            sel = sd == side;
            if ~any(sel), continue, end
            if side == 1, pc = T.plus; else, pc = T.minus; end
            ss = s(sel); ts = tn(sel); sas = sa(sel); i1s = i1(sel);
            act = trn(sel).*(pc(:,1).^2).' <= 50;          % the pieces with e^{-Re(t) a^2} not negligible at their lower end
            [M0,M1,M2,M3] = pmoments(ts,ss,sas,i1s,T.c0,T.q,T.A,pc,act,U.nk);
            [pp,pm] = psi_pm(ss,T.A);
            if side == 1, kF = T.q*M0 + 2*T.c0*M1 + pp; else, kF = T.q*M0 + 2*T.c0*M1 + pm; end
            if T.bnd && isfinite(T.vmax)
                vm = T.vmax; [bp,bm] = psi_pm(ss,T.c0*vm*vm + T.q*vm + T.A);
                if side == 1, kF = kF - exp(-ts*vm*vm).*bp; else, kF = kF - exp(-ts*vm*vm).*bm; end
            end
            tF(sel) = tF(sel) + T.w*kF;
            tf(sel) = tf(sel) + T.w*2*(T.dc0*M3 + T.dq*M2 + T.dA*M1);
            if v1 > 0
                % the part of the kernel at the end v1 (the root of X; Filon's rule) without its factor e^{-t v1^2}: the
                % terms of the moments at v1, with the sign of v1 as an upper (-) or lower (+) end of the active pieces
                sg = act*((pc(:,1) == v1) - (pc(:,2) == v1));
                if any(sg ~= 0)
                    q = sg ~= 0; Ev = sg(q).*exp(-ss(q)*(T.c0*v1*v1 + T.q*v1 + T.A));
                    [e0,e1,e2,e3] = rend_part(v1,sas(q),i1s(q),T.q*ss(q),Ev,U.nk);
                    a = find(sel); a = a(q);
                    jF(a) = jF(a) + T.w*(T.q*e0 + 2*T.c0*e1);
                    if strcmp(U.mode,'tratio') && T.dc0 ~= 0
                        jf(a) = jf(a) + 2*e3;
                    else
                        jf(a) = jf(a) + T.w*2*(T.dq*e2 + T.dA*e1);
                    end
                end
            end
        end
    end
    tF = m.*tF; tf = m.*tf;
    IF = accumarray(g,tF,[nk 1]); If = accumarray(g,tf,[nk 1]);
    aF = accumarray(g,E.*abs(tF),[nk 1]); af = accumarray(g,E.*abs(tf),[nk 1]);
    if v1 > 0
        JJF = accumarray(g,m.*jF,[nk 1]); JJf = accumarray(g,m.*jf,[nk 1]);
        ebn = exp(-t*v1*v1);
        IF = IF - ebn.*JJF; If = If - ebn.*JJf;
        JF(:,ix) = JJF.*dt./t; Jf(:,ix) = JJf.*dt;
    end
    GF(:,ix) = imag(IF.*dt./t); Gf(:,ix) = imag(If.*dt);
    AF(:,ix) = aF.*abs(dt./t); Af(:,ix) = af.*abs(dt);
end
end

% =====================================================================================================================
% the U-shapes
% =====================================================================================================================
function B = chunk_nodes(CP,kk,ev,U)
% the nodes s of the U-shapes of the outer nodes kk (g: the outer node within the chunk, side: +1 on G+, -1 on G-), Mt
% times the weights there (m) and the factors 1 + |log Mt| of the error weights (E)
Mo = CP.Mo; fixed = strcmp(Mo.start,'fixed');
nk = numel(kk); Sc = cell(nk,1); Gc = Sc; Dc = Sc; Mc = Sc; Ec = Sc;
for i = 1:nk
    k = kk(i); tr = real(CP.t(k));
    rsp = 0; rsm = 0;
    if strcmp(U.mode,'tratio')
        if U.c0s ~= 0
            for ix = 1:numel(ev)
                pl = active(ev{ix}.plus,tr); mi = active(ev{ix}.minus,tr);
                rsp = max(rsp,xscale(tr,U.c0s,ev{ix}.q,U.A,pl)); rsm = max(rsm,xscale(tr,U.c0s,ev{ix}.q,U.A,mi));
            end
            rsp = min(rsp,1e4); rsm = min(rsm,1e4);
        end
    else
        rAp = 0; rAm = 0;
        for ix = 1:numel(ev)
            for j = 1:numel(ev{ix})
                T = ev{ix}(j);
                % (the kernel psi(s, A) of every term carries e^{-A s}, which decays only on the scale 1/|A|: the rays reach
                % 1/|A| for every term)
                if T.A > 0, rAp = max(rAp,1/T.A); elseif T.A < 0, rAm = max(rAm,-1/T.A); end
                if T.c0 == 0, continue, end
                rsp = max(rsp,xscale(tr,T.c0,T.q,T.A,active(T.plus,tr))); rsm = max(rsm,xscale(tr,T.c0,T.q,T.A,active(T.minus,tr)));
            end
        end
        rsp = max(min(rsp,1e4),min(rAp,1e8)); rsm = max(min(rsm,1e4),min(rAm,1e8));
    end
    lm = logmt_fun(CP,k);
    [Sp,Wp,Sm,Wm] = ushapes(CP,k,lm,U.c0s,rsp,rsm,U.alpha,U.A,U.n,U.ntail,U.xq,fixed);
    s = [Sp; Sm]; W = [Wp; Wm];
    Sc{i} = s; Gc{i} = i*ones(numel(s),1); Dc{i} = [ones(numel(Sp),1); -ones(numel(Sm),1)];
    if numel(kk) == 1 && numel(s) > 2e5
        % (a large U-shape: Mt and the weights of its nodes in parts, by the caller)
        B = struct('s',s,'g',Gc{i},'side',Dc{i},'W',W,'lm',lm); return
    end
    L = lm(s); Mc{i} = exp(L).*W; Ec{i} = 1 + abs(L);
end
B = struct('s',vertcat(Sc{:}),'g',vertcat(Gc{:}),'side',vertcat(Dc{:}),'m',vertcat(Mc{:}),'E',vertcat(Ec{:}));
end

function B = catnodes(Bs)
% the nodes of the U-shapes of consecutive outer nodes (each built alone) as those of one batch
n = cellfun(@(b) numel(b.s),Bs); g = cell(size(Bs));
for i = 1:numel(Bs), g{i} = i*ones(n(i),1); end
B = struct('s',cell2mat(cellfun(@(b) b.s,Bs(:),'UniformOutput',false)),'g',vertcat(g{:}), ...
    'side',cell2mat(cellfun(@(b) b.side,Bs(:),'UniformOutput',false)),'m',cell2mat(cellfun(@(b) b.m,Bs(:),'UniformOutput',false)), ...
    'E',cell2mat(cellfun(@(b) b.E,Bs(:),'UniformOutput',false)));
end

function lm = logmt_fun(CP,k)
% log Mt(s) at the outer node k: the algebraic factors and, with a fixed start, the term g0^2 (cf + E1(s))
rts = CP.R(:,k); rp = CP.Mo.rp; lP = CP.lP(k);
if strcmp(CP.Mo.start,'fixed')
    fx = struct('H',CP.fx.H(:,:,k),'w1',CP.fx.w1(:,k),'p1',CP.fx.p1(:,k),'w2',CP.fx.w2(:,k),'p2',CP.fx.p2(:,k));
    g2 = CP.Mo.g0^2; cf = CP.fx.cf(k);
    lm = @(s) logmtilde(rts,rp,lP,s) + g2*(cf + dfdet_fixe1(fx,s));
else
    lm = @(s) logmtilde(rts,rp,lP,s);
end
end

function L = logmtilde(rts,rp,lP,s)
% log of e^{s/2} M(s,t) = P(t) prod_{k <= rp} (s_k - s)^{-1/2} prod_{k > rp} (s - s_k)^{-1/2} (principal logarithms)
s = s(:); L = lP + zeros(size(s));
for j = 1:rp, L = L - log(rts(j) - s)/2; end
for j = rp+1:numel(rts), L = L - log(s - rts(j))/2; end
end

function ushape_guard(CP,U,tol)
% The guard on the nodes of the U-shapes, before any U-shape is built or any inner integral evaluated (USHAPE_GUARD of
% DFDet.jl, ushape_guard of the C engine): the outer path carries the roots at its nodes, and with them the count that the
% guard of USHAPES compares with 1e6 at every node t_k: the segments of the narrow U-shapes with the zero and the
% stationary start, the wide heights at the abscissas of PLACE with a fixed start (placed only where an upper bound of the
% count exceeds the cap: every abscissa lies in the interval of PLACE, every cap of the segments is at least 20/A).  The
% first node beyond the cap at which the transform has not fallen below the tolerance of the outer sum, |M_V(t_k) t'_k|
% max(1, 1/|t_k|) >= tol, refuses the call with the error of USHAPES.  The nodes where the transform is below tol are left
% to the guard of USHAPES, as the outer sum stops before them unless its terms are not small there (d(r) = 1 + r on the
% vertical line at c = 0 is beyond the cap from |t| = 3e4 on, where |M_V| is 1e-36, and its sum stops at |t| = 1.6e4).
% The design (1 + r, 1{r <= 0.4}) on the vertical line is beyond the cap from |t| = 6300 on, where |M_V| is 1e-15: it took
% seconds and some 0.8 GB (in Julia and C) before USHAPES refused, and is refused here before any inner integral.
Mo = CP.Mo; rp = Mo.rp; fixed = strcmp(Mo.start,'fixed');
A = U.A; n = U.n; alpha = U.alpha; c0s = U.c0s; c = CP.c;
% (the counts of all nodes at once: the gap, the default abscissas, the branch points, the caps, as in USHAPES)
gp = CP.R(1:rp,:); gm = CP.R(rp+1:end,:); K = size(CP.R,2);
hi = min(real(gp),[],1); if Mo.rm > 0, lo = max(real(gm),[],1); else, lo = -Inf(1,K); end
ok = lo < min(hi,0);
mid = hi > 0 & lo < 0;
lo2 = max(lo,hi - 2*max(abs(hi),1));
sp = hi - alpha*(hi - lo2); sm = lo2 + alpha*(min(0,hi) - lo2);
sp(mid) = alpha*hi(mid); sm(mid) = alpha*lo(mid);
tr = real(CP.t(:)).';
if fixed
    % (an upper bound of the count: sigma+ >= L, the lower end of the interval of PLACE, sigma- <= min(hi, 0))
    if Mo.rm > 0, L = lo; else, L = hi - 2*max(abs(hi),1); end
    hmx = max([zeros(1,K); abs(imag(gp)) + real(gp) - L; abs(imag(gm)) + min(hi,0) - real(gm)],[],1);
    cand = find(ok & 2*hmx*A/20*n > 1e6);
else
    for c0 = c0s
        if ~(c0 > 0), continue, end
        lb = -tr/c0;
        sp(~(lb < hi)) = NaN;
        mv = sp <= lb; sp(mv) = lb(mv) + alpha*(hi(mv) - lb(mv));
    end
    hseg = 20/A; hp = hseg*ones(1,K); hm = hp;
    if c < -60, q = sp > 0; hp(q) = 20./min(A,40./sp(q)); q = sm < 0; hm(q) = 20./min(A,40./abs(sm(q))); end
    % (the narrow U-shapes: from min Im - d to max Im + d, d the distance of the segment from the nearest root of the group)
    lp = max(imag(gp),[],1) - min(imag(gp),[],1) + 2*min(real(gp) - sp,[],1);
    nest = lp./hp*n; hsg = lp;
    if Mo.rm > 0
        lm_ = max(imag(gm),[],1) - min(imag(gm),[],1) + 2*min(sm - real(gm),[],1);
        nest = nest + lm_./hm*n; hsg = max(hsg,lm_);
    end
    cand = find(ok & ~isnan(sp) & nest > 1e6);
end
lt = log(tol);
for k = cand
    rts = CP.R(:,k);
    % log |M_V(t_k)| = Re log Mt(0) (with a fixed start and g0^2 Re cf(t_k): E1(0) = 0)
    lmv = real(logmtilde(rts,rp,CP.lP(k),0));
    if fixed, lmv = lmv + Mo.g0^2*real(CP.fx.cf(k)); end
    if ~(lmv + log(abs(CP.dt(k))*max(1,1/abs(CP.t(k)))) >= lt), continue, end
    t = CP.t(k);
    if fixed
        % (the count of USHAPES, with its placed abscissas)
        g1 = rts(1:rp); g2 = rts(rp+1:end);
        [s1,s2] = ufixed_place(CP,logmt_fun(CP,k),c0s,g1,g2,hi(k),lo(k),sp(k),sm(k),tr(k),alpha,A,U.xq);
        s1 = ubranch(s1,hi(k),tr(k),c0s,alpha);
        if isnan(s1), continue, end
        [h1,h2] = ucaps(c,s1,s2,A);
        [ns,hx] = ufixed_count(rts,s1,s2,h1,h2,n);
        if ns > 1e6, error('dfdet:paths',ushape_msg_fixed(),num2str(t),round(ns),hx); end
        continue
    end
    error('dfdet:paths',ushape_msg(),Mo.start,num2str(t),nest(k),hsg(k));
end
end

function m = ushape_msg()
m = 'the U-shapes of the %s start at t = %s would need some %.3g nodes (segments to %.3g long): a root runs off (more than 1e6 nodes).';
end
function m = ushape_msg_fixed()
m = 'the U-shapes of the fixed start at t = %s would need some %d nodes (heights to %.3g): the roots run off (more than 1e6 nodes).';
end

function sp = ubranch(sp,hi,tr,c0s,alpha)
% sigma+ moved right of the branch point -Re(t)/c0 of every c0 > 0 (NaN where a branch point is not left of hi)
for c0 = c0s
    if ~(c0 > 0), continue, end
    lb = -tr/c0;
    if ~(lb < hi), sp = NaN; return, end
    if sp <= lb, sp = lb + alpha*(hi - lb); end
end
end

function [hp,hm] = ucaps(c,sp,sm,A)
% the caps of the panels of the segments of G+ and G-: 20/A, and below c = -60 20/min(A, 40/|sigma|) (there the segments
% lie far from the origin, where the kernels are damped: only rates up to 40/|sigma| matter)
hseg = 20/A;
if c < -60 && sp > 0, hp = 20/min(A,40/sp); else, hp = hseg; end
if c < -60 && sm < 0, hm = 20/min(A,40/abs(sm)); else, hm = hseg; end
end

function [nest,hmx] = ufixed_count(rts,sp,sm,hp,hm,n)
% the count of the wide U-shapes of the fixed start at the abscissas sp, sm (heights to hmx)
hmx = 0;
for z = rts.'
    if real(z) > 0, sg = sp; else, sg = sm; end
    hmx = max(hmx,abs(imag(z)) + abs(real(z) - sg));
end
nest = 2*hmx/min(hp,hm)*n;
end

function [sp,sm] = ufixed_place(CP,lm,c0s,gp,gm,hi,lo,sp,sm,tr,alpha,A,xq)
% the abscissas of the U-shapes of the fixed start: where the largest |Mt| on the U-shape is smallest; above the unit
% root the objective also counts the growth of the v-kernels of the upper tail (x up to xq > 0) on G+ left of the origin,
% and with c0 > 0 the placement is searched right of the branch point -Re(t)/c0
c = CP.c; Mo = CP.Mo;
if c > 0, xg = xq; else, xg = 0; end
lbp = -Inf;
if c > 0, for c0 = c0s, if c0 > 0, lbp = max(lbp,-tr/c0); end, end, end
if isfinite(lbp), sp0 = max(sp,lbp + alpha*(hi - lbp)); else, sp0 = sp; end
c0p = 0; if c > 0 && any(c0s > 0), c0p = max(c0s(c0s > 0)); end
if Mo.rm > 0
    sp = place(lm,gp,gm,1,max(lo,lbp),hi,A,sp0,xg,tr,c0p);
    sm = place(lm,gm,gp,-1,lo,min(hi,0),A,sm,0,tr,0);
else
    sp = place(lm,gp,gm,1,max(hi - 2*max(abs(hi),1),lbp),hi,A,sp0,xg,tr,c0p);
end
end

function [Sp,Wp,Sm,Wm] = ushapes(CP,k,lm,c0s,rsp,rsm,alpha,A,n,ntail,xq,fixed)
% (the gap, the abscissas, the branch points and the caps are written out here, node by node, as USHAPE_GUARD computes
% them for all nodes at once: helper functions called at every node cost up to 12% of a call)
Mo = CP.Mo; t = CP.t(k); tr = real(t); c = CP.c;
rts = CP.R(:,k); rp = Mo.rp;
gp = rts(1:rp); gm = rts(rp+1:end);
hi = min(real(gp)); if Mo.rm > 0, lo = max(real(gm)); else, lo = -Inf; end
if ~(lo < min(hi,0)), error('dfdet:paths','the root groups are not separated by a line Re s = sigma < 0 at t = %s (lo = %g, hi = %g).',num2str(t),lo,hi); end
if hi > 0 && lo < 0
    sp = alpha*hi; sm = alpha*lo;
else
    lo2 = max(lo,hi - 2*max(abs(hi),1));
    sp = hi - alpha*(hi - lo2); sm = lo2 + alpha*(min(0,hi) - lo2);
end
exp_p = complex(zeros(0,1)); exp_m = exp_p;
if fixed, [sp,sm] = ufixed_place(CP,lm,c0s,gp,gm,hi,lo,sp,sm,tr,alpha,A,xq); end
for c0 = c0s
    if c0 == 0, continue, end
    sb = -t/c0;
    if c0 > 0
        lb = -tr/c0;
        if ~(lb < hi), error('dfdet:paths','no abscissa for G+ right of the branch point -t/c0 at t = %s.',num2str(t)); end
        if sp <= lb, sp = lb + alpha*(hi - lb); end
        exp_p(end+1,1) = sb; %#ok<AGROW>
    else
        exp_m(end+1,1) = sb; %#ok<AGROW>
    end
end
hseg = 20/A;
% (below c = -60 the segments lie far from the origin, where the kernels are damped: only rates up to 40/|sigma| matter;
% until November 2026 below c = -100, DFDet.jl FAR_CMIN)
if c < -60 && sp > 0, hp = 20/min(A,40/sp); else, hp = hseg; end
if c < -60 && sm < 0, hm = 20/min(A,40/abs(sm)); else, hm = hseg; end
% (the guard on the size of the U-shapes, which USHAPE_GUARD checks before the outer sum at the nodes where the transform
% is not negligible; here at the other nodes, which the sum reaches where its terms are not small.  The wide heights of
% the fixed start follow a root that runs off)
if fixed
    [nest,hmx] = ufixed_count(rts,sp,sm,hp,hm,n);
    if nest > 1e6, error('dfdet:paths',ushape_msg_fixed(),num2str(t),round(nest),hmx); end
    lmf = lm;
else
    % (the zero and the stationary start: the segments of the narrow U-shapes, at most hp or hm per panel with n nodes, a
    % lower bound of their nodes, within 1% of the count beyond 1e5)
    nest = 0; hsg = 0;
    [ylo,yhi] = ulevels(gp,sp,1,false); nest = nest + (yhi - ylo)/hp*n; hsg = max(hsg,yhi - ylo);
    if Mo.rm > 0, [ylo,yhi] = ulevels(gm,sm,-1,false); nest = nest + (yhi - ylo)/hm*n; hsg = max(hsg,yhi - ylo); end
    if nest > 1e6, error('dfdet:paths',ushape_msg(),Mo.start,num2str(t),nest,hsg); end
    lmf = [];
end
[Sp,Wp] = ushape(gp,gm,sp,1,n,ntail,exp_p,rsp,hp,lmf,A);
if Mo.rm > 0
    [Sm,Wm] = ushape(gm,gp,sm,-1,n,ntail,exp_m,rsm,hm,lmf,A);
else
    Sm = complex(zeros(0,1)); Wm = Sm;
end
end

function [S,W] = ushape(grp,oth,sig,side,n,ntail,extra,rs,hmax,lm,akern)
% the nodes S and weights W of (1/2 pi i) int over the U-shape around the group grp: side +1 the segment Re s = sig with
% the rays to +infinity, -1 the mirror image, traversed from below to above
[z,w] = glrule(n); [zt,wt] = glrule(ntail);
d = min(side*(real(grp) - sig));
if ~(d > 0), error('dfdet:paths','ushape: segment not separated from its roots.'); end
[ylo,yhi] = ulevels(grp,sig,side,~isempty(lm));
foc = [imag(grp(:)), max(abs(real(grp(:)) - sig),1e-6*d); imag(oth(:)), max(abs(real(oth(:)) - sig),1e-6*d); 0, max(abs(sig),1e-6*d); ...
       imag(extra(:)), max(abs(real(extra(:)) - sig),1e-6*d)];
% the panels of the segment (each at most hmax long) and of the rays, as straight pieces from za to zb
seg = graded_panels(ylo,yhi,foc);
pz = complex(zeros(0,2)); kind = zeros(0,1);         % kind 0: segment; 1, 2: ray at yhi, ylo
for p = 1:size(seg,1)
    a0 = seg(p,1); b0 = seg(p,2);
    if isfinite(hmax), mm = max(1,ceil((b0 - a0)/hmax)); else, mm = 1; end
    jj = (0:mm-1)';
    a1 = a0 + (b0 - a0)*jj/mm; b1 = a0 + (b0 - a0)*(jj + 1)/mm;
    pz = [pz; complex(sig,a1), complex(sig,b1)]; kind = [kind; zeros(mm,1)]; %#ok<AGROW>
end
Rr = zeros(1,2); yv = [yhi ylo];
for e = 1:2
    yc = yv(e);
    rf = [side*(real(grp(:)) - sig), max(abs(imag(grp(:)) - yc),1e-6*d)];
    R = max(max(rf(:,1)) + max(2*d,1),rs); Rr(e) = R;
    rp = graded_panels(0,R,rf);
    pz = [pz; complex(sig + side*rp(:,1),yc), complex(sig + side*rp(:,2),yc)]; kind = [kind; e*ones(size(rp,1),1)]; %#ok<AGROW>
end
ms = nsplit(lm,pz(:,1),pz(:,2),n,akern,side);
% the nodes of the panels, each cut into ms equal parts
np = size(pz,1); cnt = ms*n;
S = complex(zeros(sum(cnt) + 2*ntail,1)); W = S; pos = 0;
for p = 1:np
    za = pz(p,1); zb = pz(p,2); m = ms(p);
    L = (zb - za)/m;
    nodes_ = za + L*((0:m-1) + z);                       % n x m
    if kind(p) == 0
        wts = repmat(1i*imag(L)*w,1,m);
    else
        if kind(p) == 1, dir = 1; else, dir = -1; end
        if side == 1, sgn = dir; else, sgn = -dir; end
        wts = repmat(sgn*abs(real(L))*w,1,m);
    end
    S(pos+1:pos+cnt(p)) = nodes_(:); W(pos+1:pos+cnt(p)) = wts(:); pos = pos + cnt(p);
end
for e = 1:2
    yc = yv(e); dir = 3 - 2*e; if side == 1, sgn = dir; else, sgn = -dir; end
    R = Rr(e); Rt = max(R,2);
    yy = zt./(1 - zt);
    r = R + Rt*yy.*yy; dr = 2*Rt*yy./(1 - zt).^2;
    S(pos+1:pos+ntail) = complex(sig + side*r,yc); W(pos+1:pos+ntail) = sgn*dr.*wt; pos = pos + ntail;
end
W = W/(2i*pi);
end

function [ylo,yhi] = ulevels(grp,sig,side,wide)
% the heights of the rays of the U-shape: narrow (zero and stationary start) min Im - d and max Im + d, d the horizontal
% distance of the segment from the nearest root; wide (the fixed start, whose roots are essential singularities)
% min_k (Im s_k - d_k) and max_k (Im s_k + d_k), d_k = |Re s_k - sig|
if wide
    ylo = min(imag(grp) - abs(real(grp) - sig)); yhi = max(imag(grp) + abs(real(grp) - sig));
    return
end
d = min(side*(real(grp) - sig));
ylo = min(imag(grp)) - d; yhi = max(imag(grp)) + d;
end

function out = graded_panels(A,B,foc)
% panels [lo hi] of [A, B] graded toward the focus points (position, scale): a focus point whose position lies within its
% own scale of one with a smaller scale is dropped; between two kept points the interval is split at the midpoint.
% Points within 1e-12 (B - A) of each other or of an end are merged with it, a tolerance relative to the length of the
% interval: far above the unit root the segments of the U-shapes are tiny (the positive root lies within about
% 2c e^{-2c} of the branch point -t/c0 of the coefficient statistic, 1e-16 at c = 20), and with the tolerance
% 1e-12 max(B - A, 1) used before a segment shorter than 1e-12 got no panel at all (from c = 16 on the coefficient
% statistic was off by up to 0.2)
[~,o] = sort(foc(:,2));
kl = zeros(0,1); kh = kl;
for k = o(:).'
    pos = min(max(foc(k,1),A),B); h = foc(k,2);
    if any(abs(pos - kl) <= h), continue, end
    kl(end+1,1) = pos; kh(end+1,1) = h; %#ok<AGROW>
end
tol = 1e-12*(B - A);
pts = [A Inf];
[~,o] = sort(kl);
for k = o(:).'
    if abs(kl(k) - pts(end,1)) <= tol
        pts(end,2) = min(pts(end,2),kh(k));
    else
        pts(end+1,:) = [kl(k) kh(k)]; %#ok<AGROW>
    end
end
if abs(B - pts(end,1)) <= tol, pts(end,1) = B; else, pts(end+1,:) = [B Inf]; end
out = zeros(0,2);
for k = 1:size(pts,1)-1
    a = pts(k,1); ha = pts(k,2); b = pts(k+1,1); hb = pts(k+1,2);
    if b - a <= 0, continue, end
    if ~isfinite(ha) && ~isfinite(hb)
        out(end+1,:) = [a b]; %#ok<AGROW>
    elseif ~isfinite(ha)
        out = seg_panels(out,b,a,hb);
    elseif ~isfinite(hb)
        out = seg_panels(out,a,b,ha);
    else
        m = (a + b)/2;
        out = seg_panels(out,a,m,ha); out = seg_panels(out,b,m,hb);
    end
end
end

function out = seg_panels(out,st,en,h0)
% panels of [st, en] (either order) graded from st: first length h0 (Inf: one panel), doubling
len = abs(en - st); if en >= st, dir = 1; else, dir = -1; end
if len == 0, return, end
a = 0; if isfinite(h0), h = min(h0,len); else, h = len; end
while a < len*(1 - 1e-14)
    b = min(a + h,len); if (len - b) < 0.5*h, b = len; end
    u1 = st + dir*a; u2 = st + dir*b;
    out(end+1,:) = [min(u1,u2) max(u1,u2)]; %#ok<AGROW>
    a = b; h = h*2;
end
end

function m = nsplit(lm,za,zb,n,A,side)
% the number of equal parts into which each panel from za to zb is cut so that the integrand is resolved, from log Mt at
% the ends and the midpoint: the rate omega = max(|L(mid) - L(za)|, |L(zb) - L(mid)|)/(|zb - za|/2) + A, parts short
% enough that omega times their length is at most n; a panel on which the integrand is negligible (max Re L + log|zb - za|
% + max(0, -A Re s) below -40) is not cut; no fixed start (lm empty): 1; at most 64 parts
m = ones(size(za));
if isempty(lm), return, end
len = abs(zb - za); zm = (za + zb)/2;
L = lm([za; zm; zb]); np = numel(za);
La = L(1:np); Lm = L(np+1:2*np); Lb = L(2*np+1:end);
if side == 1, kb = @(z) max(0,-A*real(z)); else, kb = @(z) zeros(size(z)); end
big = max([real(La) + kb(za), real(Lm) + kb(zm), real(Lb) + kb(zb)],[],2) + log(len);
om = max(abs(Lm - La),abs(Lb - Lm))./(len/2) + A;
mm = min(max(ceil(min(om.*len/n,1e6)),1),64);
ok = len > 0 & isfinite(La) & isfinite(Lm) & isfinite(Lb) & ~(big < -40);
m(ok) = mm(ok);
end

function sig = place(lm,grp,oth,side,lo,hi,A,sig0,xq,tr,c0p)
% the abscissa of the U-shape around the group grp in the admissible interval (lo, hi) with a fixed start: the largest
% Re log Mt + kb at the probe points of the U-shape plus the grading cost 0.5 sum log(1 + d/dist) toward the origin and
% the nearest root of the other group, minimized on a grid of distances geometric from both ends (12 points each),
% refined by 12 golden-section steps; sig0 kept unless the minimum is lower by more than 2
if side == 1, edge = min(real(grp)); else, edge = max(real(grp)); end
foc = 0;
if ~isempty(oth), if side == 1, foc(2) = max(real(oth)); else, foc(2) = min(real(oth)); end, end
% (above the unit root with c0 > 0 the branch point -Re(t)/c0 of the v-kernels, the left end of the admissible interval of
% G+, is a focus of the grading cost too: the kernel of the density of the coefficient statistic has a double pole there,
% (t + c0 s)^{-2}, which the objective did not see, and with a fixed start at c from 1 to 5 the density was off by up to
% 0.12 (case DU, c = 5, g0 = 0.5, at the 90% point) without a warning; as place of DFDet.jl)
if side == 1 && c0p > 0, foc(end+1) = -tr/c0p; end
f = @(sg) objective(lm,grp,oth,side,A,xq,tr,c0p,foc,sg);
if side == 1, D = edge - lo; else, D = hi - edge; end
m = 1e-3*D; nh = 12;
g1 = m*(D/2/m).^((0:nh-1)'/(nh - 1));
ds = unique([g1; D - g1]);
vals = f(edge - side*ds);
[~,kb] = min(vals); ng = numel(ds);
a = ds(max(kb - 1,1)); b = ds(min(kb + 1,ng));
gr = (sqrt(5) - 1)/2;
x1 = b - gr*(b - a); x2 = a + gr*(b - a); f1 = f(edge - side*x1); f2 = f(edge - side*x2);
for it = 1:12
    if f1 < f2
        b = x2; x2 = x1; f2 = f1; x1 = b - gr*(b - a); f1 = f(edge - side*x1);
    else
        a = x1; x1 = x2; f1 = f2; x2 = a + gr*(b - a); f2 = f(edge - side*x2);
    end
end
if f1 < f2, sb = edge - side*x1; else, sb = edge - side*x2; end
if min(f1,f2) < f(sig0) - 2, sig = sb; else, sig = sig0; end
end

function V = objective(lm,grp,oth,side,A,xq,tr,c0p,foc,sg)
% the objective of the placement at the abscissas sg (a column): one row of probe points per abscissa
sg = sg(:); ns = numel(sg); grp = grp(:).'; oth = oth(:).';
d = min(side*(real(grp) - sg),[],2);
ylo = min(imag(grp) - abs(real(grp) - sg),[],2); yhi = max(imag(grp) + abs(real(grp) - sg),[],2);
far = max(side*(real(grp) - sg),[],2) + d;
ng = numel(grp);
Pr = [complex(sg,ylo), complex(sg,yhi), complex(sg,min(max(0,ylo),yhi)), ...
      complex(sg + 0*grp,imag(grp) + 0*sg), complex(real(grp) + 0*sg,yhi + 0*grp), complex(real(grp) + 0*sg,ylo + 0*grp), ...
      complex(sg + 0*oth,min(max(imag(oth) + 0*sg,ylo),yhi)), ...
      complex(sg + side*far,ylo), complex(sg + side*1e3*max(far,1),ylo), complex(sg + side*far,yhi), complex(sg + side*1e3*max(far,1),yhi)];
L = lm(Pr(:)); L = reshape(real(L),size(Pr));
if side == 1
    % (the kernels (1 - e^{-A s})/s and e^{-s X} grow to the left of the origin; with xq > 0 the v-kernels of the upper
    % tail by e^{Re(s)^2 xq^2/(4 Re(t + c0 s))}, without bound toward the branch point -Re(t)/c0)
    rs = real(Pr);
    kbv = max(0,-A*rs);
    if xq > 0
        dd = tr + c0p*rs; kx = rs.^2*xq^2./(4*dd); kx(dd <= 0) = Inf;
        ng_ = rs < 0; kbv(ng_) = kbv(ng_) + kx(ng_);
    end
else
    kbv = zeros(size(Pr));
end
v = L + kbv; v(~isfinite(v)) = Inf;
V = max(v,[],2);
for fq = foc
    V = V + 0.5*log(1 + d./max(abs(sg - fq),1e-300));
end
if ng == 0, V(:) = Inf; end
end

% =====================================================================================================================
% the v-kernels in closed form
% =====================================================================================================================
function [M0,M1,M2,M3] = pmoments(t,s,sa,i1,c0,q,A,pcs,act,nk)
% int over the active pieces (a, b) of v^m e^{-t v^2 - s X(v)} dv, m = 0..nk-1, X(v) = c0 v^2 + q v + A: the endpoint
% remainders and the full-line Gaussian as often as the flags at the lower ends exceed those at the upper ends
ns = numel(s);
M0 = complex(zeros(ns,1)); M1 = M0; M2 = M0; M3 = M0; nf = zeros(ns,1);
if isempty(pcs) || ns == 0, return, end
beta = q*s;
for j = 1:size(pcs,1)
    sel = act(:,j);
    if ~any(sel), continue, end
    for e = 1:2
        v = pcs(j,e);
        if ~isfinite(v), continue, end
        ss = s(sel); sas = sa(sel); i1s = i1(sel);
        Ev = exp(-t(sel)*v*v - ss*(c0*v*v + q*v + A));
        z = sas*v + beta(sel).*i1s/2;
        [n,R0,R1,R2,R3] = rend(v,z,i1s,Ev,nk);
        if e == 1
            M0(sel) = M0(sel) + R0; M1(sel) = M1(sel) + R1; M2(sel) = M2(sel) + R2; M3(sel) = M3(sel) + R3; nf(sel) = nf(sel) + n;
        else
            M0(sel) = M0(sel) - R0; M1(sel) = M1(sel) - R1; M2(sel) = M2(sel) - R2; M3(sel) = M3(sel) - R3; nf(sel) = nf(sel) - n;
        end
    end
end
m = nf ~= 0;
if any(m)
    al = sa(m).^2; be = beta(m); nm = nf(m);
    G0 = sqrt(pi)*i1(m).*exp(be.*be./(4*al) - s(m)*A); b2 = be./(2*al);
    M0(m) = M0(m) + nm.*G0; M1(m) = M1(m) - nm.*b2.*G0; M2(m) = M2(m) + nm.*(b2.*b2 + 1./(2*al)).*G0;
    if nk == 4, M3(m) = M3(m) - nm.*(b2.^3 + 3*b2./(2*al)).*G0; end
end
end

function [n,R0,R1,R2,R3] = rend(a,z,i1,Ea,nk)
% the terms of the end a: Ea times the incomplete Gaussian integrals from a, reflected where Re z < 0 (flag n)
n = double(~(real(z) >= 0));
zz = z; ng = n == 1; zz(ng) = -z(ng);
% (where the factor Ea underflows to zero the terms vanish: the h-functions only where it does not)
g0 = complex(zeros(size(z))); g1 = g0; g2 = g0; g3 = g0;
nz = Ea ~= 0;
if all(nz)
    [g0,g1,g2,g3] = hfun(zz,nk);
elseif any(nz)
    [g0(nz),g1(nz),g2(nz),g3(nz)] = hfun(zz(nz),nk);
end
g0(ng) = -g0(ng); g2(ng) = -g2(ng);
i2 = i1.*i1; i3 = i2.*i1;
R0 = Ea.*i1.*g0;
R1 = Ea.*(a*i1.*g0 + i2.*g1);
R2 = Ea.*(a*a*i1.*g0 + 2*a*i2.*g1 + i3.*g2);
if nk == 4
    R3 = Ea.*(a^3*i1.*g0 + 3*a*a*i2.*g1 + 3*a*i3.*g2 + i3.*i1.*g3);
else
    R3 = complex(zeros(size(z)));
end
end

function [e0,e1,e2,e3] = rend_part(v,sa,i1,beta,E,nk)
% the terms of the end v with the factor E = e^{-s X(v)} in place of e^{-t v^2 - s X(v)}: at a root v of X the part of the
% kernel that Filon's rule integrates against e^{-t v^2}
z = sa*v + beta.*i1/2;
[~,e0,e1,e2,e3] = rend(v,z,i1,E,nk);
end

function [h0,h1,h2,h3] = hfun(z,nk)
% h_k(z) = int_0^Inf u^k exp(-u^2 - 2 z u) du, k = 0..nk-1 (nk = 3 or 4), Re z >= 0: erfcx and the forward recurrence for
% |z| < 6, the asymptotic series to its smallest term beyond (ctr_h and ctr_h3 of DFDet.jl)
h0 = complex(zeros(size(z))); h1 = h0; h2 = h0; h3 = h0;
sm = real(z).^2 + imag(z).^2 < 36;
if any(sm)
    zs = z(sm);
    a0 = (sqrt(pi)/2)*faddeeva(1i*zs);                % erfcx(z) = w(iz)
    a1 = (1 - 2*zs.*a0)/2; a2 = (a0 - 2*zs.*a1)/2;
    h0(sm) = a0; h1(sm) = a1; h2(sm) = a2;
    if nk == 4, h3(sm) = a1 - zs.*a2; end
end
lg = ~sm;
if any(lg)
    % the asymptotic series, in bands of |z|: in each band the number of terms at which the series at the smallest |z|
    % of the band has converged to 1e-17 or reached its smallest term (ctr_h stops each element there; the elements of
    % larger |z| converge earlier and their further terms are smaller still, before their own smallest term)
    pos = find(lg); zl = z(lg); az = abs(zl);
    edges = [6 7 8.5 10 12 15 20 30 50 100 Inf];
    for b = 1:numel(edges)-1
        sel = az >= edges(b) & az < edges(b+1);
        if ~any(sel), continue, end
        J = nterms(min(az(sel)),nk);
        zb = zl(sel); w = 1./(2*zb); w2 = w.*w;
        a0 = w; a1 = w2; a2 = 2*w2.*w; a3 = 6*w2.*w2;
        s0 = a0; s1 = a1; s2 = a2; s3 = a3;
        for j = 0:J-1
            a0 = a0.*(-(2*j + 1)*(2*j + 2)/(j + 1)*w2); a1 = a1.*(-(2*j + 2)*(2*j + 3)/(j + 1)*w2); a2 = a2.*(-(2*j + 3)*(2*j + 4)/(j + 1)*w2);
            s0 = s0 + a0; s1 = s1 + a1; s2 = s2 + a2;
            if nk == 4, a3 = a3.*(-(2*j + 4)*(2*j + 5)/(j + 1)*w2); s3 = s3 + a3; end
        end
        q = pos(sel); h0(q) = s0; h1(q) = s1; h2(q) = s2; h3(q) = s3;
    end
end
end

function J = nterms(rho,nk)
% the number of terms of the asymptotic series of h_k at |z| = rho: until the terms of h_0 and of the last h_k fall below
% 1e-17 of the first, or the smallest term is reached (the next ratio at least one)
v = 1/(4*rho*rho); p0 = 1; pl = 1; J = 0;
while true
    if nk == 4, rl = (2*J + 4)*(2*J + 5)/(J + 1)*v; else, rl = (2*J + 3)*(2*J + 4)/(J + 1)*v; end
    if rl >= 1, return, end
    p0 = p0*(2*J + 1)*(2*J + 2)/(J + 1)*v; pl = pl*rl;
    J = J + 1;
    if p0 < 1e-17 && pl < 1e-17, return, end
end
end

function [pp,pm] = psi_pm(s,A)
% the constant terms of the kernels on G+ and G-: the split of P(R <= A under the tilt) at R = 0
if A > 0
    pp = -expm1(-A*s)./s; sm = real(s).^2 + imag(s).^2 < 1e-16; pp(sm) = A - A*A*s(sm)/2;
    pm = -1./s;
elseif A < 0
    pp = complex(zeros(size(s))); pm = -exp(-A*s)./s;
else
    pp = complex(zeros(size(s))); pm = -1./s;
end
end

function [plus,minus] = xpieces(x,c0,A)
% the pieces of {v > 0} with X(v) = c0 v^2 + x v + A > 0 (plus) and < 0 (minus), split at the positive roots of X (a
% double root is no split)
if c0 == 0
    if ~(x < 0), plus = [0 Inf]; minus = zeros(0,2); return, end
    v = -A/x; plus = [0 v]; minus = [v Inf]; return
end
dsc = x*x - 4*c0*A;
if dsc <= 0, r = zeros(0,1); else, r = sort([(-x - sqrt(dsc))/(2*c0); (-x + sqrt(dsc))/(2*c0)]); r = r(r > 0); end
bk = [0; r; Inf];
plus = zeros(0,2); minus = zeros(0,2);
for j = 1:numel(bk)-1
    if mod(j,2) == 1, plus(end+1,:) = [bk(j) bk(j+1)]; else, minus(end+1,:) = [bk(j) bk(j+1)]; end %#ok<AGROW>
end
end

function [plus,minus] = xpieces_gen(q,c0,A,vmax)
% the pieces of (0, vmax) on which X(v) = c0 v^2 + q v + A is positive (plus) and negative (minus), split at its positive
% roots (computed without cancellation; a double root is no split), for any sign of A
r = zeros(0,1);
if c0 == 0
    if q ~= 0, v = -A/q; if v > 0, r = v; end, end
else
    dsc = q*q - 4*c0*A;
    if dsc > 0
        if q >= 0, sg = 1; else, sg = -1; end
        qq = -(q + sg*sqrt(dsc))/2;
        for y = [qq/c0, A/qq]
            if y > 0 && isfinite(y), r(end+1,1) = y; end %#ok<AGROW>
        end
        r = sort(r);
    end
end
bk = [0; r(r < vmax); vmax];
plus = zeros(0,2); minus = zeros(0,2);
for j = 1:numel(bk)-1
    a = bk(j); b = bk(j+1);
    if ~(b > a), continue, end
    if isfinite(b), m = (a + b)/2; else, m = a + 1; end
    X = c0*m*m + q*m + A;
    if X > 0, plus(end+1,:) = [a b]; elseif X < 0, minus(end+1,:) = [a b]; end %#ok<AGROW>
end
end

function sc = xscale(tr,c0,q,A,pcs)
% the slowest decay scale of the v-kernel along the rays: per piece 1/min(max|X|, |X(vertex)|) and |c0|/X'(v_j)^2 at an
% end v_j that is a root of X; semi-infinite pieces cut where e^{-Re(t) v^2} = e^{-18}
sc = 0;
for j = 1:size(pcs,1)
    a = pcs(j,1); b = pcs(j,2);
    if isfinite(b), bb = b; else, bb = a + 6/sqrt(tr); end
    v = a + (bb - a)*(0:32)/32;
    m = max(abs(c0*v.*v + q*v + A));
    sp = 0;
    if c0 ~= 0
        v0 = -q/(2*c0);
        if a < v0 && v0 < b, m = min(m,abs(c0*v0*v0 + q*v0 + A)); end
        for e = [a b]
            if e > 0 && isfinite(e), sp = max(sp,abs(c0)/max((2*c0*e + q)^2,1e-12)); end
        end
    end
    sc = max([sc, 1/max(m,1e-12), sp]);
end
end

function p = active(pcs,tr)
p = pcs(tr*pcs(:,1).^2 <= 50,:);
end

% =====================================================================================================================
% the vertical line: Filon's rule, the refinement, the warning
% =====================================================================================================================
function vr = filon_roots_x(CP,xs,c0)
% the positive root of X(v) = c0 v^2 + x v + a where it has one (0 without), [] where an X has two or the path is not a
% line of Gauss-Legendre panels
vr = [];
if c0 == 0 || isempty(linepanels(CP)), return, end
v = zeros(size(xs)); found = false;
for i = 1:numel(xs)
    if ~isfinite(xs(i)), continue, end
    [pl,mi] = xpieces(xs(i),c0,CP.Mo.a);
    pc = [pl; mi]; rr = unique(pc(isfinite(pc(:,2)) & pc(:,2) > 0,2));
    if isempty(rr), continue, end
    if numel(rr) > 1, return, end
    v(i) = rr; found = true;
end
if found, vr = v; end
end

function vr = filon_roots(CP,ev)
% the events whose parts at a root of the parabola go to Filon's rule: one term with q = 0 (the linear form), c0 ~= 0,
% no derivative in c0, vmax = Inf: its positive root of X(v) = c0 v^2 + A (0 without one); [] otherwise
vr = [];
if isempty(linepanels(CP)), return, end
v = zeros(1,numel(ev)); found = false;
for i = 1:numel(ev)
    S = ev{i}; rr = zeros(0,1);
    for j = 1:numel(S)
        pc = [S(j).plus; S(j).minus];
        rr = [rr; pc(isfinite(pc(:,2)) & pc(:,2) > 0,2)]; %#ok<AGROW>
    end
    if isempty(rr), continue, end
    if ~(numel(S) == 1 && S(1).q == 0 && S(1).c0 ~= 0 && S(1).dc0 == 0 && ~isfinite(S(1).vmax) && numel(unique(rr)) == 1), return, end
    v(i) = rr(1); found = true;
end
if found, vr = v; end
end

function [nq,wa,L] = linepanels(CP)
% the Gauss-Legendre panels [wa(p), wa(p) + L(p)] of a vertical line, each with nq nodes (24, 32, 16 or 48), from the
% nodes themselves; [] if the nodes are not such panels
nq = []; wa = []; L = [];
if ~isnan(CP.h), return, end
w = imag(CP.t); K = numel(w);
for q = [24 32 16 48]
    np = floor(K/q);
    if ~(np > 0 && np*q == K), continue, end
    z = glrule(q);
    W = reshape(w,q,np);
    Lq = (W(q,:) - W(1,:))/(z(q) - z(1)); wq = W(1,:) - Lq*z(1);
    if all(all(abs(W - (wq + Lq.*z)) <= 1e-9*(abs(wq) + Lq)))
        nq = q; wa = wq(:); L = Lq(:); return
    end
end
end

function Fd = filondata(nq)
% the data of Filon's rule on a panel of nq Gauss-Legendre nodes z on [0, 1]: the barycentric weights lam and the
% derivatives of orders 0..nq-1 of the interpolating polynomial at x = 0 (D0) and x = 1 (D1), as rows applied to its
% values: the Taylor coefficients of the Lagrange basis polynomials lam_q prod_{j ~= q} (x - z_j) at 0 and at 1, whose
% products of linear factors with roots of one sign have coefficients of one sign (no cancellation)
persistent FC
if isempty(FC), FC = containers.Map('KeyType','double','ValueType','any'); end
if isKey(FC,nq), Fd = FC(nq); return, end
z = glrule(nq); z = z(:);
lam = zeros(nq,1); D0 = zeros(nq); D1 = zeros(nq);
fk = factorial(0:nq-1)';
for q = 1:nq
    o = z([1:q-1, q+1:nq]);
    lam(q) = 1/prod(z(q) - o);
    c0 = 1; c1 = 1;                                   % coefficients in ascending powers: prod (x - o_j), prod (y - (o_j - 1))
    for j = 1:nq-1
        c0 = [-o(j)*c0; 0] + [0; c0]; c1 = [(1 - o(j))*c1; 0] + [0; c1];
    end
    D0(:,q) = fk.*lam(q).*c0; D1(:,q) = fk.*lam(q).*c1;
end
Fd = struct('z',z,'lam',lam,'D0',D0,'D1',D1);
FC(nq) = Fd;
end

function acc = filon_int(Fd,g,om)
% int_0^1 p(x) e^{-i om x} dx, p the polynomial through the values g at the nodes: for om <= 2000 a composite
% Gauss-Legendre rule (64 nodes on each of ceil(om/20) parts, p by the barycentric formula); beyond, the expansion by
% parts, sum_k (p^(k)(0) - p^(k)(1) e^{-i om})/(i om)^(k+1)
g = g(:); nq = numel(Fd.z);
if om <= 2000
    m = max(1,ceil(om/20)); [zr,wr] = glrule(64);
    x = reshape(((0:m-1) + zr)/m,[],1); wx = repmat(wr/m,m,1);
    dd = x - Fd.z.';                                   % numel(x) x nq
    cc = Fd.lam.'./dd;
    p = (cc*g)./sum(cc,2);
    [ih,jh] = find(dd == 0); p(ih) = g(jh);
    acc = sum(wx.*p.*exp(-1i*om*x));
    return
end
io = 1/complex(0,om); e = exp(-1i*om); acc = 0; fac = io;
d0 = Fd.D0*g; d1 = Fd.D1*g;
for k = 1:nq
    acc = acc + fac*(d0(k) - d1(k)*e); fac = fac*io;
end
end

function line_edge_warn(CP,b2,c0neg)
% the warning of the vertical line for a parabola with a positive root b (b^2 = b2) whose parts Filon's rule does not take
% (ADFmax, two roots): the kernels of the pieces that end at b carry e^{-t b^2}, which oscillates in w at the rate b^2,
% damped only by e^{-tau b^2}
if ~(b2 > 0), return, end
if isnan(CP.Mo.P.cbar), gl = ''; else, gl = ', GLS detrending'; end
if c0neg && CP.tau*b2 < 40
    warning('dfdet:line',['dfdet: on the vertical line Re t = %.3g (c = %g%s: no hyperbola keeps the positive roots right of the ' ...
        'origin) with c0 < 0 the parabola has a root b with b^2 = %.3g, whose part of the outer integrand the line does not resolve: ' ...
        'the result may be off by 1e-5 to 1e-3.'],CP.tau,CP.c,gl,b2);
    return
end
[nq,~,~] = linepanels(CP);
if isempty(nq), nq = 24; end
np = floor(numel(CP.t)/nq);
w = imag(CP.t); est = 0; worst = 0;
for k = 1:np
    j1 = (k - 1)*nq + 1; j2 = k*nq;
    if k == 1, L = 2*w(j2); else, L = (w(j2) - w(j1))*1.02; end
    if ~(b2*L > 30), continue, end
    lmv = max(real(logmtilde(CP.R(:,j1),CP.Mo.rp,CP.lP(j1),0)),real(logmtilde(CP.R(:,j2),CP.Mo.rp,CP.lP(j2),0)));
    est = est + exp(lmv - CP.tau*b2)*L/(pi*max(abs(CP.t(j1)),CP.tau));
    worst = max(worst,b2*L/30);
end
if est > 1e-6
    warning('dfdet:line',['dfdet: on the vertical line Re t = %.3g (c = %g%s: no hyperbola keeps the positive roots right of the ' ...
        'origin) the parabola has a root b with b^2 = %.3g, whose part of the outer integrand (of size %.2g) oscillates up to %.2g ' ...
        'times faster than its panels resolve: the result may be off by more than 1e-5.'],CP.tau,CP.c,gl,b2,est,worst);
end
end

function CPr = line_refined(CP,xs,amax,c0s)
% the line cut to the oscillation of the outer integrand (line_refine), cached per path, the bound of |x| and the values
% |A| beyond the shift a, both rounded up to powers of two, and the local null values c0s ~= 0 of the kernels with a
% slope; CP itself unless it is a vertical line
CPr = CP;
if ~isnan(CP.h), return, end
if nargin < 3 || isempty(amax) || isstruct(amax), amax = CP.Mo.a; end
if nargin < 4, c0s = []; end
fin = abs(xs(isfinite(xs)));
if isempty(fin), xm = 0; else, xm = max(fin); end
if xm > 0, xb = 2^ceil(log2(xm)); else, xb = 0; end
av = amax(isfinite(amax) & amax > CP.Mo.a);
ab = unique(2.^ceil(log2(av(:)')));
if ~(xb > 0 || ~isempty(ab)), return, end
c0s = c0s(:)'; c0s = reshape(unique(c0s(isfinite(c0s) & c0s > 0)),1,[]);   % (only c0 > 0 refines further: LINE_REFINE)
if xb == 0, c0s = []; end
key = sprintf('%s|%.17g|%s',CP.key,xb,sprintf('%.17g,',ab));
if ~isempty(c0s), key = [key '|c0' sprintf('%.17g,',c0s)]; end
hit = dfdet_cache('get','line',key);
if ~isempty(hit), CPr = hit; return, end
CPr = line_refine(CP,xb,[CP.Mo.a ab],c0s);
dfdet_cache('put','line',key,CPr);
end

function CPr = line_refine(CP,xmax,amax,c0s)
% each panel of the line cut into parts of length at most 30/Omega, Omega the largest rate kappa^2 x^2/4 + kappa |A| over
% the roots at the panel's ends and over |x| <= xmax with an amplitude above e^{-40}; the roots continued along the new
% nodes; CP itself when no panel needs cutting.  A local null value c0 ~= 0 (c0s): X(v) = c0 v^2 + x v + a, the exponent
% -alpha v^2 - s_k x v - s_k a with alpha = t + c0 s_k has the stationary value phi_s = (s_k x)^2/(4 alpha) - s_k a, which
% changes with w at the rate Re dphi_s/dt (s_k' = Im s_k/w the drift of the root along the line); for c0 near c in case
% DK alpha ~ t (1 - c0/c) and the rate is that of c0 = 0 times 1/|1 - c0/c| (at c = 10, c0 = 9 the panels cut to the
% rate of c0 = 0 left the distribution function off by 1.1e-3 at x = 1000 and 3.8e-3 at x = 3000); the largest |x| <=
% xmax at which the term is not negligible from log|M_V| + Re phi_s; the cut is to the larger of this rate and that of
% c0 = 0 (line_refine of DFDet.jl).  Only c0 > 0: with c0 < 0, alpha ~ t (1 + |c0|/c) and the rate is below that of
% c0 = 0, which the cut covers (for c0 < 0 the stationary values at the roots of a fixed start gave rates far above the
% oscillation: 1e6 nodes in case DU at c = 10, g0 = 1, c0 = -2, x = 132)
if nargin < 4, c0s = []; end
CPr = CP;
if ~(isnan(CP.h) && (xmax > 0 || max(amax) > CP.Mo.a)), return, end
[nq,wa,Ls] = linepanels(CP);
if isempty(nq), return, end
% only a line with the default panels (the first of length min(tau, 1)/4, doubling) is refined
np = numel(wa); cuts = zeros(np + 1,1); hh = min(CP.tau,1)/4;
for k = 1:np, cuts(k+1) = cuts(k) + hh; hh = hh*2; end
if ~(all(abs(wa - cuts(1:np)) <= 1e-9*cuts(2:end)) && all(abs(Ls - diff(cuts)) <= 1e-9*cuts(2:end))), return, end
w = imag(CP.t); tau = CP.tau;
newcuts = 0; changed = false;
for k = 1:np
    a = cuts(k); b = cuts(k+1); Om = 0;
    for j = [(k - 1)*nq + 1, k*nq]
        wj = w(j); if ~(wj > 0), continue, end
        lmv = real(logmtilde(CP.R(:,j),CP.Mo.rp,CP.lP(j),0));
        for r = CP.R(:,j).'
            kap = abs(imag(r))/wj; ar = abs(real(r));
            for am = amax
                num = 40 + lmv - ar*am; den = tau*kap^2/4 + ar*kap/2;
                if ~(num > 0), continue, end
                if den > 0, xu = min(xmax,sqrt(num/den)); else, xu = xmax; end
                Om = max(Om,kap^2*xu^2/4 + kap*am);
            end
        end
        if xmax > 0
            for c0 = reshape(c0s,1,[])          % (a row: an empty c0s of any shape gives no pass)
                for r = CP.R(:,j).'
                    rp = imag(r)/wj; al = CP.t(j) + c0*r;
                    for am = amax
                        L0 = 40 + lmv - real(r)*am;
                        if ~(L0 > 0), continue, end
                        Q = real(r*r/(4*al));
                        if Q < 0, xu = min(xmax,sqrt(L0/(-Q))); else, xu = xmax; end
                        Om = max(Om,abs(real(xu^2*(2*r*rp*al - r*r*(1 + c0*rp))/(4*al^2))) + abs(rp)*am);
                    end
                end
            end
        end
    end
    m = max(1,ceil(min((b - a)*Om/30,1e6)));
    if m > 1, changed = true; end
    newcuts = [newcuts; a + (b - a)*(1:m)'/m]; %#ok<AGROW>
end
if ~changed, return, end
nn = (numel(newcuts) - 1)*nq;
if nn > 400000, error('dfdet:line','line_refine: the line would need %d nodes for |x| <= %g at c = %g, more than 400000 (the oscillation of the outer integrand grows with x).',nn,xmax,CP.c); end
if numel(newcuts) > 4000, warning('dfdet:line','line_refine: the line needs %d panels for |x| <= %g at c = %g (the oscillation of the outer integrand grows with x); this is slow.',numel(newcuts) - 1,xmax,CP.c); end
CPr = dfdet_cpath(CP.Mo,CP.c,struct('cuts',newcuts,'nq',nq,'ts',CP.tstar,'tau',tau));
CPr.key = sprintf('%s|refined%.17g',CP.key,xmax);
if ~isempty(c0s), CPr.key = [CPr.key '|c0' sprintf('%.17g,',c0s)]; end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
