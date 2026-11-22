function CP = dfdet_cpath(Mo,c,opts)
%DFDET_CPATH  The outer contour of the general engine in t, with the roots and log P(t) continued along it.
%
%   CP = DFDET_CPATH(Mo,c,opts) for a model Mo (DFDET_STARTMODEL) at c: the
%   hyperbola t(u) = tau + a (cosh u - 1) + i b sinh u, asymptotic angle phi,
%   y = min(0.5, 0.9 (pi/2 - phi)) the half-width of the strip in u,
%   b = 0.9 min(tau, t_* - tau)/sin(y), a = b cot(phi), step h = hfac 2 pi y/
%   (37 + 37 y tan(phi)), up to |t - tau| = 1e5 (20 c^2 below c = -100),
%   tau = min(1, t_*/2): the nodes t_k = t(k h), t'(u_k), the weights of the
%   trapezoidal rule (h/pi, h/(2 pi) at u = 0), the roots s_k(t_k) (the first
%   Mo.rp the group positive at tau) and log P(t_k), P(t)^-2 = D_V(t)
%   prod_+ mu_k prod_- (-mu_k), so that e^{s/2} M(s,t) = P(t) prod_+ (s_k -
%   s)^-1/2 prod_- (s - s_k)^-1/2, continued along the path on 8 points per
%   step, evaluated only where a step is not safe (a root moves by 0.45 times
%   its distance to the nearest other root, or the increment of the phase of
%   log P differs from the prediction by more than pi/4).  phi: pi/4 by
%   default; with a fixed start and below c = -100 the smallest of the angles
%   45 to 88 degrees at which every term |M_V(t_k) e^{g0^2 cf(t_k)} t'/t_k|
%   h/pi of the outer rule is at most 2.  Where the positive group leaves
%   Re s > 0 the hyperbola is steepened by 5 degrees up to 88; beyond, the
%   vertical line t = tau + i w with Gauss-Legendre panels in w (the first of
%   length min(tau,1)/4, doubling; 24 nodes each), on which the root groups
%   stay on their sides of the gap at the real point.  opts: phi, hfac (1),
%   tau, and for the line nq (24), cuts (the panel boundaries in w) and ts.
%   The paths are cached per model, c, phi and hfac; the classical cases DK,
%   DU, DT at c = 0 with a zero start load from the tables dfdet_paths.mat
%   (DFDET_MAKE_PATHS).
%
%   See also DFDET_STARTMODEL, DFDET_CONTOUR, DFDET_TSTAR.
if nargin < 3 || isempty(opts), opts = struct; end
phi = getf(opts,'phi',[]); hfac = getf(opts,'hfac',1); tau = getf(opts,'tau',[]);
if isfield(opts,'fixedphi')
    % the angle of the rule of the fixed start alone (the law of V where the t-ratio takes the vertical line)
    CP = fixed_phi(Mo,c,opts.ts,tau,hfac,default_rmax(c));
    return
end
if isfield(opts,'cuts') || isfield(opts,'nq')
    % the vertical line with given panels or nodes per panel (line_refine, the variants of the checked computation)
    rmax = default_rmax(c); if ~isempty(getf(opts,'rmax',[])), rmax = opts.rmax; end
    CP = linecpath(Mo,c,opts.ts,tau,rmax,getf(opts,'nq',24),8,getf(opts,'cuts',[]));
    CP.key = sprintf('%s|%.17g|line%d|%d',Mo.key,c,getf(opts,'nq',24),numel(CP.t));
    return
end
if isempty(phi), pk = -1; else, pk = phi; end
key = sprintf('%s|%.17g|%.17g|%.17g',Mo.key,c,pk,hfac);
CP = dfdet_cache('get','cpath',key);
if ~isempty(CP), return, end
if isempty(phi) && hfac == 1 && isempty(tau), CP = table_path(Mo,c); end
if isempty(CP), CP = cpath_(Mo,c,phi,hfac,8,[],tau); end
CP.key = key;
dfdet_cache('put','cpath',key,CP);
end

% ---------------------------------------------------------------------------------------------------------------------
function r = default_rmax(c)
% the end |t - tau| of the outer contour: 1e5, and 20 c^2 below c = -100 (far from the unit root the transform decays
% along the contour only once |t| exceeds c^2)
if c < -100, r = max(1e5,20*c*c); else, r = 1e5; end
end

function [a,b,h,K] = hyperbola(ts,phi,hfac,rmax,tau)
% (at most 1e5 nodes, PATH_MAXNODES of DFDet.jl: the default rules need at most some 5500, 11000 in the finer variant of
% the checked computation; a factor hfac of the step that is not positive or far below one is refused)
y = min(0.5,0.9*(pi/2 - phi));
b = 0.9*min(tau,ts - tau)/sin(y); a = b/tan(phi);
h = hfac*2*pi*y/(37 + 37*tan(phi)*y);
K = asinh(rmax/hypot(a,b))/h;
if ~(h > 0 && K <= 1e5), error('dfdet:paths','the outer contour would need %.3g nodes (hfac = %g, the angle %.3g degrees), more than 1e5: hfac must be a positive factor of order one.',K,hfac,phi*180/pi); end
K = ceil(K);
end

function CP = cpath_(Mo,c,phi,hfac,sub,rmax,tau)
if isempty(rmax), rmax = default_rmax(c); end
ts = dfdet_tstar(Mo,c);
if isempty(tau), tau = min(1,ts/2); end
auto = isempty(phi);
if auto
    if strcmp(Mo.start,'fixed') || c < -100, phi = fixed_phi(Mo,c,ts,tau,hfac,rmax); else, phi = pi/4; end
end
while true
    [a,b,h,K] = hyperbola(ts,phi,hfac,rmax,tau);
    ud = (0:sub*K)'*h/sub;
    tt = complex(tau + a*(cosh(ud) - 1),b*sinh(ud));
    % (with the automatic angle a path on which the positive group leaves Re s > 0 is discarded: its walk stops there)
    [R,lP,FX,minpos,maxneg] = rootpath(Mo,c,tt,sub,auto && phi < 88*pi/180);
    % GLS detrending, c > 0, designs without a proper intercept: far out on a flat hyperbola the positive roots can move
    % into Re s < 0, where the kernels of G+ grow; steepen the hyperbola until the positive group stays right of the
    % origin (on the line Re t = tau it does)
    if ~(auto && minpos <= 0 && phi < 88*pi/180), break, end
    phi = min(phi + 5*pi/180,88*pi/180);
end
if auto && minpos <= 0
    CP = linecpath(Mo,c,ts,tau,rmax,24,8,[]);
    return
end
wt = repmat(h/pi,K + 1,1); wt(1) = wt(1)/2;
k = (0:K)';
CP = struct('c',c,'tau',tau,'tstar',ts,'a',a,'b',b,'h',h,'phi',phi,'t',tt(1:sub:end),'dt',complex(a*sinh(k*h),b*cosh(k*h)), ...
    'R',R,'lP',lP,'minpos',minpos,'maxneg',maxneg,'Mo',Mo,'fx',FX,'wt',wt);
end

function CP = linecpath(Mo,c,ts,tau,rmax,nq,sub,cuts)
% the outer integral on the vertical line t = tau + i w itself, (1/pi) int_0^inf Im{I_F(t) i/t} dw, by Gauss-Legendre
% panels (nq nodes) in w, the first of length min(tau, 1)/4 and doubling up to rmax, the roots continued on sub points
% between nodes
[z,w] = glrule(nq);
if isempty(cuts)
    cuts = 0; hh = min(tau,1)/4;
    while cuts(end) < rmax, cuts(end+1) = cuts(end) + hh; hh = hh*2; end %#ok<AGROW>
end
cuts = cuts(:);
L = diff(cuts);
ws = reshape(cuts(1:end-1)' + L'.*z,[],1); wq = reshape(L'.*w,[],1)/pi;
% the path for the continuation: 0, the nodes, and sub - 1 points between consecutive ones
prev = [0; ws(1:end-1)];
pathw = [0; reshape([prev' + (ws - prev)'.*((1:sub-1)'/sub); ws'],[],1)];
tt = complex(tau,pathw);
[R,lP,FX,minpos,maxneg] = rootpath(Mo,c,tt,sub,false);
R = R(:,2:end); lP = lP(2:end); FX = fxsel(FX,2:size(R,2)+1);
nn = numel(ws);
CP = struct('c',c,'tau',tau,'tstar',ts,'a',0,'b',1,'h',NaN,'phi',pi/2,'t',complex(tau,ws),'dt',complex(zeros(nn,1),ones(nn,1)), ...
    'R',R,'lP',lP,'minpos',minpos,'maxneg',maxneg,'Mo',Mo,'fx',FX,'wt',wq);
end

function F = fxsel(F,idx)
if isempty(F), return, end
F.cf = F.cf(idx); F.H = F.H(:,:,idx); F.w1 = F.w1(:,idx); F.p1 = F.p1(:,idx);
F.w2 = F.w2(:,idx); F.p2 = F.p2(:,idx); F.mt = F.mt(:,idx);
end

function phi = fixed_phi(Mo,c,ts,tau,hfac,rmax)
% the asymptotic angle of the outer hyperbola with a fixed start (and below c = -100): the smallest of the angles at which
% the largest term |M_V(t_k) e^{g0^2 cf(t_k)} t'(u_k)/t_k| h/pi of the outer rule is finite and at most 2 (the steepest
% one if none qualifies)
PHI = [45 50 55 60 65 70 75 78 80 82 84 85 86 87 88];
g2 = Mo.g0^2; fixed = strcmp(Mo.start,'fixed');
for phd = PHI
    phi = phd*pi/180;
    [a,b,h,K] = hyperbola(ts,phi,hfac,rmax,tau);
    ok = true;
    for k = 0:K
        t = complex(tau + a*(cosh(k*h) - 1),b*sinh(k*h)); dt = complex(a*sinh(k*h),b*cosh(k*h));
        if fixed
            [~,ld,~,fx] = dfdet_mgf(Mo,c,t); e = g2*real(fx.cf);
        else
            [~,ld] = dfdet_mgf(Mo,c,t,true); e = 0;
        end
        v = exp(-real(ld)/2 + e)*abs(dt/t)*h/pi;
        if ~(isfinite(v) && v <= 2), ok = false; break, end
    end
    if ok, return, end
end
phi = PHI(end)*pi/180;
end

% ---------------------------------------------------------------------------------------------------------------------
function [R,lP,FX,minpos,maxneg] = rootpath(Mo,c,tt,sub,stopneg)
% the roots and log P along the fine path tt (sub points per step), evaluated at the nodes tt(1:sub:end) and at the
% sub - 1 points between two nodes only when the step between them is not safe for the continuation (rootpath_sparse
% of DFDet.jl); minpos and maxneg over the points evaluated.  stopneg: stop as soon as a root of the positive group has
% Re s <= 0 (minpos is then <= 0 and the path is discarded by the caller)
nt = numel(tt); r = Mo.r; rp = Mo.rp; rm = Mo.rm;
idx = 1:sub:nt; K = numel(idx);
fixed = strcmp(Mo.start,'fixed');
MU = complex(zeros(r,K)); LD = complex(zeros(K,1));
if fixed
    FX = struct('cf',complex(zeros(K,1)),'H',complex(zeros(r,r,K)),'w1',complex(zeros(r,K)),'p1',complex(zeros(r,K)), ...
        'w2',complex(zeros(r,K)),'p2',complex(zeros(r,K)),'mt',complex(zeros(r,K)));
else
    FX = [];
end
% (the transform at the nodes is computed as the walk reaches them, so that a discarded path stops early)
done = false(K,1);
    function node(k)
        if fixed
            [~,ld,mu,fx] = dfdet_mgf(Mo,c,tt(idx(k)));
            FX.cf(k) = fx.cf; FX.H(:,:,k) = fx.H; FX.w1(:,k) = fx.w1; FX.p1(:,k) = fx.p1;
            FX.w2(:,k) = fx.w2; FX.p2(:,k) = fx.p2; FX.mt(:,k) = fx.mt;
        else
            [~,ld,mu] = dfdet_mgf(Mo,c,tt(idx(k)));
        end
        MU(:,k) = mu; LD(k) = ld; done(k) = true;
    end
node(1);
R = complex(zeros(r,K)); lP = complex(zeros(K,1));
s1 = 1./MU(:,1);
if ~(max(abs(imag(s1))) < 1e-8*max(abs(s1))), warning('dfdet:rootpath','rootpath: roots not real at the start.'); end
pos = sort(real(s1(real(s1) > 0))); neg = sort(real(s1(real(s1) < 0)));
if ~(numel(pos) == rp && numel(neg) == rm)
    error('dfdet:rootpath','rootpath: %d positive and %d negative roots at the start, inertia of A %d, %d.',numel(pos),numel(neg),rp,rm);
end
R(:,1) = complex([pos; neg]);
Lr1 = lr(LD(1),R(:,1),rp);
m = round(imag(Lr1)/(2*pi));
if ~(abs(imag(Lr1) - 2*pi*m) < 1e-7), warning('dfdet:rootpath','rootpath: D_V prod mu not positive at the start.'); end
L1 = complex(real(Lr1),0); lP(1) = -L1/2;
Lprev = L1; dprev = NaN;
mnp = min(real(R(1:rp,1))); if rm == 0, mxn = -Inf; else, mxn = max(real(R(rp+1:r,1))); end
for k = 2:K
    if ~done(k), node(k); end
    Rp = R(:,k-1);
    % the predicted increment of the imaginary part: the previous one, scaled by the ratio of the steps in t
    if isnan(dprev), dpred = NaN; else, dpred = dprev*abs(tt(idx(k)) - tt(idx(k-1)))/abs(tt(idx(k-1)) - tt(idx(k-2))); end
    nw = match(1./MU(:,k),Rp);
    Lr = lr(LD(k),nw,rp);
    dp = dpred; if isnan(dp), dp = 0; end
    Lc = Lr + 2i*pi*round((imag(Lprev) + dp - imag(Lr))/(2*pi));
    % safe: every root moved by less than 0.45 times its distance to the nearest other root (a double root counted
    % once), and the increment agrees with the prediction to pi/4
    safe = ~isnan(dpred) && abs(imag(Lc) - imag(Lprev) - dpred) <= pi/4;
    if safe && r > 1
        dd = abs(Rp - Rp.');
        dd(dd <= 1e-6*abs(Rp)) = Inf;                    % (the diagonal, and a double root counted once)
        dmin = min(dd,[],2);
        safe = all(abs(nw - Rp) < 0.45*dmin);
    end
    if ~safe
        % walk through the points between the nodes (the rule of the dense continuation for the multiple of 2 pi)
        Rw = Rp; Lw = Lprev;
        for j = idx(k-1)+1:idx(k)
            if j < idx(k)
                [~,ld,mu] = dfdet_mgf(Mo,c,tt(j));
                [Rw,Lw] = step(Rw,Lw,mu,ld,rp);
                mnp = min(mnp,min(real(Rw(1:rp)))); if rm > 0, mxn = max(mxn,max(real(Rw(rp+1:r)))); end
            else
                [Rw,Lw] = step(Rw,Lw,MU(:,k),LD(k),rp);
            end
        end
        nw = Rw; Lc = Lw;
    end
    R(:,k) = nw; lP(k) = -Lc/2;
    dprev = imag(Lc) - imag(Lprev); Lprev = Lc;
    mnp = min(mnp,min(real(R(1:rp,k)))); if rm > 0, mxn = max(mxn,max(real(R(rp+1:r,k)))); end
    if stopneg && mnp <= 0
        R = R(:,1:k); lP = lP(1:k); break
    end
end
minpos = mnp; maxneg = mxn;
end

function L = lr(ld,rts,rp)
L = ld + sum(log(1./rts(1:rp))) + sum(log(-1./rts(rp+1:end)));
end

function [nw,Lc] = step(Rp,Lp,mu,ld,rp)
nw = match(1./mu,Rp);
Lr = lr(ld,nw,rp);
Lc = Lr + 2i*pi*round((imag(Lp) - imag(Lr))/(2*pi));
end

function new = match(new,prev)
% the assignment of the new roots to the previous labels: the permutation of least total distance (exhaustive for
% r <= 6, in the order of DFDet.jl, the first of least distance), greedy beyond
r = numel(new);
if r == 1, return, end
if r <= 6
    Pm = juliaperms(r);
    d = abs(reshape(new(Pm),size(Pm)) - prev(:).');
    cost = d(:,1); for k = 2:r, cost = cost + d(:,k); end
    [~,ib] = min(cost);
    new = new(Pm(ib,:));
else
    left = 1:r; out = new;
    for k = 1:r
        [~,j] = min(abs(new(left) - prev(k))); out(k) = new(left(j)); left(j) = [];
    end
    new = out;
end
end

function Pm = juliaperms(n)
% the permutations of 1..n in the order of permutations_of of DFDet.jl (n inserted at every position of each permutation
% of 1..n-1), cached
persistent PC
if isempty(PC), PC = cell(1,6); end
if ~isempty(PC{n}), Pm = PC{n}; return, end
if n == 1, Pm = 1; PC{1} = Pm; return, end
Q = juliaperms(n - 1); Pm = zeros(size(Q,1)*n,n); i = 0;
for a = 1:size(Q,1)
    for k = 1:n
        i = i + 1; q = Q(a,:);
        Pm(i,:) = [q(1:k-1) n q(k:end)];
    end
end
PC{n} = Pm;
end

% ---------------------------------------------------------------------------------------------------------------------
function CP = table_path(Mo,c)
% the outer paths of cases DK, DU, DT at c = 0 with a zero start, shipped as tables (dfdet_paths.mat, written by
% DFDET_MAKE_PATHS from the tables of DFDet.jl); [] for any other model, or when the file is missing
CP = [];
if ~(c == 0 && strcmp(Mo.start,'zero') && isnan(Mo.P.cbar)), return, end
if ~isfield(Mo,'key') || isempty(strfind(Mo.key,'|ls|')), return, end   % (the least squares statistic, not the LM one)
D = Mo.P.D;
if ~(strcmp(D.detrend,'onestep') && ~D.loc && isequal(D.br,[0 1])), return, end
p = -2;
for q = -1:1
    T = dfdet_trend(q);
    if isequal(size(D.C{1}),size(T.C{1})) && isequal(D.C{1},T.C{1}), p = q; break, end
end
if p == -2, return, end
persistent TB
if isempty(TB)
    f = fullfile(fileparts(mfilename('fullpath')),'dfdet_paths.mat');
    if ~exist(f,'file'), TB = struct('none',true); return, end
    TB = load(f);
end
nm = {'DK','DU','DT'}; nm = nm{p + 2};
if ~isfield(TB,nm), return, end
T = TB.(nm);
if ~(size(T.R,1) == Mo.r && T.rp == Mo.rp), error('dfdet:table','dfdet_cpath: the table %s has %d roots (%d positive), the design %d (%d).',nm,size(T.R,1),T.rp,Mo.r,Mo.rp); end
CP = struct('c',T.c,'tau',T.tau,'tstar',T.tstar,'a',T.a,'b',T.b,'h',T.h,'phi',T.phi,'t',T.t(:),'dt',T.dt(:),'R',T.R, ...
    'lP',T.lP(:),'minpos',T.minpos,'maxneg',T.maxneg,'Mo',Mo,'fx',[],'wt',T.wt(:));
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
