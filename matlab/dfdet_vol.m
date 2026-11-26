function varargout = dfdet_vol(op,varargin)
%DFDET_VOL  Volatility and persistence regimes in the general engine: the model, the tilt, the transform, the checks.
%
%   VM = DFDET_VOL('model',D,V,cr): the transform model of the design D with
%   the regimes of the volatility path V (VOLPATH) at the c of the regimes cr
%   (cached): the design refined at the regime boundaries (its prepared
%   design P and its zero-start model Mi, the inner model of the paths), np
%   pieces of lengths L with omega^2 (om2) and c on each, the vector ell =
%   (J(1), J at the jumps of d, z = int g J, J at the regime boundaries that
%   are not jumps of d) of dimension Ne, bpos(j) the position in ell of J at
%   the boundary between pieces j and j + 1 (0: not in ell), the padded
%   matrices A, B, Lx of the design and the eigenvalues Dv; the fields start
%   ('zero'), g0, r, rp, rm and key that the outer contour (DFDET_CPATH) reads.
%   DFDET_MGF, DFDET_SIGMA and DFDET_TSTAR take VM in place of an engine model.
%   [th,gam,ex] = DFDET_VOL('tilt',VM,t): the theta_j at t, the diagonal tilt
%   Gamma and the exponent sum_j (c_j + theta_j) L_j.
%   d = DFDET_VOL('dval',VM,s,t): D(s,t) = exp{sum_j (c_j + theta_j) L_j}
%   det{I - 2 (s A + t B + Gamma) Sigma}, M(s,t) = e^{-s/2} D(s,t)^{-1/2}.
%   [VM,c1] = DFDET_VOL('args',fn,D,c,opts,start,g0,gls): the checks of the
%   arguments with opts.vol (vol_args of DFDet_vol.jl) and the model; VM = []
%   for one regime, which the engine itself computes at c1.
%   K = DFDET_VOL('kernel',VM,r,u): the covariance of the path with the
%   regimes, zero start (the Nystrom estimate of t_*).
%
%   See also VOLPATH, DFDET_CDF.
switch op
    case 'model', varargout = {vol_model(varargin{:})};
    case 'tilt', [a,b,c] = vol_tilt(varargin{:}); varargout = {a,b,c};
    case 'dval', varargout = {vol_dval(varargin{:})};
    case 'args', [a,b] = vol_args(varargin{:}); varargout = {a,b};
    case 'kernel', varargout = {vol_kernel(varargin{:})};
    otherwise, error('dfdet:args','dfdet_vol: unknown operation %s.',op);
end
end

function VM = vol_model(D,V,cr)
cr = double(cr(:).');
key = sprintf('%s|vol%s|%s|%s',dfdet_key(D),sprintf('%.17g,',V.lam),sprintf('%.17g,',V.om),sprintf('%.17g,',cr));
VM = dfdet_cache('get','volmodel',key);
if ~isempty(VM), return, end
% the design refined at the regime boundaries (vol_refine): a point within 1e-12 of a piece boundary is that boundary; a
% piece keeps the coefficients of the piece that contains its midpoint, shifted to its left end when they are local
br = D.br;
for x = V.lam, if ~any(abs(br - x) < 1e-12), br(end+1) = x; end, end %#ok<AGROW>
br = sort(br); np = numel(br) - 1;
C = cell(1,np);
for j = 1:np
    mid = (br(j) + br(j+1))/2; jo = 1;
    while jo < numel(D.br) - 1 && mid > D.br(jo+1), jo = jo + 1; end
    if D.loc, C{j} = taylorc(D.C{jo},br(j) - D.br(jo)); else, C{j} = D.C{jo}; end
end
Dr = dfdet_determ(br,C,D.name,D.loc,D.detrend);
Mi = dfdet_startmodel(Dr,0,'zero',0,[]);
P = Mi.P;
np = numel(P.br) - 1; L = diff(P.br);
reg = zeros(1,np);
for j = 1:np
    mid = (P.br(j) + P.br(j+1))/2; i = 1;
    while i < numel(V.om) && mid > V.lam(i+1), i = i + 1; end
    reg(j) = i;
end
om2 = V.om(reg).^2; c = cr(reg);
bpos = zeros(1,max(np - 1,0)); Ne = P.N;
for j = 1:np-1
    b = P.br(j+1);
    i = find(abs(P.jl - b) < 1e-12,1);
    if ~isempty(i), bpos(j) = 1 + i;
    elseif om2(j) ~= om2(j+1) || c(j) ~= c(j+1), Ne = Ne + 1; bpos(j) = Ne; end
end
A = zeros(Ne); B = zeros(Ne); Lx = zeros(Ne,P.r);
A(1:P.N,1:P.N) = P.A; B(1:P.N,1:P.N) = P.B; Lx(1:P.N,:) = P.L;
VM = struct('vm',true,'P',P,'Mi',Mi,'lam',V.lam,'om',V.om,'cr',cr,'np',np,'L',L,'om2',om2,'c',c,'Ne',Ne,'bpos',bpos, ...
    'A',A,'B',B,'Lx',Lx,'Dv',P.Dv(:),'r',P.r,'rp',P.rp,'rm',P.rm,'start','zero','g0',0,'a',0.5,'key',key);
dfdet_cache('put','volmodel',key,VM);
end

function [th,gam,ex] = vol_tilt(VM,t)
th = sqrt(complex(VM.c.^2 - 2*t*VM.om2));
g = (VM.c + th)./(2*VM.om2);
gam = complex(zeros(VM.Ne,1));
gam(1) = gam(1) + g(VM.np);
for j = 1:VM.np-1
    pos = VM.bpos(j);
    if pos > 0, gam(pos) = gam(pos) + g(j) - g(j+1); end
end
ex = sum((VM.c + th).*VM.L);
end

function d = vol_dval(VM,s,t)
[th,gam,ex] = vol_tilt(VM,complex(t));
S = dfdet_sigma(VM,th);
M = eye(VM.Ne) - 2*((s*VM.A + t*VM.B + diag(gam))*S);
d = exp(ex)*det(M);
end

function K = vol_kernel(VM,r,u)
% e^{C(r) + C(u)} int_0^{min(r,u)} omega^2 e^{-2C(v)} dv, C(r) = int_0^r c, piece by piece in closed form, each term with
% its exponents combined so that nothing overflows far from the unit root (vol_kernel of DFDet_vol.jl); r a column, u a row
R = r + 0*u; U = u + 0*r;
br = VM.P.br; np = VM.np;
Cb = [0 cumsum(VM.c.*VM.L)];
E = cat_(R) + cat_(U); m = min(R,U);
K = zeros(size(R));
for j = 1:np
    a = br(j); sel = a < m;
    if ~any(sel(:)), break, end
    b = min(br(j+1),m(sel)); L = b - a; cj = VM.c(j); Ca = Cb(j); Cbb = Ca + cj*L; Es = E(sel);
    z = 2*cj*L; sm = abs(z) < 1e-3;
    v = zeros(size(L));
    zs = z(sm); e0 = 1 - zs/2 + zs.*zs/6 - zs.^3/24;
    v(sm) = exp(Es(sm) - 2*Ca).*L(sm).*e0;
    v(~sm) = (exp(Es(~sm) - 2*Cbb(~sm)) - exp(Es(~sm) - 2*Ca))/(-2*cj);
    K(sel) = K(sel) + VM.om2(j)*v;
end
    function C = cat_(x)
        C = zeros(size(x));
        for q = 1:numel(x)
            jj = 1; while jj < np && x(q) > br(jj+1), jj = jj + 1; end
            C(q) = Cb(jj) + VM.c(jj)*(x(q) - br(jj));
        end
    end
end

function [VM,c1] = vol_args(fn,D,c,opts,start,g0,gls)
% the checks of the arguments of the volatility extension (vol_args of DFDet_vol.jl) and its model
V = opts.vol;
if ~(isstruct(V) && isfield(V,'lam') && isfield(V,'om')), error('dfdet:args','%s: vol must be a volatility path, volpath(breaks, omega)',fn); end
k = numel(V.om);
if isscalar(c), cr = repmat(double(c),1,k); else, cr = double(c(:).'); end
if numel(cr) ~= k, error('dfdet:args','%s: c must be a number or a vector with one value per regime of the volatility path (%d), not %d',fn,k,numel(cr)); end
if ~all(isfinite(cr)), error('dfdet:args','%s: c must be finite',fn); end
if ~(strcmp(start,'zero') || (strcmp(start,'fixed') && g0 == 0))
    error('dfdet:args','%s: a volatility path (vol) is implemented for the zero start only (J(0) = 0); a fixed or a stationary start is not supported with it',fn);
end
if ~isempty(gls), error('dfdet:args','%s: a volatility path (vol) is implemented for least squares detrending only; GLS detrending (gls) is not supported with it',fn); end
if ~strcmp(D.detrend,'onestep'), error('dfdet:args','%s: a volatility path (vol) is implemented for the one-step regression only; the statistic that detrends first (detrend_first) is not supported with it',fn); end
m = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), m = lower(char(opts.method)); end
if ~any(strcmp(m,{'auto','paths'}))
    error('dfdet:args','%s: with a volatility path (vol) the method is the inversion on paths: method must be ''auto'' or ''paths'' (the line integral is not implemented for the volatility transform)',fn);
end
if isscalar(c), cs = sprintf('%g',c); else, cs = ['[' strjoin(arrayfun(@(v) sprintf('%g',v),cr,'UniformOutput',false),', ') ']']; end
if ~all(cr <= 0), error('dfdet:args','%s: with a volatility path (vol) every c must be at most 0 (above the unit root the extension is not implemented); c = %s',fn,cs); end
if ~all(cr >= -480), error('dfdet:args','%s: with a volatility path (vol) every c must be at least -480 (the range of the general engine); c = %s',fn,cs); end
c1 = cr(1);
if k == 1, VM = []; return, end
VM = vol_model(D,V,cr);
end

function a = taylorc(C,s)
% the Taylor coefficients at s of the columns of C (repeated synthetic division)
a = C; D = size(C,1) - 1;
for k = 0:D-1
    for j = D-1:-1:k
        a(j+1,:) = a(j+1,:) + s*a(j+2,:);
    end
end
end
