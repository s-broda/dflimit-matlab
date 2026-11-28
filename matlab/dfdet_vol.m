function varargout = dfdet_vol(op,varargin)
%DFDET_VOL  Volatility and persistence regimes in the general engine: the model, the tilt, the transform, the checks.
%
%   VM = DFDET_VOL('model',D,V,cr,start,g0,gls,v0): the transform model of the
%   design D (with DFDET_DETREND_FIRST the statistic that detrends first) with
%   the regimes of the volatility path V (VOLPATH) at the c of the regimes cr
%   (cached; vol_model of DFDet_vol.jl): the design refined at the regime
%   boundaries (its prepared design P, by least squares or, with gls = cbar,
%   by GLS detrending, and Mi, the engine's model of P with the start, the
%   inner model of the paths), np pieces of lengths L with omega^2 (om2) and
%   c on each, the vector ell = (J(1), J at the jumps of d, z = int g J, J at
%   the regime boundaries that are not jumps of d) of dimension Ne, bpos(j)
%   the position in ell of J at the boundary between pieces j and j + 1 (0:
%   not in ell), the padded matrices A, B, Lx of the design and the
%   eigenvalues Dv; the fields start, g0, r, rp, rm and key that the outer
%   contour (DFDET_CPATH) reads.  start 'zero' (the default) or 'stationary':
%   a normal start J(0) ~ N(0, v0) independent of W, v0 the variance of the
%   stationary law of the first regime, om_1^2/(-2 c_1) with c_1 < 0, unless
%   given (v0 = NaN, the default); J(0) is then one more functional, N =
%   Ne + 1, and A, B, Lx, Dv are the matrices of the engine's stationary start
%   (A bordered by b/2 and -1/2, B by zeros); otherwise N = Ne.  GLS
%   detrending and the statistic that detrends first take the zero start
%   (DFDET_MODEL refuses another).  A fixed start (g0 ~= 0) is refused: its
%   transform under regimes is DFDet.vol_qf of DFLimit.jl, and its inversion
%   on the engine's paths is not dependable under regimes.  DFDET_MGF,
%   DFDET_SIGMA and DFDET_TSTAR take VM in place of an engine model.
%   [th,gam,ex,S] = DFDET_VOL('tilt',VM,t): the theta_j at t, the diagonal
%   tilt Gamma (g_last at J(1), g_j - g_{j+1} at the boundary between pieces j
%   and j + 1, g_j = (c_j + theta_j)/(2 omega_j^2); with the stationary start
%   -g_1 at J(0)), the exponent sum_j (c_j + theta_j) L_j and, when asked, the
%   covariance S of ell under the tilted process (DFDET_SIGMA; with the
%   stationary start that of (ell, J(0)), [Sigma + v0 h h', v0 h; v0 h', v0],
%   h the mean of ell under the tilted process started at one).
%   d = DFDET_VOL('dval',VM,s,t): D(s,t) = exp{sum_j (c_j + theta_j) L_j}
%   det{I - 2 (s A + t B + Gamma) Sigma}, M(s,t) = e^{-s/2} D(s,t)^{-1/2}
%   (with the stationary start the matrices of (ell, J(0))).
%   [VM,c1] = DFDET_VOL('args',fn,D,c,opts,start,g0,gls): the checks of the
%   arguments with opts.vol (vol_args of DFDet_vol.jl) and the model; VM = []
%   for one regime, which the engine itself computes at c1 with the start and
%   the detrending of the call.  The variants: least squares detrending with
%   a zero or a stationary start (that of the first regime), and with a zero
%   start GLS detrending (opts.gls = cbar) and the statistic that detrends
%   first (DFDET_DETREND_FIRST); a fixed start is refused with more than one
%   regime.
%   K = DFDET_VOL('kernel',VM,r,u): the covariance of the path with the
%   regimes, with the term of the stationary start (the Nystrom estimate of
%   t_*).
%   phi = DFDET_VOL('startphi'): the angle at which the outer hyperbola of the
%   volatility model starts (70 degrees, where the engine's designs start at
%   45; DFDET_CPATH steepens it from there); hf = DFDET_VOL('hfac',hf,r): the
%   factor of the outer step of the rule r under regimes, 2 for the rules
%   with the thin outer step ('thin' and 'outer', 1.5 for the engine's
%   designs) and hf otherwise (VOL_PHI, VOL_HFAC, start_phi and vol_hfac of
%   DFDet_vol.jl; the reason in the comment of vol_consts).
%
%   See also VOLPATH, DFDET_CDF.
switch op
    case 'startphi', varargout = {vol_consts()};
    case 'hfac'
        [~,hf] = vol_consts(); r = varargin{2};
        if ~any(strcmp(r,{'thin','outer'})), hf = varargin{1}; end
        varargout = {hf};
    case 'model', varargout = {vol_model(varargin{:})};
    case 'tilt'
        if nargout > 3, [a,b,c,d] = vol_tilt(varargin{:}); varargout = {a,b,c,d};
        else, [a,b,c] = vol_tilt(varargin{:}); varargout = {a,b,c}; end
    case 'dval', varargout = {vol_dval(varargin{:})};
    case 'args', [a,b] = vol_args(varargin{:}); varargout = {a,b};
    case 'kernel', varargout = {vol_kernel(varargin{:})};
    otherwise, error('dfdet:args','dfdet_vol: unknown operation %s.',op);
end
end

function [phi,hfac] = vol_consts()
% VOL_PHI, VOL_HFAC of DFDet_vol.jl: the outer hyperbola of the volatility model starts at 70 degrees (start_phi; the
% engine's designs at 45, and DFDET_CPATH steepens either while a root group leaves its side of the imaginary axis), and
% the thin outer step of the rules 'thin' and 'outer' (hfac 1.5 for the engine's designs) is 2 times the step of that
% hyperbola (vol_hfac).  Under regimes a root of the transform can come close to the imaginary axis far out on a flat
% hyperbola (an early regime of high variance with a break in the deterministic terms or a persistence that changes: the
% smallest real part of the positive group fell from 15 at the vertex to 1.2 at |t| = 1200 with a break in the trend at
% 0.5, three regimes of standard deviations 2, 1, 0.5 and c = -5, 0, -10); its part of the outer integrand is then not
% damped while its phase turns faster than the outer step resolves, and the trapezoidal rule loses the same amount at
% every x.  On 1456 cases (14 designs, 13 paths of two or three regimes, c constant and changing, zero and stationary
% start; November 2026, against the 80-degree hyperbola with half the step and the inner nodes of the present rule) the
% 45-degree hyperbola with the thin step was off by more than 1e-5 in 17 cases and more than 1e-8 in 71 (3.8e-4 at most;
% in the far upper tail the range check of the distribution function then refused values above 1 + 1e-8), with the
% present step by more than 1e-5 in 2 (6.8e-5); starting at 70 degrees with the factor 2 the largest error was 1.7e-8 on
% 195 nodes on average, against 140 of the 45-degree hyperbola with the thin step and 209 with the present one (the
% factor 2.5 lost 2.6e-6 in one case; 65 degrees with the factor 2 lost 2.3e-5).  Without regimes the thin step on the
% 45-degree hyperbola is within 1e-13 of the same reference on these designs.
phi = 70*pi/180; hfac = 2;
end

function VM = vol_model(D,V,cr,start,g0,gls,v0)
if nargin < 4 || isempty(start), start = 'zero'; end
if nargin < 5 || isempty(g0), g0 = 0; end
if nargin < 6, gls = []; end
if nargin < 7 || isempty(v0), v0 = NaN; end
start = lower(char(start)); g0 = double(g0); v0 = double(v0);
if ~any(strcmp(start,{'zero','fixed','stationary'})), error('dfdet:args','dfdet_vol: start must be ''zero'', ''fixed'' or ''stationary''.'); end
if strcmp(start,'fixed') && g0 == 0, start = 'zero'; end
if strcmp(start,'fixed'), error('dfdet:args','dfdet_vol: %s',dfdet_msg('vol_fixed')); end
cr = double(cr(:).');
if isempty(gls), gk = 'NaN'; else, gls = double(gls); gk = sprintf('%.17g',gls); end
key = sprintf('%s|vol%s|%s|%s|%s|%s|%.17g',dfdet_key(D),sprintf('%.17g,',V.lam),sprintf('%.17g,',V.om),sprintf('%.17g,',cr),start,gk,v0);
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
% the engine's model of the refined design with the start, prepared by least squares or by GLS detrending (DFDET_MODEL
% refuses a start other than zero with GLS detrending or with the statistic that detrends first; the model of the
% stationary start does not depend on c, and the c of the first regime enters through v0 below)
Mi = dfdet_startmodel(Dr,-1,start,0,gls);
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
N = Ne; Dv = P.Dv(:); r = P.r; rp = P.rp; rm = P.rm; vs = 0;
if strcmp(Mi.start,'stationary')
    % a normal start: J(0) the last functional, the matrices of the engine's stationary start (DFDET_MODEL)
    vs = v0; if isnan(vs), vs = om2(1)/(-2*c(1)); end
    if ~(isfinite(vs) && vs > 0), error('dfdet:args','dfdet_vol: the variance of a normal start must be positive and finite (the stationary start of the first regime needs c < 0 there)'); end
    N = Ne + 1;
    Ax = zeros(N); Bx = zeros(N);
    Ax(1:Ne,1:Ne) = A; Bx(1:Ne,1:Ne) = B;
    Ax(1:P.N,N) = P.b/2; Ax(N,1:P.N) = P.b'/2; Ax(N,N) = -0.5;
    [U,E] = eig((Ax + Ax')/2); [ev,ix] = sort(diag(E)); U = U(:,ix);
    keep = abs(ev) > 1e-10*max(abs(ev));
    Dv = ev(keep); Lx = U(:,keep); A = Ax; B = Bx;
    if ~(numel(Dv) == Mi.r && nnz(Dv > 0) == Mi.rp && nnz(Dv < 0) == Mi.rm)
        error('dfdet:internal','dfdet_vol: internal error: the inertia of the numerator with the start differs from that of the engine''s model');
    end
    r = Mi.r; rp = Mi.rp; rm = Mi.rm;
end
VM = struct('vm',true,'P',P,'Mi',Mi,'lam',V.lam,'om',V.om,'cr',cr,'np',np,'L',L,'om2',om2,'c',c,'Ne',Ne,'N',N,'bpos',bpos, ...
    'A',A,'B',B,'Lx',Lx,'Dv',Dv,'r',r,'rp',rp,'rm',rm,'start',Mi.start,'g0',0,'v0',vs,'a',Mi.a,'key',key);
dfdet_cache('put','volmodel',key,VM);
end

function [th,gam,ex,S] = vol_tilt(VM,t)
% (vol_tilt! of DFDet_vol.jl; the covariance only when asked)
th = sqrt(complex(VM.c.^2 - 2*t*VM.om2));
g = (VM.c + th)./(2*VM.om2);
gam = complex(zeros(VM.N,1));
gam(1) = gam(1) + g(VM.np);
for j = 1:VM.np-1
    pos = VM.bpos(j);
    if pos > 0, gam(pos) = gam(pos) + g(j) - g(j+1); end
end
ex = sum((VM.c + th).*VM.L);
stat = strcmp(VM.start,'stationary');
if stat, gam(VM.N) = -g(1); end
if nargout < 4, return, end
if stat
    % the covariance of (ell, J(0)) from Sigma and the mean h of ell under the tilted process started at one (statsigma!)
    [Sg,h] = dfdet_sigma(VM,th); n = VM.Ne; N = VM.N; v = VM.v0;
    S = complex(zeros(N));
    S(1:n,1:n) = Sg + v*(h*h.'); S(N,1:n) = v*h.'; S(1:n,N) = v*h; S(N,N) = v;
else
    S = dfdet_sigma(VM,th);
end
end

function d = vol_dval(VM,s,t)
[~,gam,ex,S] = vol_tilt(VM,complex(t));
M = eye(VM.N) - 2*((s*VM.A + t*VM.B + diag(gam))*S);
d = exp(ex)*det(M);
end

function K = vol_kernel(VM,r,u)
% e^{C(r) + C(u)} (v0 + int_0^{min(r,u)} omega^2 e^{-2C(v)} dv), C(r) = int_0^r c and v0 the variance of the stationary
% start (zero otherwise), piece by piece in closed form, each term with its exponents combined so that nothing overflows
% far from the unit root (vol_kernel of DFDet_vol.jl); r a column, u a row
R = r + 0*u; U = u + 0*r;
br = VM.P.br; np = VM.np;
Cb = [0 cumsum(VM.c.*VM.L)];
E = cat_(R) + cat_(U); m = min(R,U);
v0 = 0; if strcmp(VM.start,'stationary'), v0 = VM.v0; end
K = v0*exp(E);
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
% the checks of the arguments of the volatility extension (vol_args of DFDet_vol.jl, in its order) and its model
V = opts.vol;
if ~(isstruct(V) && isfield(V,'lam') && isfield(V,'om')), error('dfdet:args','%s: vol must be a volatility path, volpath(breaks, omega)',fn); end
k = numel(V.om);
if isscalar(c), cr = repmat(double(c),1,k); else, cr = double(c(:).'); end
if numel(cr) ~= k, error('dfdet:args','%s: c must be a number or a vector with one value per regime of the volatility path (%d), not %d',fn,k,numel(cr)); end
if ~all(isfinite(cr)), error('dfdet:args','%s: c must be finite',fn); end
if ~any(strcmp(start,{'zero','fixed','stationary'})), error('dfdet:args','%s: start must be ''zero'', ''fixed'' or ''stationary''',fn); end
if ~(isnumeric(g0) && isscalar(g0) && isreal(g0) && isfinite(g0)), error('dfdet:args','%s: g0 must be a finite real scalar',fn); end
zero0 = strcmp(start,'zero') || (strcmp(start,'fixed') && g0 == 0);
if ~any(strcmp(D.detrend,{'onestep','first'})), error('dfdet:args','%s: the statistic of the design must be ''onestep'' or ''first''',fn); end
first = strcmp(D.detrend,'first');
if ~isempty(gls) && ~(isnumeric(gls) && isscalar(gls) && isreal(gls) && isfinite(gls) && gls < 0)
    error('dfdet:args','%s: GLS detrending needs opts.gls = cbar < 0 (a finite number)',fn);
end
if ~isempty(gls) && first, error('dfdet:args','%s: %s',fn,dfdet_msg('first_gls')); end
if ~isempty(gls) && ~zero0, error('dfdet:args','%s: %s',fn,dfdet_msg('gls_start')); end
if first && ~zero0, error('dfdet:args','%s: %s',fn,dfdet_msg('first_start')); end
m = 'auto'; if isfield(opts,'method') && ~isempty(opts.method), m = lower(char(opts.method)); end
if ~any(strcmp(m,{'auto','paths'}))
    error('dfdet:args','%s: with a volatility path (vol) the method is the inversion on paths: method must be ''auto'' or ''paths'' (the line integral is not implemented for the volatility transform)',fn);
end
if isscalar(c), cs = sprintf('%g',c); else, cs = ['[' strjoin(arrayfun(@(v) sprintf('%g',v),cr,'UniformOutput',false),', ') ']']; end
if ~all(cr <= 0), error('dfdet:args','%s: with a volatility path (vol) every c must be at most 0 (above the unit root the extension is not implemented); c = %s',fn,cs); end
if ~all(cr >= -480), error('dfdet:args','%s: with a volatility path (vol) every c must be at least -480 (the range of the general engine); c = %s',fn,cs); end
if strcmp(start,'stationary') && ~(cr(1) < 0), error('dfdet:args','%s: %s; c = %s',fn,dfdet_msg('vol_stat'),cs); end
c1 = cr(1);
if k == 1, VM = []; return, end
if ~(zero0 || strcmp(start,'stationary')), error('dfdet:args','%s: %s',fn,dfdet_msg('vol_fixed')); end
if zero0, st = 'zero'; else, st = 'stationary'; end
VM = vol_model(D,V,cr,st,0,gls);
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
