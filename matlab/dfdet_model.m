function Mo = dfdet_model(P,start,g0)
%DFDET_MODEL  A prepared design with the start of the Ornstein-Uhlenbeck process.
%
%   Mo = DFDET_MODEL(P,start,g0): start 'zero' (J(0) = 0, the default),
%   'fixed' (J(0) = g0; y_0 = g0 sigma sqrt(n)) or 'stationary' (J(0) ~
%   N(0, -1/(2c)) independent of W, c < 0).  With J(0) = x the functionals gain
%   two terms, U + 1/2 = ell' A ell - x^2/2 + x b' ell, b = Pm d(0) on the block
%   of z, and V is unchanged.
%      'stationary': ell_ext = (ell, J(0)) is Gaussian with mean zero under the
%   tilt, U + 1/2 = ell_ext' A_ext ell_ext with A_ext = [A b/2; b'/2 -1/2],
%   B_ext = [B 0; 0 0] and the extended covariance [Sigma + v h h', v h;
%   v h', v], v = -1/(2c); Girsanov's exponent puts -gamma on the last
%   diagonal entry.  The fields N, A, B, L, Dv, r, rp, rm are the extended
%   ones.
%      'fixed': b lies in the range of A for every design; with 2 A yb = b and
%   zeta = ell + x yb, U + 1/2 + x^2 kappa = zeta' A zeta, kappa = 1/2 +
%   yb' A yb, so that Mt = e^{a s} M(s,t), a = 1/2 + g0^2 kappa, is a bounded
%   function times algebraic factors on the paths.  kappa is exactly 0 or 1/2:
%   0 when the span of d contains a continuous phi with phi(0) ~= 0 and phi' in
%   the span (the constant, for instance: a = 1/2), 1/2 otherwise (no
%   deterministic terms: a = (1 + g0^2)/2).
%   Mo.a is the shift a of Mt = e^{a s} M(s,t), Mo.yb, Mo.kappa those of the
%   fixed start.
%
%   See also DFDET_PREP, DFDET_MGF.
if nargin < 2 || isempty(start), start = 'zero'; end
if nargin < 3 || isempty(g0), g0 = 0; end
start = lower(char(start)); g0 = double(g0);
if ~isnan(P.cbar) && ~(strcmp(start,'zero') || (strcmp(start,'fixed') && g0 == 0)), error('dfdet:args','%s',dfdet_msg('gls_start')); end
if strcmp(P.detrend,'first') && ~(strcmp(start,'zero') || (strcmp(start,'fixed') && g0 == 0)), error('dfdet:args','%s',dfdet_msg('first_start')); end
if strcmp(start,'zero') || (strcmp(start,'fixed') && g0 == 0)
    Mo = mk(P,'zero',0,P.N,P.A,P.B,P.L,P.Dv,0.5,zeros(P.N,1),0.5);
elseif strcmp(start,'fixed')
    % (b lies in the range of A and kappa is exactly 0 or 1/2 for every design; the two checks below are internal
    % assertions, which only a numerical failure of the decomposition of A could trigger)
    cb = P.L'*P.b;
    res = max([0; abs(P.b - P.L*cb)])/max(1,max([0; abs(P.b)]));
    yb = P.L*(cb./P.Dv)/2;
    kraw = 0.5 + yb'*(P.A*yb);
    if res > 1e-8, error('dfdet:internal','dfdet_model: internal error: b is not in the range of A to 1e-8 (residual %g), which contradicts the shift identity of the fixed start.',res); end
    if abs(kraw) < abs(kraw - 0.5), kappa = 0; else, kappa = 0.5; end
    if abs(kraw - kappa) > 1e-6, error('dfdet:internal','dfdet_model: internal error: kappa = %g is not 0 or 1/2 to 1e-6, which contradicts the shift identity of the fixed start.',kraw); end
    Mo = mk(P,'fixed',g0,P.N,P.A,P.B,P.L,P.Dv,0.5 + g0^2*kappa,yb,kappa);
elseif strcmp(start,'stationary')
    n = P.N; N = n + 1;
    A = zeros(N); B = zeros(N);
    A(1:n,1:n) = P.A; A(1:n,N) = P.b/2; A(N,1:n) = P.b'/2; A(N,N) = -0.5;
    B(1:n,1:n) = P.B;
    [U,E] = eig((A + A')/2); [ev,ix] = sort(diag(E)); U = U(:,ix);
    keep = abs(ev) > 1e-10*max(abs(ev));
    Mo = mk(P,'stationary',0,N,A,B,U(:,keep),ev(keep),0.5,zeros(N,1),0.5);
else
    error('dfdet:args','start must be ''zero'', ''fixed'' or ''stationary''.');
end
end

function Mo = mk(P,start,g0,N,A,B,L,Dv,a,yb,kappa)
Dv = Dv(:);
Mo = struct('P',P,'start',start,'g0',g0,'N',N,'A',A,'B',B,'L',L,'Dv',Dv,'r',numel(Dv),'rp',nnz(Dv > 0),'rm',nnz(Dv < 0), ...
    'a',a,'yb',yb,'kappa',kappa);
end
