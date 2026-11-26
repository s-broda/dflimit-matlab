function out = dfdet_lmpath(D,c)
%DFDET_LMPATH  The outer path of the denominator of the LM statistic, and P(V <= v) on it.
%   CP = DFDET_LMPATH(D,c): the path of the law of V alone for the design of the LM statistic (DFDET_PREP, kind 'lm'; an
%   error unless the span of the trend functions contains the constant, for which alone the limit -1/(2 sqrt V) holds,
%   and unless J_0(1) vanishes, which needs a linear trend among them): the path of the t-ratio, except
%   that above the unit root, where that is the vertical line, the hyperbola of V alone at pi/4 (on the line P(V <= v) was
%   off by 6e-3 at c = 5 with a linear trend, and F exceeded 1).  c at most 10: above, the inversion of this V moves between
%   contours by more than 1e-8 (2e-7 at c = 12, 3e-5 at c = 15), and an error is raised.  P = DFDET_LMPATH(CP,v) =
%   P(V <= v), the single inversion of M_V on it.  Shared by LM_CDF, LM_QUANTILE, LM_POWER.
if isfield(D,'Mo')
    CP = D; v = c;
    Mo = CP.Mo; r = size(CP.R,1); rp = Mo.rp;
    Z = CP.R; Z(rp+1:r,:) = -Z(rp+1:r,:);
    lg = log(Z);
    for k = 2:size(Z,2), lg(:,k) = lg(:,k) + 2i*pi*round((imag(lg(:,k-1)) - imag(lg(:,k)))/(2*pi)); end
    lM = CP.lP(:) - sum(lg,1).'/2;
    out = 1 - sum(CP.wt.*imag(exp(lM - CP.t*v).*CP.dt./CP.t));
    return
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('dfdet:args','the LM statistic: c must be a finite real scalar.'); end
if c > 10, error('dfdet:args','the LM statistic: c must be at most 10 (above, its distribution is not computed accurately: off by 2e-7 at c = 12 and by 3e-5 at c = 15).'); end
pk = [dfdet_key(D) '|lm'];
P = dfdet_cache('get','prep',pk);
if isempty(P)
    % (the limit -1/(2 sqrt V) holds only for trend functions whose span contains the constant: the LM regression estimates
    % the level from the first observation; r alone, r with e^r, a level dummy with r gave the limit of the linear trend)
    if ~span_has_constant(D)
        error('dfdet:args',['the LM statistic: the limit computed here holds only for trend functions whose span contains the ' ...
              'constant (the LM regression estimates the level from the first observation, and the limit of the statistic is ' ...
              'otherwise a different one); the design %s does not contain it.'],D.name);
    end
    [P,a1] = dfdet_prep(D,'lm');
    if ~(max(abs(a1)) < 1e-9), error('dfdet:args','the LM statistic: the law of -1/(2 sqrt V) needs a linear trend among the trend functions (J_0(1) does not vanish).'); end
    dfdet_cache('put','prep',pk,P);
end
Mo = dfdet_model(P); Mo.key = [pk '|zero|0'];
out = dfdet_cpath(Mo,c);
if isnan(out.h)
    % (the vertical line of the t-ratio does not resolve e^{-t v}: the hyperbola of V alone, at pi/4, or at the angle of the
    % fixed start's rule below c = -60, as DFDET_V_CDF)
    if c < -60, ph = dfdet_cpath(Mo,c,struct('fixedphi',true,'ts',out.tstar,'tau',out.tau,'hfac',1)); else, ph = pi/4; end
    out = dfdet_cpath(Mo,c,struct('phi',ph));
end
end

function tf = span_has_constant(D)
% whether the constant function lies in the span of the regressors of D (span_has_constant of DFDet.jl), decided on values:
% the regressors at 2 (deg + 1) Chebyshev points of the first kind inside each piece (each with its own piece's
% coefficients, so that no break is evaluated), the least squares fit of the constant 1 by them (QR with column pivoting)
% and its largest residual, against 1e-6.  For polynomial trends, breaks and dummies the residual is zero to rounding; for
% a function design and a design from a matrix the check is made on their interpolants, so that a span that contains the
% constant through a combination (sin^2 and cos^2, say) gives a residual at their level.  (The coefficients themselves are
% not compared: the high-order local Taylor coefficients of a function design carry rounding of up to 1e-2, the intercept
% with a sine and a cosine 1.3e-2, that cancels in the values.)  No regressors: no.
K = size(D.C{1},2);
if K == 0, tf = false; return; end
np = numel(D.br) - 1; m = 2*size(D.C{1},1);
X = zeros(np*m,K);
for j = 1:np
    a = D.br(j); b = D.br(j+1); off = 0; if D.loc, off = a; end
    u = (a + b)/2 + (b - a)/2*cos(pi*((1:m)' - 0.5)/m) - off;
    for i = 1:K, X((j-1)*m+(1:m),i) = polyval(flipud(D.C{j}(:,i)),u); end
end
[Q,R,p] = qr(X,0);
d = abs(diag(R)); rk = sum(d > min(size(X))*eps*d(1));
beta = zeros(K,1); beta(p(1:rk)) = R(1:rk,1:rk) \ (Q(:,1:rk)'*ones(np*m,1));
tf = max(abs(X*beta - 1)) <= 1e-6;
end
