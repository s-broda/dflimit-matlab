function out = dfdet_lmpath(D,c)
%DFDET_LMPATH  The outer path of the denominator of the LM statistic, and P(V <= v) on it.
%   CP = DFDET_LMPATH(D,c): the path of the law of V alone for the design of the LM statistic (DFDET_PREP, kind 'lm'; an
%   error unless J_0(1) vanishes, which needs a linear trend among the trend functions): the path of the t-ratio, except
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
    [P,a1] = dfdet_prep(D,'lm');
    if ~(max(abs(a1)) < 1e-9), error('dfdet:args','the LM statistic: the law of -1/(2 sqrt V) needs a linear trend among the trend functions (J_0(1) does not vanish).'); end
    dfdet_cache('put','prep',pk,P);
end
Mo = dfdet_model(P); Mo.key = [pk '|zero|0'];
out = dfdet_cpath(Mo,c);
if isnan(out.h)
    % (the vertical line of the t-ratio does not resolve e^{-t v}: the hyperbola of V alone, at pi/4, or at the angle of the
    % fixed start's rule below c = -100, as DFDET_V_CDF)
    if c < -100, ph = dfdet_cpath(Mo,c,struct('fixedphi',true,'ts',out.tstar,'tau',out.tau,'hfac',1)); else, ph = pi/4; end
    out = dfdet_cpath(Mo,c,struct('phi',ph));
end
end
