function out = dfdet_lmpath(D,c)
%DFDET_LMPATH  The outer path of the denominator of the LM statistic, and P(V <= v) on it.
%   CP = DFDET_LMPATH(D,c): the path of the design of the LM statistic (DFDET_PREP, kind 'lm'; an error unless J_0(1)
%   vanishes, which needs a linear trend among the trend functions); P = DFDET_LMPATH(CP,v) = P(V <= v), the single
%   inversion of M_V on it.  Shared by LM_CDF, LM_QUANTILE, LM_POWER.
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
pk = [dfdet_key(D) '|lm'];
P = dfdet_cache('get','prep',pk);
if isempty(P)
    [P,a1] = dfdet_prep(D,'lm');
    if ~(max(abs(a1)) < 1e-9), error('dfdet:args','the LM statistic: the law of -1/(2 sqrt V) needs a linear trend among the trend functions (J_0(1) does not vanish).'); end
    dfdet_cache('put','prep',pk,P);
end
Mo = dfdet_model(P); Mo.key = [pk '|zero|0'];
out = dfdet_cpath(Mo,c);
end
