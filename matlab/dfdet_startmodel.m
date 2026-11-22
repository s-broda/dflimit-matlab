function Mo = dfdet_startmodel(D,c,start,g0,gls)
%DFDET_STARTMODEL  The model of a design with its start (the prepared design cached by its coefficients).
%   Mo = DFDET_STARTMODEL(D,c,start,g0,gls): start 'zero', 'fixed' (J(0) = g0) or 'stationary' (c < 0); gls = cbar
%   for GLS-detrended data ([] for least squares).  Mo.key identifies the model for the caches of the outer paths.
if nargin < 5, gls = []; end
if strcmp(start,'stationary') && c >= 0, error('dfdet:args','%s',dfdet_msg('stat_c')); end
dk = dfdet_key(D);
if isempty(gls), pk = [dk '|ls']; else, pk = sprintf('%s|gls%.17g',dk,gls); end
P = dfdet_cache('get','prep',pk);
if isempty(P)
    if isempty(gls), P = dfdet_prep(D); else, P = dfdet_prep(D,'gls',gls); end
    dfdet_cache('put','prep',pk,P);
end
Mo = dfdet_model(P,start,g0);
Mo.key = sprintf('%s|%s|%.17g',pk,Mo.start,Mo.g0);
end
