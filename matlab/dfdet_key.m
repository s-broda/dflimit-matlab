function k = dfdet_key(D)
%DFDET_KEY  A key of a design for the caches of the general engine: its coefficients, breaks, the form of its
%   coefficients and its statistic, written out in full precision.
c = cellfun(@(M) sprintf('%.17g,',size(M),M(:)),D.C,'UniformOutput',false);
k = [sprintf('%.17g,',D.br) '|' [c{:}] '|' sprintf('%d|',D.loc) D.detrend];
end
