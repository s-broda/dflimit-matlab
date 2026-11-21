function P = tr_linexact(b,a,cas,c,opts,isx)
%TR_LINEXACT  P(U - aV <= b) where the inversion of Gil-Pelaez is not accurate, by a conditional integral.
%
%   P = TR_LINEXACT(b,a,cas,c,opts) (elementwise in b, a finite), zero start, clamped to [0, 1]: in case 'DK' the
%   probability of the region {U <= aV + b} by the conditional integral given V (TR_REGION), in case 'DG' (detrended at
%   opts.cbar, default -13.5) the conditional integral given R = U + 1/2 (TR_LINCOND).  TR_LINEXACT(x,a,cas,c,opts,true)
%   takes x = b + 1/2 in place of b (POWER_ENVELOPE: in case DG near cbar = 0 the scale of U + 1/2 - aV falls far
%   below the spacing of the doubles near b = -1/2).  The route of TR_GILPELAEZ where its transform has not decayed by
%   the end of its grid or c > 0.  As lin_exact and lin_exact_x of julia/DFExact.jl.
%
%   See also TR_GILPELAEZ, TR_LINCOND, TR_REGION, LTRATIO_LINCDF.
if nargin < 5 || isempty(opts), opts = struct; end
if nargin < 6 || isempty(isx), isx = false; end
cbar = -13.5; if isfield(opts,'cbar') && ~isempty(opts.cbar), cbar = opts.cbar; end
P = zeros(size(b));
for k = 1:numel(b)
    bk = b(k);
    if ~isfinite(bk), if isnan(bk), P(k) = NaN; else, P(k) = double(bk > 0); end, continue, end
    if strcmpi(cas,'DK')
        if isx, b0 = bk - 0.5; else, b0 = bk; end
        P(k) = min(1,max(0,tr_region(struct('type','line','a',a,'b',b0),c)));
    else
        if isx, x = bk; else, x = bk + 0.5; end
        P(k) = tr_lincond(x,a,'DG',c,cbar);
    end
end
end
