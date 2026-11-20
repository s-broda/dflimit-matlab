function lam = lam_gls(cb)
%LAM_GLS  The weight lambda = (1 - cbar)/(1 - cbar + cbar^2/3) of the GLS detrending of Elliott, Rothenberg and Stock
%   (1996), trend case, at the detrending parameter cbar (shared by DFX_CORE and GLS_PARTS).
lam = (1 - cb) / (1 - cb + cb^2 / 3);
end
