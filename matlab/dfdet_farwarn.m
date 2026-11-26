function dfdet_farwarn(fn,D,c,start,g0,gls)
%DFDET_FARWARN  The warning of the brittle region far above the unit root (far_above_warn of DFDet.jl).
%
%   DFDET_FARWARN(fn,D,c,start,g0,gls) warns (identifier dfdet:brittle) where
%   one evaluation of the inversion on paths is known to be brittle: above
%   c = 5.5, beyond the range in which the rules of one evaluation were swept,
%   for a design beyond the classical cases (trend(-1), trend(0), trend(1) in
%   the one-step regression with least squares detrending, and trend(1)
%   GLS-detrended at cbar = -13.5 with a zero start) and for a fixed start with
%   deterministic terms.  There one evaluation was off by more than 1e-5 in
%   the upper tail (November 26, 2026, against the paths with a quarter of the
%   outer step: Perron's model C at 0.3 by 1.5e-5 in F at c = 7, the linear
%   trend with a fixed start, g0 = 1, by 8.7e-5 at c = 7); opts.verified = true
%   checks the values.
%
%   See also DFDET_DIST, DFDET_ICDF.
if ~(c > 5.5), return, end
fixed = strcmp(start,'fixed') && g0 ~= 0;
classical = false;
for p = -1:1
    T = dfdet_trend(p);
    if isequal(size(D.C),[1 1]) && isequal(D.br,T.br) && isequal(size(D.C{1}),size(T.C{1})) && isequal(D.C{1},T.C{1}) && ~D.loc && strcmp(D.detrend,'onestep')
        classical = isempty(gls) || (p == 1 && gls == -13.5 && strcmp(start,'zero'));
        break
    end
end
if classical && ~(fixed && size(D.C{1},2) >= 1), return, end
warning('dfdet:brittle',['%s: far above the unit root (c = %g > 5.5) the inversion on paths by one evaluation may be off by more than 1e-5 in the upper ' ...
    'tail for a design beyond the classical cases or a fixed start with deterministic terms (Perron''s model C at c = 7 by 1.5e-5 in F, the linear ' ...
    'trend with a fixed start at c = 7 by 8.7e-5); opts.verified = true checks it.'],fn,c);
end
