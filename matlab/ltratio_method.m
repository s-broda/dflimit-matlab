function [dbl,tails,ctr] = ltratio_method(cas,c,c0,g0,x,verified,cond)
%LTRATIO_METHOD  Automatic choice of the representation (opts.method = 'auto').
%
%   With verified = false (the default of the programs since November 2026: one
%   evaluation by one method) the rule is the fastest representation that one
%   evaluation gives to well below 1e-5, by the timings of these programs on one
%   thread (November 2026).  With a fixed start the contour integral in cases
%   DU, DT, except where it fails (c > 1 with c0 < 0; in case DU c > 1.5 with
%   c0 > 3), and the conditional integral there and in case DK (the only route
%   here).  With a zero or stationary start the line integral at c0 = 0 (but
%   the conditional integral in cases DK, DU for -30 <= c <= -5, 1.7 to 2.3
%   times faster here), with a zero start and c0 > 0 (where one evaluation of
%   the conditional integral is off by up to 1.5e-2 at x < 0 for c < 0; but the
%   conditional integral at x < 0 for c >= 0 in cases DU, DT), in case DK with
%   c0 < 0 for c > -5, in case SK for c > -0.5, in cases SU, ST for c0 <= -1 at
%   c >= -20 (1.7 to 3.4 times faster here in both tails, F and f at the 5% and
%   the 95% points, c from -0.01 to -20, November 27, 2026; at c0 = -0.5 1.29 to
%   3.3 times, with c0 > 0 within a factor of 1.3 in one tail or the other, in
%   case SK the conditional integral the faster) and far from the unit root
%   (c < -100 with a zero start, c < -50 with a stationary one); the
%   conditional integral elsewhere.  With cond = true (the seventh argument)
%   the rule before November 27, 2026, the conditional integral in cases SU, ST
%   with c0 ~= 0: LTRATIO_DIST sends the general engine's values below 1e-3
%   there to the conditional integral, whose error is relative.  The rule
%   differs from that of the C and Julia engines (automethod, dfl_automethod)
%   in case DK with a fixed start (they take the line integral, which these
%   programs do not have), at c0 = 0, where their line integral is the faster
%   one, and in cases SU, ST (the Julia engine takes the line integral for c0
%   <= -1/2 at c >= -20 and in case SU for 1/2 <= c0 <= 10 at c >= -2, the C
%   engine at every c0).  tails is false (no recomputation).  The
%   rule below is that of verified = true (and the one before November 2026):
%
%   [dbl,tails,ctr] = LTRATIO_METHOD(cas,c,c0,g0,x): dbl is a logical array of
%   the size of x, true where the conditional integral (TR_LIMDIST) is the
%   faster representation and false where the line integral (DFX_CORE) is;
%   ctr, of the same size, is true where the contour integral (TR_CONTOUR)
%   serves, that is at every x with a fixed start g0 ~= 0 in cases DU, DT
%   (dbl is false there); tails is false where the conditional integral is
%   not used at all, not even for the tails of the line integral or of the
%   contour integral (LTRATIO_DIST).  The rule, from the
%   per-point timings of the MATLAB programs themselves on one thread
%   (validation/method_sweep.m, November 2026; c in [-5,5], c0 in {-5,-1,1,3},
%   x in [-4,32], the six cases):
%
%     g0 ~= 0, DU, DT    contour (the contour integral, a double integral with
%                        an explicit integrand, the choice of the C and Julia
%                        engines too; its values below 1e-3 and its flagged
%                        values are recomputed by the conditional integral,
%                        LTRATIO_DIST, but only where that is checked here,
%                        c >= -20, or c >= -40 with (c g0)^2/2 <= 100, tails
%                        false elsewhere: at c = -30, g0 = 0.6 in case DU the
%                        conditional integral is off by 2e-9 relative at
%                        F = 1e-4 without a flag);
%     g0 ~= 0, DK        double (the only route with a fixed start here);
%     c0 == 0            single: 1.2 to 60 times faster everywhere;
%     cases DK, DU, DT, SK with c0 ~= 0
%                        double, except for c0 > 0 at x >= 0, where the
%                        line integral has no breakpoint and wins by 3 to
%                        30 times; with c0 < 0 the double wins by 2 to 10
%                        times at x < 0 and the line integral's cost
%                        spikes to 5-11 seconds at x of 4 to 8 (its saddle
%                        factor), while the double stays below 0.5;
%     cases SU, ST with c0 ~= 0
%                        single, except for c0 < 0 at x >= 0, where the
%                        line integral's saddle factor costs up to 2.7
%                        seconds against 0.5 by the double; at x < 0 the
%                        single wins by 1 to 4 times in these cases, whose
%                        conditional integral costs 0.5 seconds a value.
%
%   In case DK with c0 > 0 the point x = -sqrt(2 c0), where the density is
%   infinite, is sent to the line integral, which returns Inf there
%   without a warning.
%
%   Far from the unit root the line integral serves alone (dbl and tails
%   false): with a stationary start for c < -50, where the conditional
%   integral flags most points (two of five at the mean and 2 and 5 standard
%   deviations either side at c = -80, all at c = -200) at 10 to 100 times
%   the cost; with c0 ~= 0 and a zero start (g0 = 0) for c < -100, where
%   the density of case DK is off by up to 1e-9 at some points without a
%   flag (at c = c0 = -300 and -480); and with c0 > 0 and a zero start for
%   c < -20, where the checked conditional integral was off by up to 1e-5 in
%   cases DU, DT and verified = true returned such a value at c = -100,
%   c0 = 3 (2.1e-8 in case DU, 4.2e-8 in case DT, November 26, 2026).  The
%   same rule is in the C and Julia engines (dfl_autodouble, autodouble).
%
%   The rule of the C and Julia engines (dfl_automethod, automethod) is the
%   one of the first four lines for every start (with a fixed start in case
%   DK they take the line integral, which these programs do not have); it
%   differs here in cases
%   SU and ST at x < 0, where the MATLAB conditional integral is the slower one,
%   and only there.  Where the two representations are within a factor of
%   two of each other the rule follows the C engine (case DT with c0 = -5 at
%   x >= 0 and c <= 2, case DK with c0 = -5 at c <= 0 and x >= 4, where the
%   line integral is 1.2 to 2.8 times faster; case SK with c0 = -1 at
%   x >= 16).
%
%   See also LTRATIO_DIST, LTRATIO_CDF, LTRATIO_PDF.
tails = true; ctr = false(size(x));
if strcmp(cas,'DG'), dbl = false(size(x)); tails = false; return; end   % the GLS-detrended case: the line integral only
if nargin > 5 && ~verified
    % the one-evaluation rule (automethod of julia/DFExact.jl, dfl_automethod of the C engine)
    tails = false; dbl = false(size(x));
    if g0 ~= 0
        fails = (c > 1 && c0 < 0) || (strcmp(cas,'DU') && c > 1.5 && c0 > 3);   % (where the contour integral fails)
        if any(strcmp(cas,{'DU','DT'})) && ~fails, ctr = true(size(x)); else, dbl = true(size(x)); end
        return
    end
    if c0 == 0
        dbl(:) = any(strcmp(cas,{'DK','DU'})) && c >= -30 && c <= -5;
    elseif (cas(1) == 'S' && c < -50) || c < -100 || (c0 > 0 && cas(1) == 'D' && c < -20)
        % far from the unit root: the line integral (as dfl_autodouble, autodouble)
    elseif cas(1) == 'S'
        % (cases SU, ST with c0 <= -1 at c >= -20: the line integral, 1.7 to 3.4 times faster here than the conditional
        % integral in both tails, November 27, 2026; cond = true, the rule before, keeps the conditional integral there:
        % the tails that LTRATIO_DIST sends back from the general engine)
        su = any(strcmp(cas,{'SU','ST'})) && c0 <= -1 && c >= -20 && ~(nargin > 6 && cond);
        dbl(:) = ~((strcmp(cas,'SK') && c > -0.5) || su);
    elseif strcmp(cas,'DK') && c > -5
        % case DK with c0 ~= 0 near or above the unit root: the line integral
    elseif c0 > 0
        dbl = (c >= 0) & (x < 0);   % (a zero start with c0 > 0: the conditional integral at x < 0 for c >= 0)
    else
        dbl(:) = true;
    end
    return
end
if g0 ~= 0 && any(strcmp(cas,{'DU','DT'}))                % the contour integral; the conditional integral for its
    dbl = false(size(x)); ctr = true(size(x));             % tails and flagged values only where it is checked here
    tails = c >= -20 || (c >= -40 && (c*g0)^2/2 <= 100);   % (LTRATIO_ARGS)
    return
end
if g0 ~= 0, dbl = true(size(x)); return; end
if (cas(1) == 'S' && c < -50) || (c0 ~= 0 && c < -100) || (c0 > 0 && cas(1) == 'D' && c < -20)   % far from the unit root: the line integral alone
    dbl = false(size(x)); tails = false; return
end
if c0 == 0, dbl = false(size(x)); return; end
if any(strcmp(cas,{'SU','ST'}))
    dbl = (c0 < 0) & (x >= 0);
else
    dbl = ~((c0 > 0) & (x >= 0));
end
dbl = logical(dbl) & true(size(x));
% case DK with c0 > 0 at x = -sqrt(2 c0), where the density is infinite: the
% line integral returns Inf there without a warning, the conditional integral
% with one
if strcmp(cas,'DK') && c0 > 0, dbl(abs(x + sqrt(2*c0)) <= 1e-12*max(1,abs(x))) = false; end
end
