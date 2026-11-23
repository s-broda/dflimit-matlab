function gen = ltratio_engine(f,S,c,c0,x)
%LTRATIO_ENGINE  The choice of opts.engine = 'auto': the classical programs or the general engine.
%
%   gen = LTRATIO_ENGINE(f,S,c,c0,x) for a computation that the classical
%   programs accept (a case string, or a design with the classical case
%   S.case, LTRATIO_SPEC, within the checks of LTRATIO_ARGS; elsewhere the
%   general engine computes): true where the general engine (DFDET_*) is
%   taken.  f = 'dist' (the distribution function and the density, one
%   value per point of x), 'icdf', 'mom', 'lincdf' (x = a), 'envelope',
%   'adfmax', 'region'.  The rule is that of the Julia package
%   (DFLimit.autoengine) with its cells set by the timings of the two MATLAB
%   engines on the cells of its table and on further cells at the boundaries
%   (validation/engine_table_matlab.txt, the script
%   validation/engine_table_matlab.m): the faster engine where the two agree
%   to 1e-10, the more accurate one where they do not.  The general engine is
%   taken for the t-ratio
%     zero start, c0 = 0: cases DU, DT for c <= 1 (1.4 to 6 times faster),
%       case DK for c <= -50 (1.6 to 1.9 times; for -30 <= c <= 1 the
%       classical line integral is 1.3 to 2.2 times faster);
%     stationary start, c0 = 0, c <= -1 (1.5 to 4.4 times);
%     zero start, c0 < 0: case DU for -10 <= c <= 0, cases DK, DT for
%       -5 <= c <= 0 (1.2 to 2.8 times; from c = -10 in case DK, -20 in cases DU,
%       DT the classical conditional integral is 1.4 to 3 times faster);
%     zero start, c0 > 0, c <= 1 (2 to 6 times; in case DK just above the
%       tangency x* = -sqrt(2 c0) the classical line integral, the more
%       accurate there; at x = 0, the points of COEF_CDF and COEF_ICDF, the
%       classical line integral for c > -100 in cases DK, DU, DT, 3 to 8 times
%       faster there), and c0 ~= 0 for c < -100 (3.7 to 8 times);
%     stationary start, c0 ~= 0, c <= -3e-3 (2.4 to 70 times; nearer zero the
%       general engine is off by 1e-10 to 1e-7, as in the Julia package);
%     fixed start, case DK: c < -10 (the classical conditional integral is off
%       by 1.8e-7 at c = -20, g0 = 2, the general engine within 3e-14 of the
%       Julia engine) and, with c0 = 0, -1 <= c <= 0 with |g0| <= 1 (1.3 times;
%       with c0 < 0 near the unit root, c > -1, the classical programs, as in the
%       Julia and the C rule since November 2026: there the outer hyperbola of the
%       general engine does not resolve the upper tail);
%     fixed start, cases DU, DT far from the unit root: with c0 > 0 below c = -20
%       (below c = -10 where |c g0| > 20) and with any c0 below c = -20 where
%       |c g0| > 60, as the Julia and the C rule (November 2026): there the
%       contour integral's one evaluation was off, in the density by up to
%       2.2e-2 of its peak (case DU, c = -60, g0 = 1.5, c0 = 5) and in the
%       distribution function by 3.9e-5 (c0 = -2), its checked computation in the
%       density by up to 8.4e-6, and the general engine's density is the
%       derivative of the distribution function on which the engines agree;
%   the classical programs elsewhere (the cells with -100 <= c < -60, where
%   the general engine was off by up to 1.4e-9 on the 45-degree hyperbola,
%   took the classical programs until November 2026; its rules far from the
%   unit root now start at c = -60, within 5e-13 there): case DG, above the unit root (c > 1, and
%   c >= 4 with c0 < 0, where the general engine is the less accurate), a fixed
%   start in cases DU, DT otherwise (within 15% of each other).  Quantiles: the general
%   search at c0 = 0 except in case DG and for c > 1 (1.6 to 10 times faster,
%   with a fixed start too).  Moments: the classical programs (1.1 to 360
%   times faster).  The linear form: the general engine for the law of V
%   (a = NaN, 4 to 350 times faster, and the more accurate for c >= 5), in
%   case DK for c <= 0 with |a| <= 0.1 (where the classical programs take a
%   conditional integral; 1.9 to 3.4 times) and in case DG for c > 3.5, where
%   the classical inversion refuses.  The envelope: the general engine in
%   case DK for -1 < c < 0 (4 times faster), the classical programs otherwise
%   (4 to 28 times faster).  ADFmax: the general engine with a stationary
%   start for c <= -0.1 (as the Julia package) and with a zero start for
%   c <= -10 (4 times faster).  Regions: the classical programs.  With the
%   folder of the MEX file of the C engine ahead of these programs on the
%   path (LTRATIO_CROUTE('on'): LTRATIO_MEX and its drop-ins, which carry the
%   classical programs at 20 to 50 times their speed) the user-level routines
%   hand their calls with a design, or with a case string and opts.engine,
%   without options of the MATLAB engines to the C router, whose rule
%   (dfl_autoengine of c/dfl_det_route.c) is set by the timings of the two C
%   engines; for the other calls this function gives the classical engine
%   (on the C engine through the drop-ins), as before.
%
%   See also LTRATIO_SPEC, LTRATIO_CDF, DFDET_CDF.
cas = S.case; fixed = strcmp(S.start,'fixed'); stat = strcmp(S.start,'stationary');
if strcmp(f,'dist'), gen = false(size(x)); else, gen = false; end
if isempty(cas) || ltratio_croute('on'), return, end
% (above the range of the general engine, DFDET_CMAX: 20 without deterministic terms, 16 with them, 10 with them and a
% fixed start, where it refuses, the classical programs)
if c > dfdet_cmax(S.D,S.start,S.g0), return, end
ac = abs(S.g0);
switch f
    case 'dist'
        if fixed && strcmp(cas,'DK') && c < -10, gen(:) = true; return, end
        if fixed && any(strcmp(cas,{'DU','DT'})) && c < -10 && ((c0 > 0 && (c < -20 || abs(c)*ac > 20)) || (c < -20 && abs(c)*ac > 60)), gen(:) = true; return, end
        if stat && c0 ~= 0 && c <= -3e-3, gen(:) = true; return, end
        if ~fixed && ~stat && c0 > 0 && c <= 1
            gen(:) = true;
            if strcmp(cas,'DK')
                % (just above the tangency x* = -sqrt(2 c0) of case DK the classical line integral is the more accurate)
                xs = -sqrt(2*c0);
                gen(xs < x & x <= xs + 1e-3*max(1,abs(xs))) = false;
            end
            % (x = 0, the points of the coefficient statistic, F at x = 0 with the null value c0 + y: the classical line
            % integral, 3 to 8 times faster there, where X(v) = c0 v^2 + 1/2 has no root; as the Julia and the C rule)
            if c > -100 && any(strcmp(cas,{'DK','DU','DT'})), gen(x == 0) = false; end
            return
        end
        if ~fixed && ~stat && c0 ~= 0 && c < -100, gen(:) = true; return, end
        if ~fixed && ~stat && c0 < 0 && c <= 0 && ((strcmp(cas,'DU') && c >= -10) || (any(strcmp(cas,{'DK','DT'})) && c >= -5)), gen(:) = true; return, end
        if ~fixed && ~stat && c0 == 0 && ((any(strcmp(cas,{'DU','DT'})) && c <= 1) || (strcmp(cas,'DK') && c <= -50)), gen(:) = true; return, end
        if stat && c0 == 0 && c <= -1, gen(:) = true; return, end
        if fixed && strcmp(cas,'DK') && c0 == 0 && -1 <= c && c <= 0 && ac <= 1, gen(:) = true; return, end
    case 'icdf'
        gen = c0 == 0 && ~strcmp(cas,'DG') && c <= 1;
    case 'mom'
        gen = false;
    case 'lincdf'
        a = x;
        gen = isnan(a) || (strcmp(cas,'DK') && c <= 0 && abs(a) <= 0.1) || (strcmp(cas,'DG') && a ~= 0 && c > 3.5);
    case 'adfmax'
        gen = (stat && c <= -0.1) || (~stat && c <= -10);
    case 'envelope'
        gen = strcmp(cas,'DK') && c > -1 && c < 0;
    case 'region'
        gen = false;
end
end
