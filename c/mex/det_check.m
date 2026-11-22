% DET_CHECK  checks of the general engine (any deterministic terms) through the MEX file: designs, the general engine's
% functions, the router over the two engines, against the reference values of c/test/det_ref.txt and the classical
% engine; the route operation's ADFmax, linear form, envelope and regions; with matlab/ on the path behind this folder,
% the user-level routines of the MATLAB programs on the C router (ltratio_croute).  Run with the folder of the MEX file
% on the path (from a folder other than matlab/); about 10 seconds.
nf = 0; nc = 0;
function [nf, nc] = chk(nf, nc, got, want, tol, what)
    nc = nc + 1; d = max(abs(got(:) - want(:)));
    if ~(d <= tol), nf = nf + 1; fprintf('FAIL %s: %s vs %s (%.2e)\n', what, mat2str(got, 12), mat2str(want, 12), d);
    else, fprintf('ok   %s (%.1e)\n', what, d); end
end
DT = ltratio('design', 'trend', 1);
PC = ltratio('design', 'perron', 'C', 0.3);
PB = ltratio('design', 'perron', 'B', 0.3); PBf = ltratio('design', 'first', PB);
xs = [-4 -3 -2.2 -1.5 -0.5 0.8];
% reference values of c/test/det_ref.txt (Julia): DT at c = 0, perron C 0.3 at c = 0
[nf, nc] = chk(nf, nc, ltratio('det', 'cdf', DT, xs, 0), [0.0087334270022758325 0.13223088506362687 0.48909981798479413 0.82997874129168081 0.98365126236633393 0.99979130205429756], 1e-9, 'DT cdf c = 0');
F1 = ltratio('cdf', xs, DT, 0);                                     % the router (auto: a classical case)
F2 = ltratio('cdf', xs, 'DT', 0);                                   % the classical engine
[nf, nc] = chk(nf, nc, F1, F2, 1e-9, 'router auto against the classical engine, DT');
F3 = ltratio('cdf', xs, DT, 0, struct('engine', 'general'));
[nf, nc] = chk(nf, nc, F3, F2, 1e-9, 'router general against the classical engine, DT');
F4 = ltratio('cdf', xs, 'DT', 0, struct('engine', 'general'));
[nf, nc] = chk(nf, nc, F4, F2, 1e-9, 'case string with engine general');
q = ltratio('icdf', [0.05 0.5], PC, -5);
Fq = ltratio('det', 'cdf', PC, q, -5);
[nf, nc] = chk(nf, nc, Fq, [0.05 0.5], 1e-9, 'quantiles of Perron C (router: general)');
m = ltratio('mom', 4, DT, -5); m2 = ltratio('mom', 4, 'DT', -5);
[nf, nc] = chk(nf, nc, m, m2, 1e-8 * max(abs(m2)), 'moments DT through the router');
% the starts, GLS, c0, the statistic that detrends first
Ffix = ltratio('cdf', [-3 -1.5], DT, -5, struct('g0', 2, 'engine', 'general')); Ffc = ltratio('cdf', [-3 -1.5], 'DT', -5, struct('g0', 2));
[nf, nc] = chk(nf, nc, Ffix, Ffc, 1e-9, 'fixed start DT g0 = 2');
Fst = ltratio('det', 'cdf', ltratio('design', 'trend', 0), [-3 -1.5], -5, struct('start', 'stationary')); Fsc = ltratio('cdf', [-3 -1.5], 'SU', -5);
[nf, nc] = chk(nf, nc, Fst, Fsc, 1e-9, 'stationary start DU');
Fg = ltratio('det', 'cdf', DT, [-3 -2 -1], -5, struct('gls', -13.5)); Fgc = ltratio('cdf', [-3 -2 -1], 'DG', -5);
[nf, nc] = chk(nf, nc, Fg, Fgc, 1e-9, 'GLS detrending DT at -13.5 against case DG');
Ff = ltratio('det', 'cdf', PBf, [-4.5 -3.8 -3 -2], 0);
[nf, nc] = chk(nf, nc, Ff, [0.0053667464914887832 0.045349287493156108 0.25097442724109131 0.74622199188979299], 1e-9, 'detrended first, Perron B 0.3');
% companions
P = ltratio('det', 'lincdf', ltratio('design', 'trend', -1), -0.2639, -3.5, 0);
[nf, nc] = chk(nf, nc, P, ltratio('lincdf', -0.2639, -3.5, 'DK', 0), 1e-9, 'linear form DK');
Q = ltratio('det', 'pticdf', DT, 0.05, -13.5);
[nf, nc] = chk(nf, nc, ltratio('det', 'ptcdf', DT, Q, -13.5, 0), 0.05, 1e-9, 'point-optimal quantile');
qlm = ltratio('det', 'lmicdf', DT, 0.05, 0);
[nf, nc] = chk(nf, nc, ltratio('det', 'lmcdf', DT, qlm, 0), 0.05, 1e-9, 'LM quantile');
% a design from a function and from a matrix
Df = ltratio('design', 'fun', @(r) [1; r; double(r > 0.4); max(r - 0.4, 0)], 4, 0.4);
[nf, nc] = chk(nf, nc, ltratio('det', 'cdf', Df, [-3.5 -2.5], -5), ltratio('det', 'cdf', ltratio('design', 'perron', 'C', 0.4), [-3.5 -2.5], -5), 1e-10, 'fdesign against perron C 0.4');
t = (1:500)'/500; X = [ones(500,1) t (t > 0.3) max(t - 0.3, 0)];
Dm = ltratio('design', 'matrix', X);
[nf, nc] = chk(nf, nc, ltratio('det', 'cdf', Dm, [-3.5 -2.5], 0), ltratio('det', 'cdf', PC, [-3.5 -2.5], 0), 1e-10, 'design from a matrix against perron C 0.3');
% the route operation's ADFmax, linear form, envelope and regions, and the warning of a brittle regime
x3 = [-3 -2.4 -1.5]; D0 = ltratio('design', 'trend', 0);
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'adfcdf', 0, [2 1 1.5], x3, -5), ltratio('det', 'adfcdf', D0, x3, -5, struct('g0', 1.5)), 1e-12, 'ADFmax DU with a fixed start: the general engine');
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'adfcdf', 0, [2 0 0], x3, 0), ltratio_mex('adfmax', 2, x3, 2, 0, 0), 1e-9, 'ADFmax DU, zero start (router auto)');
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'lincdf', 0, [1 0 0], [-0.6 -0.2639], 0, -3.5, NaN), ltratio('lincdf', [-0.6 -0.2639], -3.5, 'DK', 0), 1e-9, 'the linear form through the router');
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'envelope', 0, [1 0 0], [-5 -10], 0.05), ltratio('envelope', [-5 -10], 'DK'), 1e-8, 'the envelope through the router');
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'region', 0, [1 0 0], 1, -5, -3.5, -0.2639), ltratio_mex('lincdf', -0.2639, -3.5, 1, -5), 1e-9, 'the region below a line through the router');
[nf, nc] = chk(nf, nc, ltratio_mex('route', 'region', 0, [1 1 1.5], 2, -5, -2, 0), ltratio('det', 'cdf', ltratio('design', 'trend', -1), -2, -5, struct('g0', 1.5)), 1e-10, 'the region of case DK with a fixed start: the general engine');
lastwarn('');
ws = warning('off', 'ltratio:brittle'); ltratio_mex('route', 'cdf', 0, [2 1 0.5], -5, -20, 2); warning(ws);
[~, wid] = lastwarn; [nf, nc] = chk(nf, nc, double(strcmp(wid, 'ltratio:brittle')), 1, 0, 'the warning ltratio:brittle of the router (DU, fixed start, c0 > 0 at c = -20)');
% the general engine's method: the line integral with the explicit kernel ('line') and the inversion on paths ('paths')
T2 = ltratio('design', 'trend', 2); xt = [-4 -3.2 -2.5 -1.8];
Fl = ltratio('det', 'cdf', T2, xt, -5, struct('method', 'line')); Fp = ltratio('det', 'cdf', T2, xt, -5, struct('method', 'paths'));
[nf, nc] = chk(nf, nc, Fl, Fp, 1e-10, 'quadratic trend: the line integral against the paths');
[nf, nc] = chk(nf, nc, ltratio('det', 'cdf', T2, xt, -5), Fl, 0, 'quadratic trend at c0 = 0: auto is the line integral');
[nf, nc] = chk(nf, nc, ltratio('det', 'cdf', T2, xt, -50), ltratio('det', 'cdf', T2, xt, -50, struct('method', 'paths')), 0, 'quadratic trend at c = -50: auto is the paths');
[nf, nc] = chk(nf, nc, ltratio('cdf', xt, T2, -5, struct('method', 'paths')), Fp, 0, 'the router with a design and opts.method paths');
[nf, nc] = chk(nf, nc, ltratio('cdf', xs, DT, -5, struct('method', 'line')), ltratio('cdf', xs, 'DT', -5, struct('method', 'single')), 0, 'a classical design with line: the classical line integral');
[nf, nc] = chk(nf, nc, ltratio('icdf', 0.05, T2, -5, struct('method', 'line')), ltratio('icdf', 0.05, T2, -5, struct('method', 'paths')), 1e-9, 'quantile of the quadratic trend by both methods');
msg = ''; try, ltratio('det', 'cdf', T2, xt, -5, struct('method', 'line', 'gls', -13.5)); catch e, msg = e.message; end
[nf, nc] = chk(nf, nc, double(contains(msg, 'GLS detrending is not covered')), 1, 0, 'method line refuses GLS detrending');
% with this folder ahead of matlab/ on the path: the user-level routines on the C router
if exist('ltratio_croute', 'file') == 2 && ltratio_croute('on')
    [nf, nc] = chk(nf, nc, ltratio_cdf(xs, 'DT', -5), ltratio('cdf', xs, 'DT', -5), 1e-9, 'ltratio_cdf on the C router');
    [nf, nc] = chk(nf, nc, ltratio_cdf([-3.5 -2.5], dfdet_perron('C', 0.3), -5), ltratio('det', 'cdf', PC, [-3.5 -2.5], -5), 1e-12, 'ltratio_cdf with a design on the C router');
    [nf, nc] = chk(nf, nc, adfmax_cdf(x3, 'DU', -5, struct('g0', 1.5)), ltratio('det', 'adfcdf', D0, x3, -5, struct('g0', 1.5)), 1e-12, 'adfmax_cdf with a fixed start on the C router');
    [nf, nc] = chk(nf, nc, ltratio_lincdf(-0.2639, -3.5, dfdet_trend(-1), -5, struct('g0', 1.5)), ltratio('det', 'lincdf', ltratio('design', 'trend', -1), -0.2639, -3.5, -5, struct('g0', 1.5)), 1e-12, 'ltratio_lincdf of a design with a fixed start on the C router');
    [nf, nc] = chk(nf, nc, power_envelope([-5 -10], 'DG'), ltratio('envelope', [-5 -10], 'DG'), 1e-8, 'power_envelope on the C router');
    id = ''; try, ltratio_cdf(-20, 'DU', -551); catch e, id = e.identifier; end
    [nf, nc] = chk(nf, nc, double(strcmp(id, 'ltratio:args')), 1, 0, 'the checks of the MATLAB routines before the C router');
end
fprintf('%d checks, %d failed\n', nc, nf);
if nf, exit(1); end
