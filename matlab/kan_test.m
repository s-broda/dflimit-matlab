function nfail = kan_test(quick)
%KAN_TEST  The devices of Bao and Kan in these programs against the values of DFLimit.jl.
%
%   nfail = KAN_TEST compares KAN_CORE (the closed-form transform of the six classical cases at points and in Taylor-mode
%   arithmetic, the scale, the single inversion of the coefficient statistic, of the linear form and of the power
%   envelope, the moments from jets, the closed form at zero), the user-level routines with their defaults (COEF_CDF,
%   COEF_PDF, COEF_ICDF, COEF_MOM, LTRATIO_MOM, LTRATIO_LINCDF, POWER_ENVELOPE), the general engine's transform with the
%   branch from the pivots and the point-optimal statistic by the single inversion (DFDET_KAN, PT_CDF, PT_PDF,
%   PT_QUANTILE, PT_POWER, PT_ENVELOPE) and the branch at the point of the classical engine (DFX_CORE('tparts')) with the
%   values of the Julia package in validation/kan_ref.txt (the file c/test/kan_ref.txt of the C engine's test, written by
%   c/test/kan_ref.jl), each kind of value with its tolerance, and prints the largest difference per kind.  KAN_TEST(true)
%   leaves out the moments from jets and half of the single inversions.  With the folder of the MEX file ahead on the
%   path the user-level routines run on the C router (the lines of DFX_CORE('tparts') are then left out).  About 30
%   seconds either way.
%
%   See also KAN_CORE, DFDET_KAN, LTRATIO_TEST, DFDET_TEST.
if nargin < 1, quick = false; end
fn = fullfile(fileparts(mfilename('fullpath')), 'validation', 'kan_ref.txt');
L = strsplit(fileread(fn), newline);
tol = struct('lcf',1e-12,'jet',1e-11,'scale',1e-12,'coef',1e-11,'coefpdf',1e-11,'coefq',1e-8,'lin',1e-11,'env',1e-11,'tmom',1e-10,'cmom',1e-10, ...
    'ray',1e-12,'rcoef',1e-11,'rcoefpdf',1e-11,'rcoefq',1e-8,'rcoefmom',1e-10,'rmom',1e-10,'rlin',1e-11,'renv',1e-11,'dcf',1e-10, ...
    'ptcdf',1e-11,'ptpdf',1e-11,'ptq',1e-10,'ptpow',1e-10,'ptenv',1e-11,'tpp',1e-10,'tppfar',2e-8);
kinds = fieldnames(tol); W = struct(); N = struct(); NF = struct();
for k = 1:numel(kinds), W.(kinds{k}) = 0; N.(kinds{k}) = 0; NF.(kinds{k}) = 0; end
nfail = 0; nskip = 0; t0 = tic; ws = warning('off','ltratio:inversion');
cl = onCleanup(@() warning(ws));
for i = 1:numel(L)
    l = strtrim(L{i}); if isempty(l), continue, end
    p = strsplit(l); tag = p{1}; v = str2double(p(2:end));
    try
        [d, kind] = one(tag, v, quick);
    catch ME
        d = Inf; kind = tag; fprintf('ERROR %s: %s\n', l(1:min(end,100)), ME.message);
    end
    if isempty(kind), nskip = nskip + 1; continue, end
    if ~isfield(tol, kind), nskip = nskip + 1; continue, end
    N.(kind) = N.(kind) + numel(d); W.(kind) = max([W.(kind), d(:)']);
    bad = ~(d <= tol.(kind));
    if any(bad), NF.(kind) = NF.(kind) + nnz(bad); nfail = nfail + nnz(bad); fprintf('FAIL %s: difference %.3g (tolerance %.1g)\n', l(1:min(end,120)), max(d), tol.(kind)); end
end
for k = 1:numel(kinds)
    if N.(kinds{k}), fprintf('  %-9s %5d checks, %d failed, worst difference %.2e (tolerance %.0e)\n', kinds{k}, N.(kinds{k}), NF.(kinds{k}), W.(kinds{k}), tol.(kinds{k})); end
end
fprintf('kan_test: %d checks, %d failed, %d lines left out (%.0f s)\n', sum(cellfun(@(k) N.(k), kinds)), nfail, nskip, toc(t0));
end

function [d, kind] = one(tag, v, quick)
d = []; kind = tag;
persistent cnt
if isempty(cnt), cnt = 0; end
cnt = cnt + 1; skipq = quick && mod(cnt, 2) == 0;
cplx = @(k) complex(v(k), v(k + 1));
switch tag
    case 'lcf'
        F = kan_core('limcf', v(1), v(2), v(3), v(4));
        Lc = kan_core('logcf', F, cplx(5), cplx(7)); Lr = cplx(9);
        d = abs(Lc - Lr) / max(1, abs(Lr));
    case 'jet'
        F = kan_core('limcf', v(1), v(2), v(3), v(4)); n = v(5); o = 5;
        a = complex(v(o + 1:2:o + 2 * n), v(o + 2:2:o + 2 * n)).'; o = o + 2 * n;
        b = complex(v(o + 1:2:o + 2 * n), v(o + 2:2:o + 2 * n)).'; o = o + 2 * n;
        r = complex(v(o + 1:2:o + 2 * n), v(o + 2:2:o + 2 * n)).';
        Lc = kan_core('logcf', F, a, b);
        d = max(abs(Lc - r) ./ max(1, abs(r)));
    case 'scale'
        F = kan_core('limcf', v(1), v(2), v(3), v(4)); s = kan_core('scale', F);
        d = max(abs(s - v(5:8)) ./ max(1, abs(v(5:8))));
    case {'coef', 'coefpdf'}
        if skipq, kind = ''; return, end
        F = kan_core('limcf', v(1), v(2), v(3), v(4));
        if strcmp(tag, 'coef'), P = kan_core('coef_cdf1', F, v(6), v(5), 1e-12); d = abs(P - v(7));
        else, f = kan_core('coef_pdf1', F, v(6), v(5), 1e-12); d = abs(f - v(7)) / max(1, abs(v(7))); end
    case 'coefq'
        if skipq, kind = ''; return, end
        F = kan_core('limcf', v(1), v(2), v(3), v(4));
        y = kan_core('coef_icdf1', F, v(6), v(5), 1e-9, 1e-12); d = abs(y - v(7)) / max(1, abs(v(7)));
    case 'lin'
        if skipq, kind = ''; return, end
        F = kan_core('limcf', v(1), v(2), v(3), v(4));
        P = kan_core('lin_cdf1', F, v(6), v(5), 1e-12); d = abs(P - v(7));
    case 'env'
        P = kan_core('envelope1', v(2), v(1), 1e-12); d = abs(P - v(3));
    case {'tmom', 'cmom'}
        if quick, kind = ''; return, end
        F = kan_core('limcf', v(1), v(2), v(3), v(4)); k = v(6); pw = 1; if tag(1) == 't', pw = 0.5; end
        m = kan_core('ratio_mom', F, k, v(5) - v(3), pw, 1e-10); r = v(7:6 + k);
        sc = max(abs(r), abs(r(min(2, k))).^((1:k) / 2)); d = max(abs(m - r) ./ sc);
    case 'ray'
        P = kan_core('ray_cdf', v(1), v(2), v(3), v(4)); d = abs(P - v(5));
    case {'rcoef', 'rcoefpdf', 'rcoefq', 'rcoefmom', 'rmom'}
        cas = casestr(v(1), v(2)); o = struct('c0', v(5)); if v(2) == 1, o.g0 = v(3); end
        c = v(4); n = v(6);
        if any(strcmp(tag, {'rcoefmom', 'rmom'})), x = []; r = v(7:6 + n); else, x = v(7:6 + n); r = v(7 + n:6 + 2 * n); end
        switch tag
            case 'rcoef', g = coef_cdf(x, cas, c, o);
            case 'rcoefpdf', g = coef_pdf(x, cas, c, o);
            case 'rcoefq', o.tol = 1e-9; g = coef_icdf(x, cas, c, o);
            case 'rcoefmom', g = coef_mom(n, cas, c, o);
            case 'rmom', g = ltratio_mom(n, cas, c, o);
        end
        if any(strcmp(tag, {'rcoefmom', 'rmom'})), d = abs(g - r) ./ max(abs(r), abs(r(min(2, n))).^((1:n) / 2));
        elseif strcmp(tag, 'rcoef'), d = abs(g - r); d(isnan(g) & isnan(r)) = 0;
        else, d = abs(g - r) ./ max(1, abs(r)); end
    case 'rlin'
        st = {'zero', 'fixed', 'stationary'}; D = dfdet_trend(v(1)); o = struct('start', st{v(2) + 1});
        if v(2) == 1, o.g0 = v(3); end
        n = v(6); b = v(7:6 + n); r = v(7 + n:6 + 2 * n);
        g = ltratio_lincdf(b, v(5), D, v(4), o); d = abs(g - r);
    case 'renv'
        n = v(2); c = v(3:2 + n); r = v(3 + n:2 + 2 * n);
        g = power_envelope(c, 'DK', struct('alpha', v(1))); d = abs(g - r);
    case {'dcf', 'ptcdf', 'ptpdf', 'ptq', 'ptpow', 'ptenv'}
        [D, k] = design(v); w = v(k:end);
        switch tag
            case 'dcf'
                st = {'zero', 'fixed', 'stationary'}; gls = []; if ~isnan(w(2)), gls = w(2); end
                Mo = dfdet_startmodel(D, w(3), st{w(1) + 1}, 0, gls);
                Lc = dfdet_kan('logcf', dfdet_kan('cf', Mo, w(3)), w(4), w(5)); Lr = complex(w(6), w(7));
                d = abs(Lc - Lr) / max(1, abs(Lr));
            case {'ptcdf', 'ptpdf'}
                n = w(3); q = w(4:3 + n); r = w(4 + n:3 + 2 * n);
                if strcmp(tag, 'ptcdf'), g = pt_cdf(q, D, w(1), w(2)); d = abs(g - r);
                else, g = pt_pdf(q, D, w(1), w(2)); d = abs(g - r) ./ max(1, abs(r)); end
            case 'ptq'
                g = pt_quantile(w(2), D, w(1)); d = abs(g - w(3)) / max(1, abs(w(3)));
            case 'ptpow'
                [b, P] = pt_power(w(2), D, w(1), struct('alpha', w(3))); d = max(abs(b - w(4)) / max(1, abs(w(4))), abs(P - w(5)));
            case 'ptenv'
                n = w(2); c = w(3:2 + n); r = w(3 + n:2 + 2 * n);
                g = pt_envelope(c, D, struct('alpha', w(1))); d = abs(g - r); d(isnan(g) & isnan(r)) = 0;
        end
    case 'tpp'
        % (with the folder of the MEX file ahead DFX_CORE is the C engine's drop-in, which has no tparts: dfl_test_kan checks
        % the C engine's)
        if ltratio_croute('on'), kind = ''; return, end
        m = v(1); st = logical(v(2)); c = v(3); t = cplx(4);
        [P, s1, s2, M0] = dfx_core('tparts', m, c, t, st);
        R = [cplx(6), cplx(8), cplx(10), cplx(12)]; Q = [P, s1, s2, M0];
        nn = 4; if (m == 1 || m == 4) && ~st, nn = 2; end
        d = max(abs(Q(1:nn) - R(1:nn)) ./ abs(R(1:nn)));
        if c < -20, kind = 'tppfar'; end
    otherwise
        kind = '';                                           % (gp: coef_cdf_gp has no counterpart here; dline: the MATLAB
                                                             % programs have no line integral for designs)
end
end

function cas = casestr(model, start)
if model == 4, cas = 'DG'; return, end
k = 'KUT'; if start == 2, cas = ['S' k(model)]; else, cas = ['D' k(model)]; end
end

function [D, k] = design(v)
% the design of a line: D kind np par.. (kind 0 trend(par), 1, 2, 3 Perron's models A, B, C), then the arguments from v(k)
% (the 'D' of the line is NaN in v(1))
kind = v(2); np = v(3); par = v(4:3 + np); k = 4 + np;
if kind == 0, D = dfdet_trend(par(1)); else, mdl = 'ABC'; D = dfdet_perron(mdl(kind), par); end
end
