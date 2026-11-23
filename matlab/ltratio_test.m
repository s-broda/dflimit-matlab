function nfail = ltratio_test(quick,strict)
%LTRATIO_TEST  Validation suite for the limiting t-ratio programs.
%
%   nfail = LTRATIO_TEST runs the checks below, prints each comparison with
%   its tolerance and a PASS/FAIL mark, and returns the number of failures.
%   LTRATIO_TEST(true) runs the quick version (fewer points; about a minute
%   in MATLAB), LTRATIO_TEST(quick,true) raises an error at the end when a
%   check has failed, for use in scripts.
%
%    1. the two representations against each other: the line integral
%       (DFX_CORE) and the conditional integral (TR_LIMDIST), distribution function
%       and density, zero and stationary starts, c on both sides of the unit
%       root and local null values c0 ~= 0 (the two share nothing but the
%       Fredholm determinant);
%    2. exact values: P[tau <= 0] = erf((2v)^(-1/2)) in case DK (c0 = 0),
%       P[tau_mu <= 0] from the distribution of the coefficient statistic
%       (Bao and Kan, 2026), and P[A' - c0 B <= 0] at c = 10, c0 = 5 from a
%       one-dimensional inversion; the 5% points against nine-digit values and
%       the asymptotic critical values of MacKinnon (2010) to their accuracy;
%       the quantile function LTRATIO_ICDF;
%    3. the initial observation: continuity at a zero start, the stationary
%       start as the Gauss-Hermite average of the fixed start (an identity the
%       programs do not use), and the fixed start of case DK against the
%       closed form of the tilted mean (TR_JBSTART);
%    4. moments: Nabeya (1999) at the unit root (eight digits), LTRATIO_MOM
%       against DFX_MOMENTS (independent implementations), and the mean
%       against a Gauss-Legendre integral of the density;
%    5. the density against central differences of the distribution function;
%    6. argument validation and non-finite x;
%    7. the tails with c0 ~= 0 (the far roots of X(v), the tail rule, no
%       negative density);
%    8. far from the unit root: c = -100 and -480 against the Julia engine's
%       line integral (julia/largec_ref.txt), both methods, with the moments
%       and a quantile;
%    9. the double zero x = -sqrt(2 c0) of a(b) (c0 > 0): finite values in
%       the six cases (the density of case DK is infinite there);
%   10. the applications against the Julia engine julia/DFExact.jl: the
%       ADFmax statistic, the coefficient statistic, the region P{U <= g(V)}
%       (a line, the boundary of the weighted-average-power test, the
%       parabola of the t-ratio), the Gil-Pelaez inversion of U - aV, the
%       power envelope;
%   11. the contour integral (TR_CONTOUR, a fixed start in cases DU, DT)
%       against the reference values of the Julia engine
%       (c/test/contour_ref.txt), the automatic choice there, and the error
%       outside its scope;
%   12. the default (one evaluation), verified and the warnings of the
%       brittle regimes;
%   13. the upper tail above the unit root (November 23, 2026): case DK with
%       c0 ~= c beyond the line integral's bound, and for c > 4 at x > 300,
%       the conditional integral with its checks (before, its flagged points
%       got one evaluation, off by up to 2.4), the refusal above c = 14, the
%       general engine's line cut to the rate with c0, and the final check of
%       a distribution function (LTRATIO_RANGECHECK), against an independent
%       computation (the law of V given J(1) by Girsanov's theorem, inverted
%       and integrated over J(1)).
%
%   See also LTRATIO_CDF, LTRATIO_PDF, LTRATIO_ICDF, LTRATIO_MOM, ADFMAX_CDF,
%   COEF_CDF, LTRATIO_REGION, POWER_ENVELOPE.
if nargin < 1 || isempty(quick), quick = false; end
if nargin < 2 || isempty(strict), strict = false; end
nfail = 0; ntest = 0; t00 = tic;
S = struct('method','single'); D = struct('method','double','verified',true);   % (verified: the checked computation)

fprintf('\n=== 1. line integral (DFX) against conditional integral (TR_LIMDIST) ===\n');
if quick
    P = {'DK',0,0,[-2.5 -0.5 1.5]; 'DU',0,0,[-3.4 -1.6 0.2]; 'DT',-5,0,[-3.5 -2.5]; 'SU',-5,0,[-2.8 -1.2]; 'DK',2,2,[-1.5 0 1.5]};
else
    P = {'DK',0,0,[-3 -2.5 -1.94 -0.5 0.9 2]; 'DU',0,0,[-3.43 -2.86 -1.57 -0.08 0.6]; 'DT',0,0,[-3.96 -3.41 -2.18 -0.94];
         'DK',-5,0,[-3 -1.7 -0.5]; 'DU',-5,0,[-3 -2.1 -1.2]; 'DT',-10,0,[-3.5 -2.86 -2.2];
         'SK',-5,0,[-2.5 -1.8 -1]; 'SU',-5,0,[-3 -2.1 -1.2]; 'ST',-10,0,[-3.5 -2.9 -2.2];
         'DK',2,2,[-1.5 0 1.5 3]; 'DU',2,1,[-3 -0.5 1.3]; 'DT',-2.5,-5,[-3.6 -1.1 0.3]; 'DU',5,-5,[57.11 126.85]; 'DK',0,-5,[-1.1 2.7 6.5]};
end
for k = 1:size(P,1)
    cas = P{k,1}; c = P{k,2}; c0 = P{k,3}; xs = P{k,4};
    o = struct('c0',c0); ov = struct('c0',c0,'verified',true);      % (the representations against each other: verified)
    [F1,i1] = ltratio_cdf(xs,cas,c,setfield(o,'method','single'));
    [F2,i2] = ltratio_cdf(xs,cas,c,setfield(ov,'method','double'));
    f1 = ltratio_pdf(xs,cas,c,setfield(o,'method','single'));
    f2 = ltratio_pdf(xs,cas,c,setfield(ov,'method','double'));
    check(sprintf('%s c=%g c0=%g: F single vs double (%d points; %.1f/%.1f s)',cas,c,c0,numel(xs),i1.time,i2.time),max(abs(F1-F2)),1e-11);
    check(sprintf('%s c=%g c0=%g: f single vs double (relative)',cas,c,c0),max(abs(f1-f2)./max(abs(f2),1e-3)),1e-10);
    [Fa,ia] = ltratio_cdf(xs,cas,c,ov); fa = ltratio_pdf(xs,cas,c,ov);        % the automatic choice (ltratio_method), verified
    check(sprintf('%s c=%g c0=%g: F auto, verified (%d of %d by the double) vs double (%.1f s)',cas,c,c0,ia.ndouble,numel(xs),ia.time),max(abs(Fa-F2)),1e-11);
    check(sprintf('%s c=%g c0=%g: f auto, verified vs double (relative)',cas,c,c0),max(abs(fa-f2)./max(abs(f2),1e-3)),1e-10);
    ws = warning('off','ltratio:brittle');
    [F0,i0] = ltratio_cdf(xs,cas,c,o); f0 = ltratio_pdf(xs,cas,c,o);          % the default: one evaluation
    warning(ws);
    check(sprintf('%s c=%g c0=%g: F, f auto, one evaluation (%d of %d by the double) vs double, verified (%.2f s)',cas,c,c0,i0.ndouble,numel(xs),i0.time), ...
          max([abs(F0-F2) abs(f0-f2)./max(abs(f2),1e-3)]),1e-8);
end

fprintf('\n=== 2. exact values, critical values and quantiles ===\n');
for cc = [-10 0 2 5 20]
    v = (exp(2*cc) - 1)/(2*cc); if cc == 0, v = 1; end
    F = ltratio_cdf(0,'DK',cc,S); G = ltratio_cdf(0,'DK',cc,D);
    check(sprintf('DK c=%g: F(0) = erf((2v)^(-1/2)), single',cc),abs(F-erf(1/sqrt(2*v))),1e-12);
    check(sprintf('DK c=%g: F(0) = erf((2v)^(-1/2)), double',cc),abs(G-erf(1/sqrt(2*v))),1e-12);
end
check('DU c=0: F(0) = 0.95756559 (Bao and Kan, 2026, eight digits)',abs(ltratio_cdf(0,'DU',0)-0.95756559),5e-9);
% the GLS-detrended trend case (DG): P(tau <= 0) at c = 0 in closed form, the 5% critical value, the coefficient form and
% the power at c = -13.5 against the Julia implementation (10 digits)
lamg = (1 + 13.5)/(1 + 13.5 + 13.5^2/3);
check('DG c=0: F(0) = erf((0.4)^(-1/2) (1-lambda)^(-1)) (closed form)',abs(ltratio_cdf(0,'DG',0)-erf(1/sqrt(0.4*(1-lamg)^2))),1e-11);
check('DG c=0: F(-2.8456308262) = 0.05 (the 5% critical value)',abs(ltratio_cdf(-2.8456308262,'DG',0)-0.05),1e-10);
check('DG c=-13.5: F(-2.89) = 0.4817352773 (Julia)',abs(ltratio_cdf(-2.89,'DG',-13.5)-0.4817352773219187),1e-10);
check('DG c=-5 c0=-8: F(0) = 0.5677921222 (the coefficient form, Julia)',abs(ltratio_cdf(0,'DG',-5,struct('c0',-8))-0.567792122181),1e-10);
check('DG: LTRATIO_ICDF(0.05) = -2.8456308262',abs(ltratio_icdf(0.05,'DG',0)+2.8456308262),1e-8);
ref = [0.000302655553129 0.000414190811913]; cs = {'DK','DU','DT'};
for j = 2:3
    F = ltratio_cdf(0,cs{j},10,struct('c0',5));
    check(sprintf('%s c=10 c0=5: F(0) vs a one-dimensional inversion',cs{j}),abs(F-ref(j-1)),1e-11);
end
q5 = [-1.94087325 -2.86156886 -3.41045363];           % nine-digit values (independent route)
mck = [-2.56574 -1.94100 -1.61682; -3.43035 -2.86154 -2.56677; -3.95877 -3.41049 -3.12705];   % MacKinnon (2010), 1%, 5%, 10%
for j = 1:3
    F = ltratio_cdf(q5(j),cs{j});
    check(sprintf('%s: F at the 5%% point',cs{j}),abs(F-0.05),2e-9);
    Fm = ltratio_cdf(mck(j,:),cs{j});
    check(sprintf('%s: F at MacKinnon''s critical values vs 0.01, 0.05, 0.10',cs{j}),max(abs(Fm-[0.01 0.05 0.10])),1.5e-4);
    [q,inf1] = ltratio_icdf([0.05 0.5],cs{j});
    check(sprintf('%s: LTRATIO_ICDF(0.05) vs the nine-digit value',cs{j}),abs(q(1)-q5(j)),6e-9);
    check(sprintf('%s: F(LTRATIO_ICDF(p)) - p, convergence',cs{j}),max(abs(inf1.F-[0.05 0.5])) + ~all(inf1.conv),1e-9);
end
[q,inf1] = ltratio_icdf([0.01 0.05 0.5],'SU',-5);
check('SU c=-5: quantiles increasing and converged',double(~all(diff(q) > 0)) + ~all(inf1.conv),0.5);

fprintf('\n=== 3. the initial observation ===\n');
cc = -5; x0 = -2;
for cas = cs
    F0 = ltratio_cdf(x0,cas{1},cc); F1 = ltratio_cdf(x0,cas{1},cc,struct('g0',1e-6));
    check(sprintf('%s c=%g: start 1e-6 vs zero start',cas{1},cc),abs(F1-F0),1e-9);
end
[gz,gw] = ghnodes(24); s0 = sqrt(-1/(2*cc)); ip = find(gz > 0);   % (16 nodes leave 7e-11)
sts = {'SK','DK'; 'SU','DU'}; if ~quick, sts = [sts; {'ST','DT'}]; end
for j = 1:size(sts,1)
    Fs = ltratio_cdf(x0,sts{j,1},cc);
    Fm = 0;
    for k = ip(:).', Fm = Fm + 2*gw(k)*tr_limdist(x0,sts{j,2},cc,struct('g0',s0*gz(k),'c0',0)); end   % (the engine: nodes beyond the region of LTRATIO_ARGS)
    check(sprintf('%s c=%g: stationary start vs the Gauss-Hermite average of %s over J(0)',sts{j,1},cc,sts{j,2}),abs(Fs-Fm),2e-11);
end
% fixed start in case DK: U = (J(1)^2 - J(0)^2)/2 - 1/2, so P[tau <= 0] = P[J(1)^2 <= 1 + g0^2]
% with J(1) ~ N(g0 e^c, (e^(2c)-1)/(2c))
g0 = 1.5; v = (exp(2*cc) - 1)/(2*cc); m = g0*exp(cc); r = sqrt(1 + g0^2);
F = ltratio_cdf(0,'DK',cc,struct('g0',g0));
Pex = (erf((r - m)/sqrt(2*v)) + erf((r + m)/sqrt(2*v)))/2;
check(sprintf('DK c=%g g0=%g: F(0) = P[J(1)^2 <= 1 + g0^2]',cc,g0),abs(F-Pex),1e-11);
% fixed start in cases DU, DT: the closed-form initial-value term Q(lam,s) of TR_NSSTART
% is, as a function of lam, Q(0,s) + (lam/2) [m1^2/(1 - 2 lam kap_1) - m2^2/(1 - 2 lam kap_2)]
% with m_k^2 >= 0 at real s (a noncentral form of rank two): fit m1^2, m2^2 at two values
% of lam and predict a third
for cas = {'DU','DT'}
    s = -1; [~,~,kap] = tr_cflim(0,-1i*s,cas{1},cc,'B'); k1 = kap(1); k2 = kap(2);
    la = [0.1 0.2 0.3]; Q = tr_nsstart(la,s,cas{1},cc) - tr_nsstart(0,s,cas{1},cc);
    A = [la(1)/2/(1-2*la(1)*k1), -la(1)/2/(1-2*la(1)*k2); la(2)/2/(1-2*la(2)*k1), -la(2)/2/(1-2*la(2)*k2)];
    mm = A\Q(1:2).';
    pred = la(3)/2*(mm(1)/(1-2*la(3)*k1) - mm(2)/(1-2*la(3)*k2));
    check(sprintf('%s c=%g: TR_NSSTART is a rank-two noncentral form in lam (prediction)',cas{1},cc),abs(pred-Q(3)) + any(mm < 0) + max(abs(imag(mm))),1e-12);
end

fprintf('\n=== 4. moments ===\n');
nab = [-0.42309565 1.14159508 -1.06171218 3.52380298;
       -1.53296244 3.05599535 -6.71991749 16.34553684;
       -2.18135582 5.32067592 -14.03580010 39.55443050];   % Nabeya (1999, Table 4)
for j = 1:3
    m = ltratio_mom(4,cs{j});
    check(sprintf('%s: moments 1-4 vs Nabeya (1999)',cs{j}),max(abs(m-nab(j,:))),1e-7);
end
for cfg = {{'DU',0},{'DU',-5},{'SU',-5},{'DT',-10},{'SK',-2.5}}
    cas = cfg{1}{1}; cc = cfg{1}{2};
    m1 = ltratio_mom(4,cas,cc);
    model = find('KUT' == cas(2)); o = struct(); if cas(1) == 'S', o.start = 'stationary'; end
    m2 = dfx_moments(model,cc,o);
    check(sprintf('%s c=%g: LTRATIO_MOM vs DFX_MOMENTS (relative)',cas,cc),max(abs(m1(:)-m2(:))./abs(m2(:))),1e-9);
end
% the mean against a Gauss-Legendre integral of the density (DU, c = -5, c0 = 0; 192
% nodes: with 96 the rule itself is off by 8e-7, the single and the conditional integral
% agreeing on both)
m = ltratio_mom(1,'DU',-5,struct('c0',0));
[xg,wg] = glnodes(192,-12,6);
f = ltratio_pdf(xg,'DU',-5,struct('c0',0,'method','double'));
check('DU c=-5 c0=0: mean vs integral of x f(x) over [-12, 6]',abs(sum(wg.*xg.*f)-m),1e-9);
check('DU c=-5 c0=0: mass of the density over [-12, 6]',abs(sum(wg.*f)-1),1e-9);

fprintf('\n=== 5. density against the distribution function ===\n');
h = 1e-3;
for cfg = {{'DK',0,0,-1.5},{'DU',-5,0,-2.2},{'ST',-10,0,-3},{'DT',2,1,-2}}
    cas = cfg{1}{1}; cc = cfg{1}{2}; c0 = cfg{1}{3}; x = cfg{1}{4}; o = struct('c0',c0);
    f = ltratio_pdf(x,cas,cc,o);
    Fp = ltratio_cdf(x + [-2 -1 1 2]*h,cas,cc,o);
    fd = (Fp(1) - 8*Fp(2) + 8*Fp(3) - Fp(4))/(12*h);
    check(sprintf('%s c=%g c0=%g: f vs fifth-order difference of F at x=%g',cas,cc,c0,x),abs(f-fd),1e-9);
end

fprintf('\n=== 6. argument validation and non-finite x ===\n');
bad = {{-2,'XX',0,struct()},{-2,'SU',0,struct()},{-2,'DU','a',struct()},{-2,'DU',0,struct('c0',[1 2])}, ...
       {-2,'SU',-5,struct('g0',1)},{-2,'DU',0,struct('method','triple')},{[1 2]*1i,'DU',0,struct()}};
nerr = 0;
for k = 1:numel(bad)
    try, ltratio_cdf(bad{k}{:}); catch, nerr = nerr + 1; end
end
check('invalid arguments rejected',double(nerr ~= numel(bad)),0);
F = ltratio_cdf([-Inf Inf NaN],'DK',0);
check('F(-Inf), F(Inf), F(NaN) = 0, 1, NaN',abs(F(1)) + abs(F(2)-1) + ~isnan(F(3)),0);
f = ltratio_pdf([-Inf Inf NaN],'DU',0,D);
check('f(-Inf), f(Inf), f(NaN) = 0, 0, NaN (double)',abs(f(1)) + abs(f(2)) + ~isnan(f(3)),0);
f = ltratio_pdf(-sqrt(4),'DK',0,struct('c0',2));
check('DK c0=2: f(-sqrt(2 c0)) = Inf',double(~isinf(f)),0);
% the range of c: below -550 an error (the same message in the C and Julia engines)
nerr = 0; ids = {};
calls = {@() ltratio_cdf(-20,'DU',-551), @() ltratio_pdf(-20,'SK',-700), @() ltratio_mom(2,'DT',-551), @() ltratio_icdf(0.05,'DG',-551)};
for k = 1:numel(calls)
    try, calls{k}(); catch e, if strcmp(e.identifier,'ltratio:args') && contains(e.message,'c must be at least -550'), nerr = nerr + 1; end, end
end
check('c below -550 rejected (cdf, pdf, mom, icdf)',double(nerr ~= numel(calls)),0);
% the region of a fixed start (an error outside it in each case, and for c < -20 in these programs), values inside
nerr = 0;
calls = {@() ltratio_cdf(-2,'DK',0,struct('g0',2.5)), @() ltratio_pdf(-2,'DU',-10,struct('g0',2)), @() ltratio_mom(2,'DT',-110,struct('g0',0.5)), @() ltratio_cdf(-5,'DK',-30,struct('g0',1)), ...
         @() ltratio_cdf(-5,'DU',-30,struct('g0',1,'method','double'))};
for k = 1:numel(calls)
    try, calls{k}(); catch e, if strcmp(e.identifier,'ltratio:args') && contains(e.message,'with a fixed start g0'), nerr = nerr + 1; end, end
end
check('outside the region of a fixed start rejected (DK, DU, DT; c < -20 for the conditional integral)',double(nerr ~= numel(calls)),0);
F = ltratio_cdf(-3.41045363,'DT',-30,struct('g0',3/sqrt(60)));     % the paper's table of power with the initial observation, xi = 3 at c = -30
check('a fixed start at c = -30 with (c g0)^2/2 <= 100: DT, xi = 3, vs the line integral',abs(F - 0.975053418863),1e-10);
F = [ltratio_cdf(-2,'DU',-5,struct('g0',1.5)) ltratio_cdf(-3,'DT',-5,struct('g0',2))];
check('inside the region of a fixed start: DU c=-5 g0=1.5, DT c=-5 g0=2 vs the line integral',max(abs(F - [0.864977240077167 0.2889970388707237])),1e-10);
% the stationary start near the unit root (November 2026): for c > -1e-4 (and c > 0) an error, the same message in the C
% and Julia engines (DFL_CSTAT, CSTAT); at c = -1e-4 the values of the Julia engine (quickcheck.jl, section 18b); and the
% coefficient statistic of case SK at c = -0.001, which the lower end EB e^-8 of the conditional integral had put at 6e-6
nerr = 0;
calls = {@() ltratio_cdf(-1.94087,'SK',-9e-5), @() ltratio_pdf(-2,'SU',-9e-5), @() ltratio_mom(2,'ST',-9e-5), @() coef_cdf(-8.039,'SK',-9e-5), ...
         @() adfmax_cdf(-2.4,'SU',-9e-5), @() ltratio_cdf(-2,'SU',1), @() tr_limdist(-1,'SK',-9e-5)};
for k = 1:numel(calls)
    try, calls{k}(); catch e, if strcmp(e.identifier,'ltratio:args') && contains(e.message,'requires c <= -1e-4'), nerr = nerr + 1; end, end
end
check('stationary start with c > -1e-4 rejected (cdf, pdf, mom, coef, ADFmax, c > 0, TR_LIMDIST)',double(nerr ~= numel(calls)),0);
F = [ltratio_cdf(-1.94087,'SK',-1e-4) coef_cdf(-8.039,'SK',-1e-4) coef_cdf(-8.039,'SK',-1e-3)];
check('stationary start: SK c=-1e-4 F(-1.94087) and coef F(-8.039), coef at c=-0.001, vs Julia',max(abs(F - [0.028490044498374834 0.00032997244558452543 0.0010438037922007037])),1e-11);

fprintf('\n=== 7. tails with c0 ~= 0 ===\n');
% a far root of X(v) = c0 v^2 + q v + 1/2 was dropped by BREAKPTS, which ran the near-side form
% past it (an error, 'Re alpha <= 0 on a ray', or a wrong value); the line integral's absolute
% error floor (1e-14 to 1e-21) is handled by the conditional integral below 1e-3 (LTRATIO_DIST)
nerr = 0; f = zeros(1,4);
try, f(1) = ltratio_pdf(9,'SU',-1.2,struct('c0',-1.2)); catch, nerr = nerr + 1; end
try, f(2) = ltratio_pdf(-9,'SU',-1.2,struct('c0',1.2)); catch, nerr = nerr + 1; end
try, f(3) = ltratio_pdf(3,'SU',-1.2,struct('c0',-0.05)); catch, nerr = nerr + 1; end
try, f(4) = ltratio_pdf(3,'SU',-1.2,struct('c0',-0.05,'method','single')); catch, nerr = nerr + 1; end
check('no error at far roots: SU c=-1.2, c0=-1.2 (x=9), c0=1.2 (x=-9), c0=-0.05 (x=3)',nerr,0);
g = ltratio_pdf([3 -9],'SU',-1.2,struct('c0',-0.05,'method','double'));
check('SU c=-1.2 c0=-0.05, x=3: single (far root kept) vs double (relative)',abs(f(4)-g(1))/g(1),1e-10);
g = ltratio_pdf(-9,'SU',-1.2,struct('c0',1.2,'method','double','verified',true));
fv = ltratio_pdf(-9,'SU',-1.2,struct('c0',1.2,'verified',true));
check('SU c=-1.2 c0=1.2, x=-9 (f = 6.2e-16): verified (tail by the conditional integral) vs double (rel)',abs(fv-g)/g,1e-10);
check('SU c=-1.2 c0=1.2, x=-9: the default (one evaluation, the line integral) to its absolute accuracy',abs(f(2)-g),1e-13);
xs = [-9 -8 -6 3 4 5]; fd = ltratio_pdf(xs,'DT',0,struct('c0',-2,'method','double')); fn = ltratio_pdf(xs,'DT',0,struct('c0',-2));
check('DT c=0 c0=-2: default vs double in both tails, f from 4e-11 to 2.5e-6 (relative)',max(abs(fn-fd)./abs(fd)),1e-10);
fs = ltratio_pdf([8.5 9 10 12],'SU',-1.2,struct('c0',-1.2,'method','single'));    % noise below the floor
fs2 = ltratio_pdf([8.5 9 10 12],'SU',-1.2,struct('c0',-1.2,'method','double'));
check('no negative density (single and double, SU c=c0=-1.2, x = 8.5 to 12)',nnz([fs fs2] < 0),0);

fprintf('\n=== 8. far from the unit root (c = -100, -480) ===\n');
% points of julia/largec_ref.txt: the line integral of the Julia engine at x = mean -/+ 2 sd (and the mean, mean - 2 sd or
% mean + 2 sd with c0 = c), stable there under node doubling and a second contour to 1e-11 and checked against the moments
% of the density; the C engine agrees to 3e-12.  Columns: case, c, c0, x, F, f.
R = {'DK',-100,0,-8.1158889037968063,0.027930433158222967,0.11252994401966537;
         'DK',-100,0,-6.0960815949861598,0.98301181585508535,0.0973083169744004;
         'DK',-480,0,-16.51010789806336,0.025271382454033042,0.11080913614668837;
         'DK',-480,0,-14.505949841882353,0.97990751069517712,0.10384302365231887;
         'DU',-100,0,-8.1977216833376261,0.027859907170031839,0.11141083483122313;
         'DU',-100,0,-6.1584167119691813,0.98293147128060032,0.096558380445607037;
         'DU',-100,-100,-0.21208273987613141,0.50006426249294456,0.40002718640069679;
         'DU',-480,0,-16.544588887079655,0.025264367409495114,0.1105736857403487;
         'DU',-480,0,-14.536285891363056,0.97990003547381099,0.10364176672137337;
         'DT',-100,0,-8.2807466580406768,0.027788420087490939,0.11032467408272117;
         'DT',-100,0,-6.2221991829084509,0.98284995658755725,0.09583545082146698;
         'DT',-480,0,-16.579197915140668,0.025257322766800039,0.11033966693700868;
         'DT',-480,0,-14.566759193817521,0.97989252741052002,0.10344185251816358;
         'DT',-480,-480,-0.1614374354442383,0.50001094860671369,0.39916003549050727;
         'SK',-100,0,-8.1345759290815689,0.027909091546056364,0.11238539842244402;
         'SK',-100,0,-6.112450529570606,0.98298680339743949,0.097252491358002485;
         'SK',-100,-100,-2.0649606382146684,0.022741842069679885,0.054115422291525717;
         'SK',-480,0,-16.518418353305258,0.025269052056467789,0.11077838860789881;
         'SK',-480,0,-14.513747602376007,0.97990497412921096,0.1038207463566452;
         'SU',-100,0,-8.2157949655590858,0.027839255294828225,0.11127772576421112;
         'SU',-100,0,-6.1743316068020047,0.98290721447614993,0.096509162881993307;
         'SU',-480,0,-16.552845104520678,0.025262052925418277,0.11054348390585496;
         'SU',-480,0,-14.544036869849391,0.97989751590232688,0.10361993724199081;
         'ST',-100,0,-8.297560492671316,0.027768836137713082,0.11020225889150781;
         'ST',-100,0,-6.2370028092040757,0.98282689190683425,0.095791024678786446;
         'ST',-480,0,-16.587333794485499,0.025255034331478576,0.11031001188710093;
         'ST',-480,0,-14.574397190915571,0.97989003602366942,0.10342044371617426;
         'ST',-480,-480,1.8375574995319197,0.97724540216381417,0.05402325138485399;
         'DG',-480,0,-15.446799628855823,0.49999999999999922,0.77592587896902709};
M = {'DK',-100,[-7.105985249391483 50.75000271236464 -364.27855144170286 2627.9374847545978];
         'DK',-480,[-15.508028869972858 240.75000002654633 -3741.3563012709815 58202.937499840336];
         'DU',-100,[-7.1780691976534037 51.784600204191079 -375.47137671803972 2736.1126747582134];
         'DU',-480,[-15.540437389221355 241.75727430597181 -3764.8607550776778 58690.974378251703];
         'SU',-100,[-7.1950632861805452 52.029408982464922 -378.12955168762568 2761.8967524120471];
         'SU',-480,[-15.548440987185034 242.0062240397057 -3770.6744408577106 58811.781106780028]};
Q = {'DU',-100,0.05,-8.0448078491061814;
         'DU',-480,0.05,-16.379729652996637};
if quick, R = R([3 8 12 14 18 22 26 29],:); M = M([4 6],:); Q = Q(2,:); end
for k = 1:size(R,1)
    cas = R{k,1}; c = R{k,2}; c0 = R{k,3}; x = R{k,4}; o = struct('c0',c0);
    Fs = ltratio_cdf(x,cas,c,setfield(o,'method','single')); fs = ltratio_pdf(x,cas,c,setfield(o,'method','single'));
    [Fa,ia] = ltratio_cdf(x,cas,c,o); fa = ltratio_pdf(x,cas,c,o);
    check(sprintf('%s c=%g c0=%g x=%.4f: F, f single vs Julia',cas,c,c0,x),max(abs([Fs-R{k,5} fs-R{k,6}])),1e-11);
    check(sprintf('%s c=%g c0=%g x=%.4f: F, f auto vs Julia (%.1f s)',cas,c,c0,x,ia.time),max(abs([Fa-R{k,5} fa-R{k,6}])),1e-10);
end
for k = 1:size(M,1)
    m = ltratio_mom(4,M{k,1},M{k,2});
    check(sprintf('%s c=%g: ltratio_mom(4) vs Julia (relative)',M{k,1},M{k,2}),max(abs(m - M{k,3})./abs(M{k,3})),1e-10);
end
for k = 1:size(Q,1)
    q = ltratio_icdf(Q{k,3},Q{k,1},Q{k,2});
    check(sprintf('%s c=%g: ltratio_icdf(%g) vs Julia',Q{k,1},Q{k,2},Q{k,3}),abs(q - Q{k,4}),1e-8);
end

fprintf('\n=== 9. the double zero x = -sqrt(2 c0) (c = -2.5, c0 = 2) ===\n');
% the two zeros of a(b) coincide at x = -2: the conditional integral of the other engines returned NaN with a stationary
% start (a node with a(b) = 0 exactly, where K0 is infinite) and an infinite density in cases DU, DT (right in case DK
% only); the default must be finite (Inf for the density of case DK at x = -2), its distribution function within 1e-9
% of the line integral's, its density within 1e-4 (the density has a cusp there: just below x = -2 the line
% integral's density is off by up to 5e-5, the conditional integral's is stable to 1e-13; at x = -2 both are good
% to about 1e-6 only)
cs = {'DK','DU','DT','SK','SU','ST'}; xs = [-2, -2 + 1e-9, -2.000001, -1.999];
if quick, cs = {'DT','SU'}; xs = [-2, -2.000001]; end
wF = 0; wf = 0; nb = 0;
xall = xs;
for j = 1:numel(cs)
    o = struct('c0',2); xs = xall; if strcmp(cs{j},'DK'), xs = xall([1 end]); end   % (case DK: the conditional integral near the double zero is slow)
    Fa = ltratio_cdf(xs,cs{j},-2.5,o); fa = ltratio_pdf(xs,cs{j},-2.5,o);
    Fs = ltratio_cdf(xs,cs{j},-2.5,setfield(o,'method','single')); fs = ltratio_pdf(xs,cs{j},-2.5,setfield(o,'method','single'));
    dkinf = strcmp(cs{j},'DK') & xs == -2;
    nb = nb + nnz(~isfinite(Fa)) + nnz(~isfinite(fa(~dkinf))) + nnz(~isinf(fa(dkinf)));
    wF = max([wF abs(Fa - Fs)]); wf = max([wf abs(fa(~dkinf) - fs(~dkinf))./fs(~dkinf)]);
end
check('the double zero: values not finite (but the density of DK at x = -2)',nb,0);
check('the double zero: F auto vs line integral',wF,1e-9);
check('the double zero: f auto vs line integral (relative)',wf,1e-4);

fprintf('\n=== 10. the applications of the paper ===\n');
% against the same functions of julia/DFExact.jl (adfmax_cdf, coef_cdf, ltratio_region, ltratio_lincdf, power_envelope): the
% ADFmax statistic (distribution function, density, a quantile), the coefficient statistic, the region P{U <= g(V)} for a
% line, the boundary of the weighted-average-power test and the parabola of the t-ratio, the Gil-Pelaez inversion of
% U - aV, the power envelope (case DK and the GLS trend case)
A = {'DU',0,-2.5,0.042709500744894657,0.10563355056449016; 'DT',-10,-3,0.34351886451279978,0.59340743849891675;
     'SU',-10,-2.5,0.45543929180199272,0.69145494355731718; 'ST',-25,-4,0.36647812431909832,0.66374680023109822};
if quick, A = A([1 3],:); end
for k = 1:size(A,1)
    F = adfmax_cdf(A{k,3},A{k,1},A{k,2}); f = adfmax_pdf(A{k,3},A{k,1},A{k,2});
    check(sprintf('ADFmax %s c=%g x=%g: F and f vs Julia',A{k,1},A{k,2},A{k,3}),max(abs([F f] - [A{k,4} A{k,5}])),1e-10);
end
check('ADFmax DU c=0: the 5% point vs Julia',abs(adfmax_icdf(0.05,'DU',0) + 2.4354164994325092),1e-8);
C = {'DU',-5,-14,0.1938495145682877; 'SU',-10,-14,0.45401150754561903; 'DG',-13.5,-16.59,0.51132404681639121};
for k = 1:size(C,1)
    check(sprintf('coefficient statistic %s c=%g: coef_cdf(%g) vs Julia',C{k,1},C{k,2},C{k,3}),abs(coef_cdf(C{k,3},C{k,1},C{k,2}) - C{k,4}),1e-10);
end
check('coefficient statistic DU c=0: the 5% point vs Julia',abs(coef_icdf(0.05,'DU',0) + 14.093590831267115),1e-7);
check('region U + 3.5V <= -0.2639, c=0: ltratio_region vs Julia',abs(ltratio_region(struct('type','line','a',-3.5,'b',-0.2639),'DK',0) - 0.049996781931320872),1e-10);
check('region U + 3.5V <= -0.2639, c=0: ltratio_lincdf (Gil-Pelaez) vs Julia',abs(ltratio_lincdf(-0.2639,-3.5,'DK',0) - 0.049996781931325174),1e-10);
check('region of the WAP test, k=3.6658, omega=10, c=-10, vs Julia',abs(ltratio_region(struct('type','wap','k',3.6658,'omega',10),'DK',-10) - 0.75806440265373276),1e-10);
check('region of the parabola x = -1.94087, c=-7: the t-ratio (line integral)', ...
      abs(ltratio_region(struct('type','parabola','x',-1.94087),'DK',-7) - ltratio_cdf(-1.94087,'DK',-7,S)),1e-11);
check('power envelope DK c=-7 vs Julia',abs(power_envelope(-7,'DK') - 0.49940198260314117),1e-10);
if ~quick, check('power envelope DG c=-13.5 vs Julia',abs(power_envelope(-13.5,'DG') - 0.51441347319065545),1e-10); end
% the linear form where the inversion of Gil-Pelaez is not trusted (small |a|, c > 0): the conditional integrals (the
% region in case DK, the conditional integral given R = U + 1/2 in case DG); at a = 0 the closed form P(U <= b) =
% erf(sqrt(kappa (b + 1/2))), kappa = 2c/(e^(2c) - 1), 0 below the support
w = 0;
for c = [-5 0 2]
    if c == 0, kap = 1; else, kap = 2*c/expm1(2*c); end
    b = -0.5 + [-0.1 1e-6 1e-3 0.05 0.5 2]/kap;
    w = max(w,max(abs(ltratio_lincdf(b,0,'DK',c) - erf(sqrt(max(kap*(b + 0.5),0))))));
end
check('linear form a = 0, case DK, c in {-5, 0, 2}: the closed form',w,1e-13);
check('linear form a = 0, case DK: 0 below the support',ltratio_lincdf(-0.6,0,'DK',-5),0);
check('linear form DG a = 0.3, c = -5, b = -0.51 vs Julia (the conditional integral given R)',abs(ltratio_lincdf(-0.51,0.3,'DG',-5) - 0.291315235195),1e-11);
P = ltratio_lincdf(-0.7:0.1:1,1e-3,'DG',2);
check('linear form DG a = 1e-3, c = 2: in [0, 1] and nondecreasing',double(any(P < 0 | P > 1) || any(diff(P) < 0)),0);
if ~quick, check('power envelope DK c = -0.01 (the conditional integral) vs Julia',abs(power_envelope(-0.01,'DK') - 0.050249919491),1e-10); end
e = 0; try, ltratio_lincdf(0,0,'DK',-600); catch err, e = contains(err.message,'at least -550'); end
check('linear form: c = -600 rejected',double(~e),0);

fprintf('\n=== 11. the contour integral (a fixed start in cases DU, DT) ===\n');
% the reference values of c/test/contour_ref.txt (DFExact.tr_contour of the Julia engine; the same file checks the C
% engine): the grid of the experiment (c in {-10, -5, -2, 1}, c0 in {0, c}), two roots of X(v) with c0 > 0, c0 < 0 far
% from c, near the double root x = -sqrt(2 c0), x = 0 and x -> 0-, c = 0 and c = 2, c = 1.9 in the upper tail, far from
% the unit root (c = -94.4, -70.5, -52.5), a large noncentrality ((c g0)^2/2 = 364 at c = -18, 7200 at c = -40 in case
% DT).  Columns: model, c, c0, g0, x, F, f.  (At c = -18 the roots of D(., t) carry a relative error of about 1e-12 in
% both engines, which the noncentrality turns into a difference of 1e-13 in F = 3.5e-9.)
R0 = [
     2 -5 0 0.31622776601683794 -2.8599999999999999 0.12484216310046495 0.28821207043794761;
     2 -5 0 0.31622776601683794 1 0.99999451778890025 2.9586738529586943e-05;
     2 -5 0 1.5 -2 0.8649772400771647 0.33128663780386036;
     2 -5 0 1 0 0.99960183161946026 0.0017314891359810492;
     2 -2 -2 1.5 0 0.8039072720576913 0.26886994985279628;
     2 -10 0 0.67082039324993692 -4 0.037477901860310615 0.11793798741289027;
     2 -10 -10 0.67082039324993692 -1.5 0.16708692728951868 0.25713086867745499;
     2 1 1 2 -2.8599999999999999 0.015980404209157282 0.04203818841266057;
     2 1 1 0.5 -1.5 0.53364170675267153 0.42371259419507401;
     2 1.8999999999999999 0 1.298 3.3039999999999998 0.42280476995264499 0.22183434368147734;
     2 -94.409999999999997 0 0.050999999999999997 -7.9290000000000003 0.038360478680122989 0.1457188199279075;
     2 -17.809999999999999 2.1899999999999999 0.249 -4.1440000000000001 0.16641711211815027 0.44049716999000105;
     2 -15.66 -2.02 0.55700000000000005 -1.9430000000000001 0.97095572323765456 0.12916868891910824;
     2 0 0 1 -1.5 0.53400079951916057 0.50879318678812302;
     2 2 0 1 0.5 0.11650278571383295 0.067496786780139109;
     2 -18 0 1.5 -9.0806812422194216 3.4899822502374152e-09 2.7716662667539702e-08;
     3 -5 0 2 -3 0.28899703887072636 0.46367336787827451;
     3 -10 -10 0.67082039324993692 -1.5 0.29337189777053652 0.3590580247130955;
     3 -70.450000000000003 -4.25 0.26800000000000002 -4.2999999999999998 0.99964458398177858 0.00297854781917506;
     3 -7.5300000000000002 1.8400000000000001 0.90000000000000002 -1.9199999999999999 0.98462193489798155 0.05623151557997863;
     3 1.26 1.95 -1.5820000000000001 -1.901 0.70593114960538439 0.39461731598347738;
     3 1 0 2 -2.8599999999999999 0.13961097469619091 0.27750350843462407;
     3 -2 -2 1.5 1 0.99819132794725107 0.0057613984523278514;
     3 -52.539999999999999 1.5600000000000001 0.17199999999999999 -3.734 0.99999011635045065 0.00011222219700637244;
     3 2 0 1 -2 0.11551965155655686 0.15070023909262359;
     3 0 0 3 -2.5 0.328229118126592 0.49645984886583056;
     3 -2 0 1.5 -4 0.009514763936226283 0.031242873748537906;
     3 -20 0 1.5 -4.5 0.89704065002340361 0.33046960953081839;
     2 -5 0 1 -1e-08 0.99960183160214677 0.0017314892071813505;
     3 -40 0 3 -13 0.49617272547857838 0.72928545720918503;
];
cs3 = {'','DU','DT'};
R = R0; if quick, R = R0([1 3 8 11 17 20 24 30],:); end
for k = 1:size(R,1)
    cas = cs3{R(k,1)}; o = struct('c0',R(k,3),'g0',R(k,4));
    [F,iF] = tr_contour(R(k,5),cas,R(k,2),o,'cdf'); [f,jf] = tr_contour(R(k,5),cas,R(k,2),o,'pdf');
    check(sprintf('%s c=%g c0=%g g0=%.4g x=%.4g: F, f contour vs Julia (%.1f/%.1f s)',cas,R(k,2),R(k,3),R(k,4),R(k,5),iF.time,jf.time), ...
          max(abs([F - R(k,6), f - R(k,7)])) + (iF.status ~= 0) + (jf.status ~= 0),1e-11);
end
% the automatic choice of the classical programs is the contour integral (values of at least 1e-3; at c = -52.54, outside
% the region where the conditional integral is checked here, also the density 1.1e-4), and opts.method = 'contour' through
% LTRATIO_CDF (engine 'classical': with the folder of the MEX file ahead on the path a case string goes to the C router
% by default, whose info has no count of the contour integral)
ia = [1 3 5 11 17 19 24 26]; if quick, ia = [1 11 24]; end
R = R0;
w = 0; nc = 0;
for k = ia
    cas = cs3{R(k,1)}; o = struct('c0',R(k,3),'g0',R(k,4),'engine','classical');
    [Fa,ja] = ltratio_cdf(R(k,5),cas,R(k,2),o); fa = ltratio_pdf(R(k,5),cas,R(k,2),o);
    Fc = ltratio_cdf(R(k,5),cas,R(k,2),setfield(o,'method','contour'));
    oc = rmfield(o,'engine');
    Ft = tr_contour(R(k,5),cas,R(k,2),oc,'cdf'); ft = tr_contour(R(k,5),cas,R(k,2),oc,'pdf');
    w = max([w abs(Fa - Ft) abs(fa - ft) abs(Fc - Ft)]); nc = nc + ja.ncontour;
end
check(sprintf('auto and ''contour'' = TR_CONTOUR with a fixed start in DU, DT (%d points, %d by the contour integral)',numel(ia),nc),w + (nc ~= numel(ia)),0);
% a density below 1e-3 at c = -5 (inside the region): recomputed by the conditional integral, within its relative accuracy
fa = ltratio_pdf(1,'DU',-5,struct('g0',1/sqrt(10))); fd = tr_limdist(1,'DU',-5,struct('g0',1/sqrt(10)),'pdf');
check('DU c=-5 g0=0.316 x=1: f = 3.0e-5 by auto = conditional integral, vs contour (relative)',abs(fa - fd) + abs(fa/2.9586738529586943e-05 - 1),1e-9);
% the scope: an error elsewhere (the same message in the C and Julia engines), in TR_CONTOUR too
nerr = 0;
calls = {@() ltratio_cdf(-2,'DK',-5,struct('g0',1,'method','contour')), @() ltratio_pdf(-2,'DU',-5,struct('method','contour')), ...
         @() ltratio_cdf(-2,'SU',-5,struct('method','contour')), @() ltratio_cdf(-2,'DG',0,struct('method','contour'))};
for k = 1:numel(calls)
    try, calls{k}(); catch e, if strcmp(e.identifier,'ltratio:args') && contains(e.message,'(the contour integral) applies to a fixed start g0 != 0 in cases DU, DT'), nerr = nerr + 1; end, end
end
try, tr_contour(-2,'DK',-5,struct('g0',1)); catch e, if strcmp(e.identifier,'tr_contour:args'), nerr = nerr + 1; end, end
try, tr_contour(-2,'DU',-5,struct('g0',0)); catch e, if strcmp(e.identifier,'tr_contour:args'), nerr = nerr + 1; end, end
check('method ''contour'' outside a fixed start in DU, DT rejected (DK, zero start, SU, DG; TR_CONTOUR)',double(nerr ~= numel(calls) + 2),0);
F = tr_contour([-Inf Inf NaN],'DT',-5,struct('g0',1)); f = tr_contour([-Inf Inf NaN],'DT',-5,struct('g0',1),'pdf');
check('contour: F(-Inf), F(Inf), F(NaN) = 0, 1, NaN; f = 0, 0, NaN',abs(F(1)) + abs(F(2)-1) + ~isnan(F(3)) + abs(f(1)) + abs(f(2)) + ~isnan(f(3)),0);

fprintf('\n=== 12. the default (one evaluation), verified and the warnings of the brittle regimes ===\n');
% verified reproduces the checked computation (values of the engine before November 2026, quickcheck.jl section 21)
VR = {'DU',-5,2,0,-3,0.1976948237625359; 'DK',-10,-5,0,-2,0.16726497409337868; 'SU',-5,-5,0,-1,0.44272919976045844; 'DU',-5,0,0,-4.5,0.0005234905959262681};
wv = 0; w1 = 0;
for k = 1:size(VR,1)
    o = struct('c0',VR{k,3},'g0',VR{k,4});
    wv = max(wv, abs(ltratio_cdf(VR{k,5},VR{k,1},VR{k,2},setfield(o,'verified',true)) - VR{k,6}));
    w1 = max(w1, abs(ltratio_cdf(VR{k,5},VR{k,1},VR{k,2},o) - VR{k,6}));
end
check('verified reproduces the checked computation at four points',wv,1e-11);
check('the default (one evaluation) at the same points',w1,1e-9);
% the warnings: once where they should fire, never elsewhere
nw = @(f) nwarn_(f);
% (the warnings of the classical programs: engine = 'classical' where the router of the two engines would take the general one)
fire = [nw(@() ltratio_cdf(-3,'DU',-5,struct('c0',2,'method','double'))), nw(@() ltratio_cdf([-3 1],'DK',2,struct('c0',-5,'method','double'))), ...
        nw(@() ltratio_cdf(-5,'DU',-20,struct('c0',2,'g0',0.5))), nw(@() ltratio_cdf(-2,'SK',-0.1,struct('c0',1,'method','double'))), ...
        nw(@() ltratio_icdf([0.05 0.5],'DU',-20,struct('c0',2,'g0',0.5))), nw(@() ltratio_cdf(-3.5,'DK',-25,struct('g0',0.3,'engine','classical'))), ...
        nw(@() tryf_(@() ltratio_cdf(-8.4,'DU',2,struct('c0',5,'g0',1.5,'method','contour'))))];
% (the failing contour integral at c = 2, c0 = 5, g0 = 1.5 gives 1.2e44 at x = -8.4: since November 23, 2026 the final
% check, LTRATIO_RANGECHECK, turns it into the error ltratio:range after the warning)
try, ltratio_cdf(-8.4,'DU',2,struct('c0',5,'g0',1.5,'method','contour')); id = ''; catch ME, id = ME.identifier; end
check('DU c=2 c0=5 g0=1.5, the contour integral at x = -8.4 (1.2e44): the error ltratio:range of the final check',double(~strcmp(id,'ltratio:range')),0);
quiet = [nw(@() ltratio_cdf(-3,'DU',-5,struct('c0',2))), nw(@() ltratio_cdf(-2,'DU',-5)), nw(@() ltratio_cdf(-2,'DU',-5,struct('g0',1.5))), ...
         nw(@() ltratio_cdf(-3,'DU',-5,struct('c0',2,'method','double','verified',true))), nw(@() ltratio_cdf(-6,'DK',-10,struct('g0',1)))];
check(sprintf('the warnings: once each where they should (%s), none elsewhere (%s)',num2str(fire),num2str(quiet)),sum(fire ~= 1) + sum(quiet),0);
% case DU with a fixed start, c0 > 3 at c > 1.5: the contour integral fails (1e44), the default takes the conditional integral
Fa = ltratio_cdf(-8.4075,'DU',2,struct('c0',5,'g0',1.5)); Fv = ltratio_cdf(-8.4075,'DU',2,struct('c0',5,'g0',1.5,'method','double','verified',true));
check('DU c=2 c0=5 g0=1.5: the default (the conditional integral, not the failing contour) vs verified',abs(Fa - Fv),1e-9);

fprintf('\n=== 13. the upper tail above the unit root ===\n');
% (the references: the independent computation by conditioning on J(1); the classical programs before: 0.8415, -1.3674,
% 1.1422, 0.6076, 1.0479)
R = [10 2 1e4 0.743623945123; 10 2 3e4 0.999338346257; 7 2 2e3 0.999999671714; 10 9 1e3 0.636120922790; 10 9 3e3 0.993550212928];
for k = 1:size(R,1)
    o = struct('c0',R(k,2),'engine','classical');
    F = ltratio_cdf(R(k,3),'DK',R(k,1),o);
    check(sprintf('DK c=%g c0=%g x=%g: F against the independent computation',R(k,1),R(k,2),R(k,3)),abs(F - R(k,4)),1e-10);
end
F = ltratio_cdf([1e3 3e3 1e4],'DK',10,struct('c0',9,'engine','general'));
check('DK c=10 c0=9 x=1e3, 3e3, 1e4: the general engine (the line cut to the rate with c0)',max(abs(F - [0.636120922790 0.993550212928 1])),1e-10);
try
    ltratio_cdf(1e7,'DK',16,struct('c0',2,'engine','classical')); ok = false;
catch ME
    ok = ~isempty(strfind(ME.message,'computes for c <= 14 only')) && ~isempty(strfind(ME.message,'c = 16, c0 = 2, x = 1e+07'));
end
check('DK c=16 c0=2 x=1e7: refused, with the arguments and the region',double(~ok),0);
try
    ltratio_cdf(1e7,'DK',16,struct('c0',-2,'verified',true,'engine','classical')); ok = false;
catch ME
    ok = ~isempty(strfind(ME.message,'computes for c <= 14 only'));
end
check('DK c=16 c0=-2 x=1e7, verified (the rule''s conditional integral): refused',double(~ok),0);
% cases DU, DT far above the unit root: the conditional integral at x > 300, against the large-c law; above c = 14 refused
F = [ltratio_cdf(4.07e5,'DU',13,struct('engine','classical')) ltratio_cdf(4.766e5,'DT',13,struct('engine','classical'))];
check('DU, DT c=13 far above the unit root: against the large-c law',max(abs(F - [0.9545202704 0.9973023556])),1e-6);
try
    ltratio_cdf(4.4e6,'DT',16,struct('engine','classical')); ok = false;
catch ME
    ok = ~isempty(strfind(ME.message,'computes for c <= 14 only'));
end
check('DT c=16 x=4.4e6 (the line integral exhausted the memory): refused',double(~ok),0);
try
    ltratio_rangecheck('ltratio_cdf','case DK, c = 1','x',[0.5 1],[0.3 1 + 2e-8]); ok = false;
catch ME
    ok = strcmp(ME.identifier,'ltratio:range') && ~isempty(strfind(ME.message,'ltratio_cdf: case DK, c = 1, x = 1: the computed distribution function is 1.00000002'));
end
ok = ok && isequal(ltratio_rangecheck('ltratio_cdf','case DK, c = 1','x',[0.5 1],[0.3 1 + 5e-9]),[0.3 1 + 5e-9]);
check('LTRATIO_RANGECHECK: an error outside [0, 1] by more than 1e-8, the values otherwise',double(~ok),0);

fprintf('\n%d of %d checks failed (%.0f s).\n',nfail,ntest,toc(t00));
if strict && nfail > 0, error('ltratio_test:fail','%d check(s) failed.',nfail); end

    function check(name,err,tol)
        ntest = ntest + 1; pass = isfinite(err) && err <= tol;
        if ~pass, nfail = nfail + 1; end
        mk = {'FAIL','pass'};
        fprintf('  %-72s %9.2e  (tol %7.1e)  %s\n',name,err,tol,mk{pass+1});
    end
end

% ---------------------------------------------------------------------------
function tryf_(f)
% f(), its error ignored (the warnings it gives before the error are counted by NWARN_)
try, f(); catch, end
end
function n = nwarn_(f)
% the number of warnings of a brittle regime (LTRATIO_BRITTLE) that f() gives
ws = warning('on','ltratio:brittle');
s = evalc('f();');
warning(ws);
n = numel(regexp(s,'checks\s+each\s+value'));   % (evalc wraps the text)
end
function [x,w] = ghnodes(m)
% Gauss-Hermite rule for the standard normal: E f(Z) ~ sum w f(x)
b = sqrt((1:m-1)); [V,D] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(D)); w = V(1,ix).^2; w = w(:);
end

function [x,w] = glnodes(n,a,b)
% Gauss-Legendre rule on [a, b]
k = (1:n-1)'; be = k./sqrt(4*k.^2 - 1);
[V,E] = eig(diag(be,1) + diag(be,-1));
[x,ix] = sort(diag(E)); w = 2*V(1,ix)'.^2;
x = a + (b - a)*(x + 1)/2; w = w*(b - a)/2;
end
