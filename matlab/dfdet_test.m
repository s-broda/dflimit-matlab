function nfail = dfdet_test(quick,strict)
%DFDET_TEST  Validation of the general engine against the Julia engine DFLimit.DFDet.
%
%   nfail = DFDET_TEST compares the general engine of these programs (DFDET_*,
%   PT_*, LM_*, MZA_*, MSB_*, and the router of LTRATIO_CDF and the other
%   user-level routines with opts.engine) with the reference values of the
%   Julia engine in validation/det_ref.txt (written by validation/det_ref.jl),
%   prints each comparison with its tolerance and a PASS/FAIL mark and
%   returns the number of failures; DFDET_TEST(true) skips the slow groups
%   (the point-optimal statistic, whose grid has 60000 points, the envelope,
%   the moments with c0 ~= 0, c = 2), DFDET_TEST(quick,true) raises an error
%   at the end when a check has failed.  The groups:
%    1. the designs (every constructor: the polynomial trends of degrees 0 to
%       4, Perron's models with one, two and three breaks, a function design,
%       Legendre trends, designs from a matrix with and without given breaks,
%       the statistic that detrends first): the sizes of the functionals, the
%       number of roots and their signs, the eigenvalues of A; GLS detrending;
%    2. the determinant of the joint moment generating function, its roots,
%       t_* (zero, fixed and stationary start, GLS detrending);
%    3. the distribution function and the density: c in {0, -5, -13.5, -20,
%       -50, 2, -300}, the three starts, GLS detrending, local null values of
%       both signs, far tails, the infinite density at the tangency;
%    4. quantiles at the levels 0.01, 0.05, 0.10, 0.5, 0.95;
%    5. moments;
%    6. the coefficient statistic, the linear form (the law of U at the edge
%       of its support), the law of V, ADFmax, the power envelope, local
%       power, regions;
%    7. the companions: the point-optimal statistic and the envelope of any
%       design, the LM statistic, MZ_alpha, MSB;
%    8. the router: LTRATIO_CDF, LTRATIO_PDF, LTRATIO_ICDF, LTRATIO_MOM,
%       COEF_CDF, LTRATIO_LINCDF, ADFMAX_CDF, POWER_ENVELOPE, LTRATIO_REGION with
%       opts.engine 'auto', 'classical' and 'general', with case strings and
%       with designs;
%    9. the path tables against the computed paths, non-finite x, argument
%       errors.
%   Tolerances: 1e-8 absolute for distribution functions, densities and
%   probabilities, 1e-6 for quantiles and critical values, 1e-8 relative for
%   moments, determinants, roots and eigenvalues, 1e-10 relative for t_*.
%   About ten minutes (two in the quick version).
%
%   See also LTRATIO_TEST, DFDET_CDF.
if nargin < 1 || isempty(quick), quick = false; end
if nargin < 2 || isempty(strict), strict = false; end
here = fileparts(mfilename('fullpath'));
nfail = 0; ntest = 0; t00 = tic;
G = containers.Map();                                  % per group: the number of comparisons and the largest deviation
ws = warning('off','dfdet:line'); wr = warning('off','dfdet:rounding'); wc = warning('off','dfdet:contour'); wb = warning('off','ltratio:brittle');
cl = onCleanup(@() warning([ws wr wc wb]));
L = strsplit(strtrim(fileread(fullfile(here,'validation','det_ref.txt'))),newline);
slow = {'ptq','ptcdf','ptpdf','ptenv','ptpow','env','renv'};
cur = '';
for i = 1:numel(L)
    l = strtrim(L{i}); if isempty(l), continue, end
    p = strtrim(strsplit(l,'|')); if numel(p) < 3, p{3} = ''; end
    h = strsplit(p{1},' '); func = h{1}; dtok = h{2}; c = str2double(h{3}); o = opts_(h{4});
    args = num_(p{2}); ref = num_(p{3});
    g = group_(func);
    if ~strcmp(g,cur), cur = g; fprintf('\n=== %s ===\n',g); end
    if quick && (any(strcmp(func,slow)) || (strcmp(func,'mom') && isfield(o,'c0')) || c == 2 || c == -300), continue, end
    t0 = tic;
    [v,tol,rel] = compute_(func,dtok,c,o,args);
    tt = toc(t0);
    v = v(:); ref = ref(:);
    if strcmp(func,'err')
        err = v;                                       % (0 where the call failed with the message of the line)
    elseif strcmp(func,'roots') && 2*numel(v) == numel(ref)
        err = 0;                                       % (each reference root against the nearest computed one)
        rr = complex(ref(1:2:end),ref(2:2:end));
        for k = 1:numel(rr), err = max(err,min(abs(v - rr(k)))/max(1,abs(rr(k)))); end
    elseif numel(v) ~= numel(ref), err = Inf;
    elseif strcmp(func,'dval')
        err = abs(complex(v(1),v(2)) - complex(ref(1),ref(2)))/abs(complex(ref(1),ref(2)));
    elseif strcmp(func,'roots')
        err = 0;                                       % (each reference root against the nearest computed one)
        rr = complex(ref(1:2:end),ref(2:2:end));
        for k = 1:numel(rr), err = max(err,min(abs(v - rr(k)))/max(1,abs(rr(k)))); end
    else
        inf_ = isinf(ref);
        err = max([0; abs(v(~inf_) - ref(~inf_))./(rel*max(1,abs(ref(~inf_))) + ~rel)]);
        if any(v(inf_) ~= ref(inf_)), err = Inf; end
    end
    check(sprintf('%-8s %-16s c = %-6g %-30s (%d values, %.2f s)',func,dtok,c,h{4},numel(ref),tt),err,tol);
    if ~isKey(G,g), G(g) = [0 0]; end
    gv = G(g); G(g) = [gv(1) + numel(ref), max(gv(2),err)];
end

fprintf('\n=== 9. the path tables, non-finite x, argument errors ===\n');
f = fullfile(here,'dfdet_paths.mat');
if exist(f,'file')
    T = load(f); nm = {'DK','DU','DT'};
    for q = -1:1
        Mo = dfdet_startmodel(dfdet_trend(q),0,'zero',0,[]);
        CP = dfdet_cpath(Mo,0,struct('phi',pi/4));      % (computed: an explicit angle bypasses the table)
        Tq = T.(nm{q + 2});
        d1 = max(abs(Tq.t(:) - CP.t(:)))/max(abs(CP.t)); d2 = max(abs(Tq.R(:) - CP.R(:)))/max(abs(CP.R(:)));
        check(sprintf('path table %s against the computed path (nodes, roots; relative)',nm{q + 2}),max(d1,d2),1e-10);
        xs = [-4 -2 0 1];
        Ft = dfdet_cdf(xs,dfdet_trend(q)); Fc = dfdet_cdf(xs,dfdet_trend(q),0,struct('phi',pi/4));
        check(sprintf('path table %s: F on the table against F on the computed path',nm{q + 2}),max(abs(Ft - Fc)),1e-13);
    end
else
    check('path tables dfdet_paths.mat present',Inf,0);
end
F = dfdet_cdf([-Inf NaN -2 Inf],dfdet_trend(1),-5); f = dfdet_pdf([-Inf NaN -2 Inf],dfdet_trend(1),-5);
check('F(-Inf), F(NaN), F(Inf) = 0, NaN, 1; f = 0, NaN, 0',abs(F(1)) + ~isnan(F(2)) + abs(F(4) - 1) + abs(f(1)) + ~isnan(f(2)) + abs(f(4)),0);
E = {@() dfdet_cdf(-2,dfdet_trend(1),-5,struct('gls',-13.5,'start','fixed','g0',1)), ...
     @() dfdet_cdf(-2,dfdet_trend(1),-5,struct('gls',-13.5,'start','stationary')), ...
     @() dfdet_cdf(-2,dfdet_detrend_first(dfdet_perron('B',0.5)),-5,struct('start','stationary')), ...
     @() dfdet_cdf(-2,dfdet_detrend_first(dfdet_perron('B',0.5)),-5,struct('gls',-13.5)), ...
     @() dfdet_mom(2,dfdet_detrend_first(dfdet_perron('B',0.5)),-5), ...
     @() dfdet_cdf(-2,dfdet_trend(1),0.5,struct('start','stationary')), ...
     @() dfdet_cdf(-2,dfdet_trend(1),-5,struct('g0',1,'start','stationary')), ...
     @() dfdet_region(struct('type','wap','k',3.6658,'omega',10),dfdet_trend(-1),0), ...
     @() dfdet_power_envelope(-7,dfdet_perron('A',0.5)), @() lm_quantile(0.05,dfdet_trend(0)), ...
     @() ltratio_cdf(-2,dfdet_perron('C',0.3),0,struct('engine','classical')), ...
     @() ltratio_cdf(-2,'DU',0,struct('method','single','engine','general')), ...
     @() ltratio_cdf(-2,'DU',0,struct('engine','fast')), ...
     @() dfdet_design_from_matrix([ones(500,1), (1:500)'/500, double((1:500)' == 250)])};
nerr = 0;
for k = 1:numel(E)
    try, E{k}(); catch ME, if strncmp(ME.identifier,'dfdet:',6) || strcmp(ME.identifier,'ltratio:args'), nerr = nerr + 1; end, end
end
check(sprintf('argument errors raised (%d of %d)',nerr,numel(E)),numel(E) - nerr,0);
% the refusals of November 2026: each call fails with an error whose message contains the text (no silently wrong value,
% no growing bracket; the refusals of the point-optimal statistic's fixed grid with opts.method = 'grid' (GR): the default
% up to c = 4, the single inversion of Bao and Kan, has none of them)
GR = struct('method','grid');
R = {{@() lm_cdf(-0.1,dfdet_trend(1),11), 'must be at most 10'}, {@() lm_quantile(0.05,dfdet_trend(1),11), 'must be at most 10'}, ...
     {@() lm_power([-5 11],dfdet_trend(1)), 'must be at most 10'}, ...
     {@() pt_cdf(10,dfdet_trend(1),-13.5,3.7,GR), 'does not resolve'}, {@() pt_cdf(10,dfdet_trend(0),-7,3.1,GR), 'does not resolve'}, ...
     {@() pt_cdf(10,dfdet_perron('B',0.5),-17.8,4.4), 'does not resolve'}, {@() pt_cdf(10,dfdet_perron('C',0.5),-20,4.3), 'does not resolve'}, ...
     {@() pt_cdf(1e5,dfdet_trend(1),-13.5,0,GR), 'beyond the range'}, {@() pt_cdf(10,dfdet_trend(1),-13.5,0,struct('method','fast')), 'opts.method must be'}, ...
     {@() pt_quantile(1.5,dfdet_trend(1),-13.5), 'must lie in (0, 1)'}, {@() dfdet_lin_icdf(1.5,-3.5,dfdet_trend(-1),0), 'must lie in (0'}, ...
     {@() dfdet_adfmax_icdf(1.5,dfdet_trend(0),0), 'must lie in (0'}, {@() dfdet_icdf(NaN,dfdet_trend(1),0), 'must lie in (0'}, ...
     {@() lm_quantile(0,dfdet_trend(1)), 'must lie in (0, 1)'}, {@() msb_quantile(1.5,dfdet_trend(1),0), 'must lie in (0, 1)'}, ...
     {@() ltratio_icdf(1.5,'DU',-5), 'must lie in (0'}, {@() coef_icdf(1.5,'DU',-5), 'must lie in (0'}, ...
     {@() pt_power(-5,dfdet_trend(1),-10,struct('alpha',1.5)), 'must lie in (0, 1)'}, {@() power_envelope(-5,'DK',struct('alpha',1.5)), 'must lie in (0, 1)'}, ...
     {@() dfdet_adfmax_cdf(-2,dfdet_trend(0),-5,struct('gls',-13.5)), 'least-squares detrending only'}, ...
     {@() dfdet_adfmax_cdf(-2,dfdet_trend(0),-5,struct('c0',1)), 'null value phi = 1'}, ...
     {@() adfmax_cdf(-2,dfdet_trend(1),-5,struct('gls',-13.5)), 'least-squares detrending only'}, ...
     {@() adfmax_cdf(-2,'DU',-5,struct('g0',1.5,'engine','classical')), 'with a zero or a stationary start'}, ...
     {@() adfmax_cdf(-2,'SU',-5,struct('g0',1.5)), 'exclude each other'}, {@() adfmax_cdf(-2,'DU',-5,struct('c0',1)), 'null value phi = 1'}, ...
     {@() ltratio_cdf(-2,'DU',-5,struct('start','stationary')), 'apply to a design'}, {@() ltratio_cdf(-2,'DU',-5,struct('gls',-7)), 'apply to a design'}, ...
     {@() ltratio_cdf(-2,dfdet_trend(0),-5,struct('start','zero','g0',1)), 'exclude each other'}, ...
     {@() ltratio_lincdf(-0.26,-3.5,'DK',-5,struct('g0',1.5)), 'a fixed start takes a design'}, ...
     {@() ltratio_lincdf(-0.26,-3.5,'DK',-5,struct('c0',1)), 'opts.c0 does not apply'}, ...
     {@() ltratio_lincdf(-0.26,-3.5,dfdet_trend(-1),-5,struct('g0',1.5,'engine','classical')), 'cases DK and DG only'}, ...
     {@() ltratio_region(struct('type','line','a',-3.5,'b',-0.26),'DK',-5,struct('g0',1.5)), 'a fixed start takes a design'}, ...
     {@() ltratio_region(struct('type','line','a',-3.5,'b',-0.26),dfdet_trend(-1),-5,struct('g0',1.5,'engine','classical')), 'with a zero start'}, ...
     {@() power_envelope(-5,dfdet_trend(-1),struct('g0',1)), 'zero start'}, {@() power_envelope(-5,'DK',struct('g0',1)), 'a fixed start takes a design'}, ...
     {@() dfdet_cdf(-2,dfdet_trend(0),-5,struct('hfac',1e-4)), 'the outer contour would need'}, {@() dfdet_mom(2,dfdet_trend(0),1e5), 'up to c = 16'}};
for k = 1:numel(R)
    msg = 'no error';
    try, R{k}{1}(); catch ME, msg = ME.message; end
    check(sprintf('refused: %s (%.60s)',func2str(R{k}{1}),msg),double(~contains(msg,R{k}{2})),0);
end
% the guard on the nodes of the U-shapes before the outer sum (USHAPE_GUARD of DFDET_CONTOUR): the design (1 + r,
% 1{r < 0.4}), whose positive root runs off on the vertical line, and d(r) = 1 + r with a fixed start are refused before
% any inner integral (the guard in the U-shapes refused after 11 s; the time is now that of the outer path, about 1 s)
Dl = dfdet_determ([0 0.4 1],{[1 1; 1 0],[1 0; 1 0]},'1 + r, 1{r < 0.4}'); Dr = dfdet_determ([0 1],{[1; 1]},'1 + r');
G2 = {{@() dfdet_cdf(-2,Dl,2), 'U-shapes of the zero start'}, {@() dfdet_cdf(-2,Dr,0,struct('start','fixed','g0',1)), 'U-shapes of the fixed start'}};
% (the time of the refusal against that of the outer path of the same design built in the same session, 4 times it and
% 1 s at most: the guard costs about one path and the refusal in the U-shapes eleven; a fixed bound of seconds failed on a
% loaded machine)
GP = {{Dl,2,'zero',0}, {Dr,0,'fixed',1}};
for k = 1:numel(G2)
    dfdet_cache('clear');
    t0 = tic; dfdet_cpath(dfdet_startmodel(GP{k}{1},GP{k}{2},GP{k}{3},GP{k}{4},[]),GP{k}{2},struct); tp = toc(t0);
    dfdet_cache('clear');
    msg = 'no error'; t0 = tic; try, G2{k}{1}(); catch ME, msg = ME.message; end; tk = toc(t0);
    check(sprintf('the guard before the outer sum: %s refused within 4 times the outer path and 1 s (%.2f s, the path %.2f s)',func2str(G2{k}{1}),tk,tp),double(~contains(msg,G2{k}{2})) + double(tk > 4*tp + 1),0);
end
% the public functions at non-finite arguments (November 2026): a limit where there is one, else an error that names the
% argument (as the C engine and DFLimit.jl); none may hang (the searches from x0 = +-Inf, the point-optimal statistic at
% c = +-Inf did not return in C and Julia)
DU = dfdet_trend(0); DT = dfdet_trend(1);
NF = {{@() dfdet_icdf(0.05,DU,-5,struct('x0',Inf)), 'start value of the search'}, ...
      {@() dfdet_icdf(0.05,DU,-5,struct('tol',Inf)), 'opts.tol must be a positive finite'}, ...
      {@() dfdet_coef_icdf(0.05,DU,-5,struct('y0',-Inf)), 'start value of the search'}, ...
      {@() dfdet_lin_icdf(0.05,-3.5,DU,-5,struct('b0',Inf)), 'start value of the search'}, ...
      {@() dfdet_adfmax_icdf(0.05,DU,-5,struct('x0',-Inf)), 'start value of the search'}, ...
      {@() dfdet_cdf(-2,DT,-5,struct('gls',-Inf)), 'needs opts.gls = cbar < 0'}, ...
      {@() dfdet_lin_cdf(-0.26,Inf,DU,-5), 'a must be a finite'}, {@() dfdet_lin_pdf(-0.26,-Inf,DU,-5), 'a must be a finite'}, ...
      {@() dfdet_region(struct('type','line','a',NaN,'b',-0.26),DU,-5), 'slope a of the line'}, ...
      {@() dfdet_local_power([-5 -5],DU,struct('alpha',NaN)), 'the level opts.alpha must lie in (0, 1)'}, ...
      {@() dfdet_local_power([-5 -5],DU,struct('c0',Inf)), 'c0 must be a finite'}, ...
      {@() pt_cdf(4,DT,-13.5,Inf), 'c must be a finite'}, {@() pt_cdf(4,DT,-Inf,-5), 'needs cbar < 0'}, ...
      {@() pt_cdf(4,DT,-13.5,-5,struct('umax',Inf)), 'umax must be'}, {@() msb_cdf(0.2,DT,NaN), 'c must be a finite'}, ...
      {@() ltratio_lincdf(-0.26,Inf,'DK',-5), 'a must be a finite'}, {@() ltratio_lincdf(-0.26,-3.5,'DG',-5,struct('cbar',-Inf)), 'cbar must be'}, ...
      {@() ltratio_region(struct('type','line','a',Inf,'b',-0.26),'DK',-5), 'slope a of the line'}};
for k = 1:numel(NF)
    msg = 'no error'; t0 = tic;
    try, NF{k}{1}(); catch ME, msg = ME.message; end
    check(sprintf('non-finite: %s refused (%.60s)',func2str(NF{k}{1}),msg),double(~contains(msg,NF{k}{2})) + double(toc(t0) > 10),0);
end
v = dfdet_power_envelope([-Inf NaN Inf],DU); w = pt_envelope([-Inf NaN Inf],DT); u = power_envelope([-Inf NaN Inf],'DK');
check('non-finite: the envelopes at c = -Inf, NaN, Inf are 1, NaN, the level',double(~(isequaln(v,[1 NaN 0.05]) && isequaln(w,[1 NaN 0.05]) && isequaln(u,[1 NaN 0.05]))),0);
r = [ltratio_region(struct('type','line','a',-3.5,'b',Inf),'DK',-5) ltratio_region(struct('type','line','a',-3.5,'b',-Inf),'DK',-5) ...
     ltratio_region(struct('type','line','a',-3.5,'b',NaN),'DK',-5) dfdet_region(struct('type','line','a',-3.5,'b',Inf),DU,-5)];
check('non-finite: the region below the line at b = Inf, -Inf, NaN is 1, 0, NaN (classical and general)',double(~isequaln(r,[1 0 NaN 1])),0);
v = [dfdet_v_cdf([-Inf -1 0 Inf NaN],DU,-5) dfdet_v_pdf([-Inf -1 Inf NaN],DU,-5)];
check('non-finite: the law of V at -Inf, -1, 0, Inf, NaN (0, 0, 0, 1, NaN; density 0, 0, 0, NaN)',double(~isequaln(v,[0 0 0 1 NaN 0 0 0 NaN])),0);
% the rules far from the unit root from c = -60 on (until November 2026 from c = -100): on -100 <= c < -60 the 45-degree
% hyperbola left errors up to 1.4e-9 against the line integral, the angle of the fixed start's rule and the caps of the
% segments from c = -60 on give 5e-13 (case DK 3e-11)
for cs = {{'DU',dfdet_trend(0),'zero'},{'DT',dfdet_trend(1),'zero'},{'SU',dfdet_trend(0),'stationary'}}
    for c = [-70 -100]
        xs = ltratio_icdf([0.01 0.5 0.99],cs{1}{1},c,struct('engine','classical'));
        Fl = ltratio_cdf(xs,cs{1}{1},c,struct('engine','classical','method','single'));
        Fg = dfdet_cdf(xs,cs{1}{2},c,struct('start',cs{1}{3}));
        check(sprintf('the rules far from the unit root: %s at c = %g, the general engine against the line integral',cs{1}{1},c),max(abs(Fg - Fl)),1e-10);
    end
end
% above the unit root (November 2026): the segments of the U-shapes shorter than 1e-12 got no panel (case DK at c = 20 gave
% 0.3128 for the 30% point), the (1,1) entry of the matrix of the transform carried an absolute rounding error of 1e-16 (the
% moments with c0 ~= 0 off far above the unit root), the circle of the moments with c0 ~= 0 continued the logarithm of the
% determinant alone (case DT at c = 0, c0 = -2 refused); the range of c is 20 without deterministic terms, 16 with them
F = dfdet_coef_cdf(19.999999940099205,dfdet_trend(-1),20);
check('above the unit root: the coefficient statistic of case DK at c = 20, the 30% point',abs(F - 0.3),2e-8);
for cs = {{'DU',dfdet_trend(0)},{'DT',dfdet_trend(1)}}
    y = coef_icdf([0.01 0.5 0.99],cs{1}{1},16,struct('engine','classical'));
    Fc = coef_cdf(y,cs{1}{1},16,struct('engine','classical')); Fg = dfdet_coef_cdf(y,cs{1}{2},16);
    check(sprintf('above the unit root: the coefficient statistic of case %s at c = 16 against the classical programs',cs{1}{1}),max(abs(Fg - Fc)),2e-8);
end
mc = ltratio_mom(4,'DT',0,struct('c0',-2,'engine','classical')); mg = dfdet_mom(4,dfdet_trend(1),0,struct('c0',-2));
check('the moments of case DT at c = 0, c0 = -2 (the circle refused; relative)',max(abs(mg - mc)./abs(mc)),1e-9);
mc = ltratio_mom(2,'DK',15,struct('c0',2,'engine','classical')); mg = dfdet_mom(2,dfdet_trend(-1),15,struct('c0',2));
check('the moments of case DK at c = 15, c0 = 2 (the (1,1) entry; relative)',max(abs(mg - mc)./abs(mc)),1e-9);
R2 = {{@() dfdet_coef_cdf(17,dfdet_trend(0),17), 'the coefficient statistic (dfdet_coef_cdf) at c = 17'}, ...
      {@() dfdet_coef_cdf(21,dfdet_trend(-1),21), 'up to c = 20'}, {@() dfdet_cdf(0,dfdet_perron('A',0.5),18,struct('c0',2)), 'the t-ratio (dfdet_cdf) at c = 18'}, ...
      {@() dfdet_mom(2,dfdet_trend(1),17), 'up to c = 16'}, {@() ltratio_cdf(0,'DU',17,struct('c0',17,'engine','general')), 'up to c = 16'}};
for k = 1:numel(R2)
    msg = 'no error';
    try, R2{k}{1}(); catch ME, msg = ME.message; end
    check(sprintf('the range of c: %s refused (%.60s)',func2str(R2{k}{1}),msg),double(~contains(msg,R2{k}{2})),0);
end
check('the range of c: the rule of the MATLAB programs takes the classical engine for case DK at c = 21 (the law of V)',double(ltratio_engine('lincdf',ltratio_spec('DK'),21,0,NaN)),0);
% the edge of the LM statistic's transform: D_V > 0 below t_* and < 0 just beyond (the estimate of the orthogonal projection
% put it up to 2.3 times too far out)
for c = [-20 2 5]
    Mo = dfdet_model(dfdet_prep(dfdet_trend(1),'lm')); ts = dfdet_tstar(Mo,c);
    dv = @(t) real(exp(ldv_(Mo,c,t)));
    ok = all(arrayfun(dv,linspace(0.01,0.99,60)*ts) > 0) && dv(1.01*ts) < 0;
    check(sprintf('LM statistic, linear trend, c = %g: D_V > 0 on (0, t_*) and < 0 at 1.01 t_*',c),double(~ok),0);
end
% the checked computation: the variant at the default angle takes the place of a steeper variant that fails with any error
% (not only those of the engine, dfdet:*), and an exhausted memory reaches the caller
CPv0 = dfdet_cpath(dfdet_startmodel(dfdet_trend(0),-5,'zero',0,[]),-5);
v = dfdet_verify(@(CPv,o) vprobe_(CPv,CPv0,'dfdet_test:other'),CPv0,1e-9);
check('verified: the fallback follows an error that is not of the engine',abs(v - 0.25),0);
id = '';
try, dfdet_verify(@(CPv,o) vprobe_(CPv,CPv0,'MATLAB:nomem'),CPv0,1e-9); catch ME, id = ME.identifier; end
check('verified: an exhausted memory (MATLAB:nomem) reaches the caller',double(~strcmp(id,'MATLAB:nomem')),0);
% the cases DK, DU, DT of the design and of the case string give the same values (a design that is its classical case)
% (with the folder of the MEX file ahead on the path the design goes to the C router and the case string to the classical
% programs on the C engine, whose choices of representation differ at some points: the same values to 1e-10)
check('a classical design is its case (DU, c = -5)',max(abs(ltratio_cdf([-3 -2],dfdet_trend(0),-5) - ltratio_cdf([-3 -2],'DU',-5))),1e-10*ltratio_croute('on'));

fprintf('\nlargest deviation by group (number of comparisons):\n');
fn = keys(G);
for k = 1:numel(fn), gv = G(fn{k}); fprintf('  %-78s %5d  %.2e\n',fn{k},gv(1),gv(2)); end
fprintf('\n%d of %d checks failed (%.0f s).\n',nfail,ntest,toc(t00));
if strict && nfail > 0, error('dfdet_test:fail','%d check(s) failed.',nfail); end

    function check(name,err,tol)
        ntest = ntest + 1; pass = ~isnan(err) && err <= tol;
        if ~pass, nfail = nfail + 1; end
        mk = {'FAIL','pass'};
        fprintf('  %-100s %9.2e  (tol %7.1e)  %s\n',name,err,tol,mk{pass+1});
    end
end

% ---------------------------------------------------------------------------------------------------------------------
function g = group_(func)
switch func
    case {'prep','dv'}, g = '1. designs: the functionals, the roots and their signs, the eigenvalues of A';
    case {'dval','roots','tstar'}, g = '2. the determinant, the roots, t_*';
    case {'cdf','pdf','err'}, g = '3. the distribution function and the density';
    case 'icdf', g = '4. quantiles';
    case 'mom', g = '5. moments';
    case {'coefcdf','coefpdf','coeficdf','lincdf','linpdf','linicdf','vcdf','vpdf','adfcdf','adfpdf','adficdf','env','lp','region'}
        g = '6. the related statistics';
    case {'ptq','ptcdf','ptpdf','ptenv','ptpow','lmq','lmp','lmcdf','lmtstar','mzacdf','mzaq','msbcdf','msbq'}, g = '7. the companions';
    otherwise, g = '8. the router';
end
end

function [v,tol,rel] = compute_(func,dtok,c,o,args)
tol = 1e-8; rel = false;
st = getf(o,'start','zero'); g0 = getf(o,'g0',0); gls = getf(o,'gls',[]);
go = struct('c0',getf(o,'c0',0),'start',st,'g0',g0,'gls',gls);
if any(strcmp(dtok,{'DK','DU','DT','SK','SU','ST','DG'})), D = dtok; else, D = design_(dtok); end
eng = getf(o,'engine','');
ro = struct; if ~isempty(eng), ro.engine = eng; end
for f = {'c0','g0','start','gls'}, if isfield(o,f{1}), ro.(f{1}) = o.(f{1}); end, end
if isfield(o,'method'), go.method = o.method; ro.method = o.method; end   % (the general engine's method: 'auto', 'paths'; 'line' refused)
switch func
    case 'prep'
        if isempty(gls), P = dfdet_prep(D); else, P = dfdet_prep(D,'gls',gls); end
        v = [P.N P.r P.rp P.rm P.nj P.p]; tol = 0;
    case 'dv'
        P = dfdet_prep(D); v = sort(P.Dv); rel = true;
    case 'dval'
        Mo = dfdet_startmodel(D,c,st,g0,gls); d = dfdet_dval(Mo,c,complex(args(1),args(2)),complex(args(3),args(4)));
        v = [real(d) imag(d)];
    case 'roots'
        Mo = dfdet_startmodel(D,c,st,g0,gls); [~,~,mu] = dfdet_mgf(Mo,c,complex(args(1),args(2))); v = 1./mu;
    case 'tstar'
        Mo = dfdet_startmodel(D,c,st,g0,gls); v = dfdet_tstar(Mo,c); tol = 1e-10; rel = true;
    case 'cdf', v = dfdet_cdf(args,D,c,go);
    case 'pdf', v = dfdet_pdf(args,D,c,go);
    case 'icdf', v = dfdet_icdf(args,D,c,go); tol = 1e-6;
    case 'mom', v = dfdet_mom(4,D,c,go); rel = true;
    case 'coefcdf', v = dfdet_coef_cdf(args,D,c,go);
    case 'coefpdf', v = dfdet_coef_pdf(args,D,c,go);
    case 'coeficdf'
        if strcmp(dtok,'pC0.3'), go.y0 = -25; end
        v = dfdet_coef_icdf(args,D,c,go); tol = 1e-6;
    case 'lincdf', v = dfdet_lin_cdf(args,o.a,D,c,go);
    case 'linpdf', v = dfdet_lin_pdf(args,o.a,D,c,go);
    case 'linicdf', v = dfdet_lin_icdf(args,o.a,D,c,go); tol = 1e-6;
    case 'vcdf', v = dfdet_v_cdf(args,D,c,go);
    case 'vpdf', v = dfdet_v_pdf(args,D,c,go);
    case 'adfcdf', v = dfdet_adfmax_cdf(args,D,c,go);
    case 'adfpdf', v = dfdet_adfmax_pdf(args,D,c,go);
    case 'adficdf', v = dfdet_adfmax_icdf(args,D,c,go); tol = 1e-6;
    case 'env', v = dfdet_power_envelope(args,D);
    case 'lp', v = dfdet_local_power(args,D);
    case 'region', v = dfdet_region(bnd_(o),D,c);
    case 'ptq', v = pt_quantile(args,D,o.cbar); tol = 1e-6;
    case 'ptcdf', v = pt_cdf(args,D,o.cbar,c);
    case 'ptpdf', v = pt_pdf(args,D,o.cbar,c);
    case 'ptenv', v = pt_envelope(args,D);
    case 'ptpow', [b,P] = pt_power(args,D,o.cbar); v = [b P]; tol = 1e-6;
    case 'lmq', v = lm_quantile(args,D); tol = 1e-6;
    case 'lmp', [~,v] = lm_power(args,D);
    case 'lmcdf', v = lm_cdf(args,D,c);
    case 'lmtstar'
        P = dfdet_prep(D,'lm'); Mo = dfdet_model(P); v = dfdet_tstar(Mo,c); tol = 1e-10; rel = true;
    case 'err'
        % a call that must fail with an error whose message contains o.txt ('_' for ' '): v = 0 if it does, Inf otherwise
        txt = strrep(o.txt,'_',' '); v = Inf; tol = 0;
        try
            switch o.fn, case 'cdf', dfdet_cdf(args,D,c,go); otherwise, error('dfdet_test: unknown err function %s.',o.fn); end
        catch ME
            if contains(ME.message,txt), v = 0; end
        end
    case 'mzacdf', v = mza_cdf(args,D,c,struct('gls',gls));
    case 'mzaq', v = mza_quantile(args,D,c,struct('gls',gls,'y0',-25)); tol = 1e-6;
    case 'msbcdf', v = msb_cdf(args,D,c,struct('gls',gls));
    case 'msbq', v = msb_quantile(args,D,c,struct('gls',gls)); tol = 1e-6;
    case 'rcdf', v = ltratio_cdf(args,D,c,ro);
    case 'rpdf', v = ltratio_pdf(args,D,c,ro);
    case 'ricdf', ro.tol = 1e-12; v = ltratio_icdf(args,D,c,ro); tol = 1e-6;
    case 'rmom', v = ltratio_mom(4,D,c,ro); rel = true;
    case 'rcoefcdf', v = coef_cdf(args,D,c,ro);
    case 'rlincdf', a = o.a; ro = rmf(ro,'a'); v = ltratio_lincdf(args,a,D,c,ro);
    case 'radfcdf', v = adfmax_cdf(args,D,c,ro);
    case 'renv', v = power_envelope(args,D,ro);
    case 'rregion', v = ltratio_region(bnd_(o),D,c,rmf(rmf(ro,'a'),'b'));
    otherwise, error('dfdet_test: unknown function %s.',func);
end
end

function l = ldv_(Mo,c,t)
[~,l] = dfdet_mgf(Mo,c,complex(t),true);
end

function v = vprobe_(CPv,CP,id)
% an evaluation of the checked computation that fails with the error id on the steeper variant
if CPv.phi > CP.phi + 1e-9, error(id,'a failure of the steeper variant (a test)'); end
v = 0.25;
end

function b = bnd_(o)
if isfield(o,'x'), b = struct('type','parabola','x',o.x); else, b = struct('type','line','a',o.a,'b',o.b); end
end

function D = design_(tok)
% the designs of validation/det_ref.jl
switch tok
    case 't-1', D = dfdet_trend(-1);
    case 'fourier', D = dfdet_fdesign(@(r) [1, sin(2*pi*r), cos(2*pi*r)]);
    case 'lt3', D = dfdet_ltrend(3);
    case 'lt8', D = dfdet_ltrend(8);
    case 'dfm'
        P = dfdet_perron('C',0.3); n = 500; X = zeros(n,4);
        for t = 1:n, r = t/n; j = 1 + (r > 0.3); X(t,:) = dfdet_peval(P.C{j},r); end
        D = dfdet_design_from_matrix(X);
    case 'dfmb'
        n = 400; t = (1:n)';
        D = dfdet_design_from_matrix([ones(n,1), t/n, double(t > 120)],struct('breaks',120));
    case 'df', D = dfdet_detrend_first(dfdet_perron('B',0.5));
    case 'lev', D = dfdet_determ([0 0.4 1],{[1 1; 1 0],[1 0; 1 0]},'1 + r, 1{r < 0.4}');
    case 'fdC', D = dfdet_fdesign(@(r) [1, r, double(r > 0.4), max(r - 0.4,0)],0.4);
    otherwise
        if tok(1) == 't', D = dfdet_trend(str2double(tok(2:end))); return, end
        if tok(1) == 'p', D = dfdet_perron(tok(2),str2double(strsplit(tok(3:end),','))); return, end
        error('dfdet_test: unknown design %s.',tok);
end
end

function o = opts_(s)
o = struct;
if strcmp(s,'-'), return, end
kv = strsplit(s,',');
for k = 1:numel(kv)
    e = strsplit(kv{k},'='); v = str2double(e{2});
    if isnan(v) && ~strcmp(e{2},'NaN'), o.(e{1}) = e{2}; else, o.(e{1}) = v; end
end
end

function v = num_(s)
s = strtrim(s);
if isempty(s), v = zeros(1,0); else, v = str2double(strsplit(s,' ')); end
end

function o = rmf(o,f)
if isfield(o,f), o = rmfield(o,f); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
