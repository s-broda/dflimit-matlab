function nfail = ext_test(strict)
%EXT_TEST  Validation of the extensions of the limit against the Julia engine DFLimit.jl (and the MEX file).
%
%   nfail = EXT_TEST compares the extensions of the general engine of these
%   programs (VOLPATH and opts.vol, opts.rho, CADF_C0, PREDICTIVE_C0,
%   HEGY_T3_CDF, HEGY_T3_PDF, HEGY_T3_QUANTILE) with the reference values of
%   the Julia engine in validation/ext_ref.txt (written by dev/ext_ref.jl of
%   DFLimit.jl; the tags in its header): the volatility transform D(s,t), its
%   roots and t_*, the distribution functions, densities and quantiles with
%   volatility regimes, with a normal component (both routes, the predictive
%   regression) and of the first seasonal t-ratio, with the tolerance of each
%   line (1e-9 for distribution functions and densities, 1e-10 for the
%   transforms and t3, 1e-7 for quantiles); the argument errors of the file;
%   the limiting cases (one regime is the engine itself, a dummy regime
%   boundary changes nothing, the transform at a constant variance is the
%   engine's, the transfer matrices without deterministic terms, the endpoints
%   and the reflection of rho); and, where the MEX file ltratio_mex is found
%   (the folder c/mex of the repository is appended to the path for the
%   test), the MATLAB programs against the C engine.  It prints each group
%   with its count and largest deviation and returns the number of failures;
%   EXT_TEST(true) raises an error at the end when a check has failed.  (The
%   second seasonal t-ratio and the F statistic are not part of the release
%   0.5.0; their lines left the reference file on November 27, 2026.)
%   About two minutes.
%
%   See also DFDET_TEST, VOLPATH, HEGY_T3_CDF.
if nargin < 1 || isempty(strict), strict = false; end
here = fileparts(mfilename('fullpath'));
mexdir = fullfile(here,'..','c','mex');
if exist(mexdir,'dir') && exist('ltratio_mex','file') ~= 3, addpath(mexdir,'-end'); end
hasmex = exist('ltratio_mex','file') == 3;
G = containers.Map(); nfail = 0; t00 = tic;
ws = warning('off','dfdet:rounding'); wc = warning('off','dfdet:contour'); wv = warning('off','dfdet:verified');
cl = onCleanup(@() warning([ws wc wv]));
    function chk(g_,got_,want_,tol_,what_)
        % (nested: its variables end in _ so that the test's own are not shared)
        got_ = got_(:); want_ = want_(:); if isscalar(want_), want_ = repmat(want_,size(got_)); end
        dd_ = max([0; abs(got_ - want_)]);
        if numel(got_) ~= numel(want_) || any(isnan(got_) ~= isnan(want_)), dd_ = Inf; end
        if ~isKey(G,g_), G(g_) = [0 0 0]; end
        gv_ = G(g_); gv_(1) = gv_(1) + 1; gv_(2) = max(gv_(2),dd_);
        if ~(dd_ <= tol_), nfail = nfail + 1; gv_(3) = gv_(3) + 1; fprintf('FAIL %-10s %s: deviation %.3g (tolerance %.1g)\n',g_,what_,dd_,tol_); end
        G(g_) = gv_;
    end
    function chkerr(g_,f_,frag_,what_)
        msg_ = '';
        try, f_(); catch ME_, msg_ = ME_.message; end
        if ~isKey(G,g_), G(g_) = [0 0 0]; end
        gv_ = G(g_); gv_(1) = gv_(1) + 1;
        if isempty(strfind(msg_,frag_)), nfail = nfail + 1; gv_(3) = gv_(3) + 1; fprintf('FAIL %-10s %s: "%s" without "%s"\n',g_,what_,msg_,frag_); end
        G(g_) = gv_;
    end
L = strsplit(strtrim(fileread(fullfile(here,'validation','ext_ref.txt'))),newline);
Ds = {}; Vs = {};
for i = 1:numel(L)
    l = strtrim(L{i}); if isempty(l) || l(1) == '#', continue, end
    w = strsplit(l,' '); tag = w{1};
    switch tag
        case 'design'
            id = str2double(w{2});
            if strcmp(w{3},'trend'), Ds{id} = dfdet_trend(str2double(w{4})); %#ok<AGROW>
            else, q = str2double(w{5}); Ds{id} = dfdet_perron(w{4},str2double(w(6:5+q))); end %#ok<AGROW>
        case 'vol'
            id = str2double(w{2}); nb = str2double(w{3});
            Vs{id} = volpath(str2double(w(4:3+nb)),str2double(w(4+nb:4+2*nb))); %#ok<AGROW>
        case {'vdval','vroots','vtstar','vval','vquant'}
            k0 = 2; fn = '';
            if strcmp(tag,'vval'), fn = w{2}; k0 = 3; end
            D = Ds{str2double(w{k0})}; V = Vs{str2double(w{k0+1})}; nc = str2double(w{k0+2}); c = str2double(w(k0+3:k0+2+nc)); k = k0 + 3 + nc;
            num = str2double(w(k:end));
            switch tag
                case 'vdval'
                    VM = dfdet_vol('model',D,V,c);
                    d = dfdet_vol('dval',VM,complex(num(1),num(2)),complex(num(3),num(4)));
                    chk('vdval',abs(d - complex(num(5),num(6))),0,num(7),l(1:min(end,60)));
                case 'vroots'
                    VM = dfdet_vol('model',D,V,c);
                    [~,ld,mu] = dfdet_mgf(VM,0,complex(num(1),num(2)));
                    r = num(4); ref = complex(num(5:2:4+2*r),num(6:2:5+2*r));
                    chk('vroots',real(ld),num(3),num(end),l(1:min(end,60)));
                    mu = sortc(mu); ref = sortc(ref(:));
                    chk('vroots',abs(mu - ref),0,num(end)*max(abs(ref)),[l(1:min(end,60)) ' roots']);
                case 'vtstar'
                    VM = dfdet_vol('model',D,V,c);
                    chk('vtstar',dfdet_tstar(VM,min(c)),num(1),num(2),l(1:min(end,60)));
                case 'vval'
                    c0 = num(1); n = num(2); x = num(3:2+n); ref = num(4+n:3+2*n); tol = num(end);
                    o = struct('vol',V,'c0',c0);
                    if strcmp(fn,'cdf'), v = dfdet_cdf(x,D,c,o); else, v = dfdet_pdf(x,D,c,o); end
                    chk(['v' fn],v,ref,tol,l(1:min(end,60)));
                    if hasmex, chk(['mexv' fn],v,mexdet(fn,D,x,c(1),[c0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN numel(V.om) V.lam(2:end-1) V.om c]),1e-12,l(1:min(end,60))); end
                case 'vquant'
                    q = dfdet_icdf(num(2),D,c,struct('vol',V,'c0',num(1)));
                    chk('vquant',q,num(3),num(4),l(1:min(end,60)));
            end
        case {'cval','cquant'}
            if strcmp(tag,'cval')
                fn = w{2}; D = Ds{str2double(w{3})}; st = str2double(w{4}); num = str2double(w([5:8 10:end]));
                g0 = num(1); c = num(2); c0 = num(3); rho = num(4); m = w{9};
                n = num(5); x = num(6:5+n); ref = num(7+n:6+2*n); tol = num(end);
                stn = {'zero','fixed','stationary'};
                o = struct('c0',c0,'start',stn{st + 1},'g0',g0,'rho',rho);
                if rho ~= 0, o.method = m; end
                if strcmp(fn,'cdf'), v = dfdet_cdf(x,D,c,o); else, v = dfdet_pdf(x,D,c,o); end
                chk(['c' fn],v,ref,tol,l(1:min(end,70)));
                if hasmex
                    mc = 0; if rho ~= 0, mc = 3 + strcmp(m,'mixture'); end
                    chk(['mexc' fn],v,mexdet(fn,D,x,c,[c0 st g0 NaN 0 mc NaN NaN NaN NaN rho]),1e-11,l(1:min(end,70)));
                end
            else
                D = Ds{str2double(w{2})}; num = str2double(w(4:end));
                q = dfdet_icdf(num(5),D,num(2),struct('c0',num(3),'rho',num(4)));
                chk('cquant',q,num(6),num(7),l(1:min(end,60)));
            end
        case 'sval'
            fn = w{2}; if ~strncmp(fn,'t3',2), continue, end
            D = Ds{str2double(w{3})}; num = str2double(w(4:end)); c = num(1); c0 = num(2); vid = num(3);
            n = num(4); x = num(5:4+n); ref = num(6+n:5+2*n); tol = num(end);
            o = struct('c0',c0); cc = c; if vid > 0, o.vol = Vs{vid}; cc = [c c]; end
            if strcmp(fn,'t3cdf'), v = hegy_t3_cdf(x,D,cc,o); else, v = hegy_t3_pdf(x,D,cc,o); end
            chk(fn,v,ref,tol,l(1:min(end,60)));
            if hasmex
                ov = [c0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]; if vid > 0, ov = [ov numel(Vs{vid}.om) Vs{vid}.lam(2:end-1) Vs{vid}.om cc]; end
                chk(['mex' fn],v,mexdet(fn,D,x,c,ov),1e-12,l(1:min(end,60)));
            end
        case 'squant'
            if ~strcmp(w{2},'t3'), continue, end
            D = Ds{str2double(w{3})}; num = str2double(w(4:end));
            chk('t3quant',hegy_t3_quantile(num(2),D,num(1)),num(3),num(4),l(1:min(end,60)));
    end
end
% the argument errors of the reference file (the fragments of the messages)
T0 = dfdet_trend(0); T1 = dfdet_trend(1); V2 = volpath(0.3,[1 3]);
chkerr('errors',@() volpath(1.2),'must lie in (0, 1)','volpath, a break date at 1.2');
chkerr('errors',@() volpath([0.5 0.3]),'strictly increasing','volpath, break dates 0.5, 0.3');
chkerr('errors',@() volpath(0.5,1),'one standard deviation omega per regime','volpath, one omega for two regimes');
chkerr('errors',@() volpath(0.5,[1 -1]),'finite and positive','volpath, omega -1');
chkerr('errors',@() dfdet_cdf(-2,T0,2,struct('vol',V2)),'must be at most 0','dfdet_cdf vol, c = 2');
chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('vol',V2,'start','stationary')),'zero start only','dfdet_cdf vol, stationary start');
chkerr('errors',@() dfdet_cdf(-2,T1,-5,struct('vol',V2,'gls',-13.5)),'GLS detrending','dfdet_cdf vol, gls = -13.5');
chkerr('errors',@() dfdet_cdf(-2,T0,[0 -1 0],struct('vol',V2)),'one value per regime','dfdet_cdf vol, three c for two regimes');
chkerr('errors',@() dfdet_cdf(-2,T0,[0 -1]),'needs a volatility path','dfdet_cdf c a vector without vol');
chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('vol',V2,'rho',0.5)),'rho and vol together','dfdet_cdf vol and rho = 0.5');
chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('rho',1.5)),'rho must lie in [-1, 1]','dfdet_cdf rho = 1.5');
chkerr('errors',@() dfdet_cdf(-2,T1,-5,struct('rho',0.5,'gls',-13.5)),'GLS detrending','dfdet_cdf rho = 0.5, gls = -13.5');
chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('rho',0.5,'g0',1,'method','lines')),'zero or a stationary start','dfdet_cdf rho = 0.5, lines, fixed start');
% (the methods of rho without it, rho = 1: refused with the message of DFLimit.jl for a method other than 'auto', 'line'
% and 'paths'; in the MEX file the message of the C engine, the same words; November 26, 2026)
LM_ = 'method must be ''auto'', ''line'' (the line integral with the explicit kernel';
for m_ = {'lines','mixture'}
    chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('method',m_{1})),LM_,['dfdet_cdf rho = 1, method ' m_{1}]);
    chkerr('errors',@() dfdet_pdf(-2,T0,-5,struct('method',m_{1})),LM_,['dfdet_pdf rho = 1, method ' m_{1}]);
    chkerr('errors',@() dfdet_icdf(0.05,T0,-5,struct('method',m_{1})),LM_,['dfdet_icdf rho = 1, method ' m_{1}]);
    if hasmex
        d_ = ltratio_croute('design',T0);
        chkerr('errors',@() ltratio_mex('det','cdf',d_{:},-2,-5,[0 0 0 NaN 0 3 + strcmp(m_{1},'mixture')]), ...
               'method must be ''auto'', ''line'' (the line integral with the explicit kernel) or ''paths'' (the inversion on paths around the roots)',['the MEX file, dist_cdf rho = 1, method ' m_{1}]);
    end
end
chkerr('errors',@() cadf_c0(-5,0),'not be 0','cadf_c0 rho = 0');
chkerr('errors',@() hegy_t3_cdf(1,T0,1),'must be at most 0','hegy_t3_cdf c = 1');
chkerr('errors',@() dfdet_cdf(-2,T0,-5,struct('vol',V2,'method','line')),'method must be ''auto'' or ''paths''','dfdet_cdf vol, method line');
chkerr('errors',@() ltratio_cdf(-2,'DU',0,struct('vol',V2,'engine','classical')),'engine = ''classical''','ltratio_cdf vol, engine classical');
chkerr('errors',@() ltratio_cdf(-2,'DU',[0 -1]),'needs a volatility path','ltratio_cdf c a vector without vol');
% the limiting cases
xs = [-4 -2.86 -1 0.5];
for D = {T0, dfdet_perron('B',0.5)}
    for c = [0 -10]
        chk('vone',dfdet_cdf(xs,D{1},c,struct('vol',volpath([]))),dfdet_cdf(xs,D{1},c),0,'one regime: F');
        chk('vone',dfdet_pdf(xs,D{1},c,struct('vol',volpath([]))),dfdet_pdf(xs,D{1},c),0,'one regime: f');
    end
end
for D = {dfdet_trend(-1), T0, T1, dfdet_perron('C',0.6)}
    for c = [0 -10]
        V = volpath(0.4,[1 1]);
        chk('vdummy',dfdet_cdf(xs,D{1},c,struct('vol',V)),dfdet_cdf(xs,D{1},c,struct('method','paths')),1e-12,sprintf('a dummy boundary, c = %g: F',c));
        chk('vdummy',dfdet_pdf(xs,D{1},c,struct('vol',V)),dfdet_pdf(xs,D{1},c,struct('method','paths')),1e-12,sprintf('a dummy boundary, c = %g: f',c));
    end
end
for D = {dfdet_trend(-1), T0, T1, dfdet_perron('A',0.3)}
    for c = [0 -5]
        VM = dfdet_vol('model',D{1},volpath(0.4,[1 1]),[c c]); Mo = dfdet_startmodel(D{1},c,'zero',0,[]);
        for st = [0.3 -2; 0.2+0.7i 0.3-4i; 0.1 -50+20i].'
            a = dfdet_vol('dval',VM,st(1),st(2)); b = dfdet_dval(Mo,c,st(1),st(2));
            chk('vconst',abs(a - b)/abs(b),0,1e-10,'the transform at a constant variance');
        end
    end
end
for k = 1:2
    if k == 1, V = volpath(0.3,[1 3]); cr = [0 0]; else, V = volpath([0.1 0.3 0.6 0.8],[1 2 0.5 3 1.5]); cr = [-1 0 -4 -2 0]; end
    VM = dfdet_vol('model',dfdet_trend(-1),V,cr);
    for st = [0.3 -2; 0.2+0.7i 0.3-4i; 0 3].'
        chk('vtransfer',abs(dfdet_vol('dval',VM,st(1),st(2))/transfer_d0(V,cr,st(1),st(2)) - 1),0,1e-11,'the transfer matrices');
    end
end
x4 = [-3 -1 0.5 2];
chk('rhoends',dfdet_cdf(x4,T0,-5,struct('rho',0)),erfc(-x4/sqrt(2))/2,1e-15,'rho = 0');
chk('rhoends',dfdet_cdf(x4,T0,-5,struct('rho',-1)),1 - dfdet_cdf(-x4,T0,-5),1e-15,'rho = -1');
chk('rhoends',dfdet_cdf(x4,T0,-5,struct('rho',-0.7,'method','lines')),1 - dfdet_cdf(-x4,T0,-5,struct('rho',0.7,'method','lines')),1e-15,'rho = -0.7 by the reflection');
chk('rhoends',dfdet_cdf([-Inf Inf],T0,-5,struct('rho',-0.7)),[0 1],0,'rho = -0.7 at -Inf, Inf');
chk('rhoends',cadf_c0(-5,sqrt(0.5)),5,1e-13,'cadf_c0(-5, sqrt(0.5))');
chk('rhoends',predictive_c0(-5,-0.95,1),-5 + 1/0.95,1e-14,'predictive_c0(-5, -0.95, 1)');
h = 1e-3;
F = dfdet_cdf([-2-h -2+h],T1,[-3 0 -10],struct('vol',volpath([0.2 0.7],[2 1 0.5])));
chk('vdensity',(F(2) - F(1))/(2*h),dfdet_pdf(-2,T1,[-3 0 -10],struct('vol',volpath([0.2 0.7],[2 1 0.5]))),1e-6,'the density against the difference quotient');
q = ltratio_icdf(0.05,'DU',0,struct('vol',V2));
chk('vclassical',ltratio_cdf(q,'DU',0,struct('vol',V2)),0.05,1e-9,'ltratio_icdf with vol: F at the quantile');
chk('vclassical',ltratio_cdf([-3 -2],'DU',-5,struct('vol',V2)),dfdet_cdf([-3 -2],T0,-5,struct('vol',V2)),0,'ltratio_cdf with vol is dfdet_cdf of the design');
% the summary
fprintf('%-12s %8s %8s %12s\n','group','values','failed','largest dev');
ks = keys(G);
for k = 1:numel(ks), v = G(ks{k}); fprintf('%-12s %8d %8d %12.2e\n',ks{k},v(1),v(3),v(2)); end
fprintf('ext_test: %d failed (%.0f s; MEX file %s)\n',nfail,toc(t00),tern(hasmex,'checked','not found'));
if strict && nfail > 0, error('ext_test:fail','ext_test: %d checks failed',nfail); end
end

function v = mexdet(fn,D,x,c,o)
d = ltratio_croute('design',D);
v = ltratio_mex('det',fn,d{:},x,c,o);
end
function z = sortc(z)
[~,i] = sortrows([real(z(:)) imag(z(:))]); z = z(i);
end
function d = transfer_d0(V,cr,s,t)
% the determinant without deterministic terms by the transfer matrices of the regimes (the backward equation)
x = [1; s];
for i = numel(V.om):-1:1
    l = V.lam(i+1) - V.lam(i); c = cr(i); w2 = V.om(i)^2;
    th = sqrt(complex(c^2 - 2*t*w2));
    ch = cosh(th*l); if abs(th*l) < 1e-6, sh = l*(1 + (th*l)^2/6); else, sh = sinh(th*l)/th; end
    x = exp(c*l)*[(ch - c*sh)*x(1) - w2*sh*x(2); 2*t*sh*x(1) + (ch + c*sh)*x(2)];
end
d = x(1);
end
function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
