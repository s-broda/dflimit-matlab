function v = dfdet_stat(stat,mode,x,D,c,opts,a)
%DFDET_STAT  The engine dispatch of the related statistics of the general engine.
%
%   v = DFDET_STAT(stat,mode,x,D,c,opts,a): stat 'coef' (the coefficient
%   statistic n(phi^ - phi0), whose limit is (U - c0 V)/V), 'lin' (the linear
%   form U - a V, a = a; a = NaN the law of V), 'v' (the law of V) or 'adfmax'
%   (the larger of the forward and the reversed t-ratio), mode 'cdf', 'pdf' or
%   'icdf' (x then the probabilities).  Non-finite points as in DFDET_CDF.
%
%   coef: P[(U - c0 V)/V <= y] = P[U - (c0 + y) V <= 0], the distribution
%   function of the t-ratio at x = 0 with the local null value c0 + y; the
%   density by the kernel 2 v^3 (the derivative in c0); the quantile by
%   Newton's method from opts.y0 (-10), steps limited to 5.
%   lin: with R = U + a0 (a0 the shift of Mt) the event is R <= X(sqrt V),
%   X(v) = a v^2 + a0 + b, a parabola term (DFDET_XSPEC) with the constant
%   a0 + b; the quantile by the safeguarded Newton's method from opts.b0 (0).
%   v: P(V > y) = (1/2 pi i) int_C M_V(t) e^{-t y} dt/t, M_V from the roots and
%   log P(t) along the contour, each logarithm continued from node to node, on
%   the hyperbola even where the t-ratio takes the vertical line.
%   adfmax: {ADFmax <= x} = {|U + 1/2| <= 1/2 + x sqrt(V)}, two parabola terms
%   R <= x v + a0 minus R < -x v + a0 - 1 on sqrt V <= -1/(2x) for x < 0, whose
%   boundary terms cancel.
%
%   See also DFDET_COEF_CDF, DFDET_LIN_CDF, DFDET_V_CDF, DFDET_ADFMAX_CDF.
if nargin < 5 || isempty(c), c = 0; end
if nargin < 6 || isempty(opts), opts = struct; end
names = struct('coef','the coefficient statistic','lin','the linear form U - aV','v','the law of V','adfmax','ADFmax');
if strcmp(stat,'lin') && nargin >= 7 && isnan(a), nm = names.v; else, nm = names.(stat); end
[c,c0,start,g0,gls,o] = dfdet_args(D,c,opts,sprintf('%s (dfdet_%s_%s)',nm,stat,mode));
if ~(isnumeric(x) && isreal(x)), error('dfdet:args','dfdet: the argument must be real.'); end
x = double(x);
if strcmp(stat,'lin') && isnan(a), stat = 'v'; end
if strcmp(stat,'adfmax')
    % (ADFmax: least-squares detrending at the null value phi = 1; GLS detrending and a local null value were dropped
    % without notice)
    if ~isempty(gls), error('dfdet:args','dfdet_adfmax: ADFmax has least-squares detrending only (opts.gls must be empty).'); end
    if c0 ~= 0, error('dfdet:args','dfdet_adfmax: ADFmax is the statistic of the null value phi = 1 (opts.c0 must be 0).'); end
end
Mo = dfdet_startmodel(D,c,start,g0,gls);
po = struct('phi',getf(o,'phi',[]),'hfac',getf(o,'hfac',1));
CP = dfdet_cpath(Mo,c,po);
verified = logical(getf(o,'verified',false)); vtol = getf(o,'vtol',1e-9);
co = struct;
for f = {'n','ntail','alpha','filon','warn'}, if isfield(o,f{1}), co.(f{1}) = o.(f{1}); end, end
if strcmp(mode,'icdf')
    if ~all(x(:) > 0 & x(:) < 1), error('dfdet:args','dfdet: p must lie in (0,1).'); end
    v = zeros(size(x));
    for k = 1:numel(x)
        switch stat
            case 'coef'
                cc = co; cc.vpow = 3;
                FF = @(y) cpair(CP,merge(cc,struct('c0',c0 + y)));
                v(k) = dfdet_newton(FF,x(k),getf(opts,'y0',-10),5,getf(opts,'tol',1e-9),100,40);
            case 'lin'
                A0 = Mo.a;
                FF = @(b) spair(dfdet_contour('refine',CP,zeros(1,0),abs(A0 + b)),{xspec(a,0,A0 + b,struct('dA',1))},co);
                v(k) = dfdet_newton(FF,x(k),getf(opts,'b0',0),1,getf(opts,'tol',1e-10));
            case 'adfmax'
                FF = @(y) spair(dfdet_contour('refine',CP,y),{adfmax_specs(y,Mo.a)},co);
                v(k) = dfdet_newton(FF,x(k),getf(opts,'x0',-2),1,getf(opts,'tol',1e-10));
            otherwise
                error('dfdet:args','dfdet: no quantile for %s.',stat);
        end
    end
    return
end
cdf = strcmp(mode,'cdf');
v = zeros(size(x)); fin = isfinite(x);
v(isnan(x)) = NaN; if cdf, v(x == Inf) = 1; end
if strcmp(stat,'v'), fin = fin & x > 0; end            % (V > 0: 0 at y <= 0, where e^{-t y} grows along the path and the sum gave NaN)
if ~any(fin(:)), return, end
xf = x(fin); xf = xf(:).';
switch stat
    case 'v'
        CPv = CP;
        if isempty(po.phi) && isnan(CP.h)
            % (the line of the t-ratio does not resolve e^{-t y}: the hyperbola at the angle the path would have without
            % the steepening for the roots)
            if strcmp(Mo.start,'fixed') || c < -60, ph = dfdet_cpath(Mo,c,struct('fixedphi',true,'ts',CP.tstar,'tau',CP.tau,'hfac',po.hfac)); else, ph = pi/4; end
            CPv = dfdet_cpath(Mo,c,struct('phi',ph,'hfac',po.hfac));
        end
        vv = zeros(size(xf));
        for k = 1:numel(xf)
            vv(k) = vsum(CPv,xf(k),cdf);
            if cdf, vv(k) = 1 - vv(k); end
        end
    case 'coef'
        if cdf, which = 1; cc = co; else, which = 2; cc = co; cc.vpow = 3; end
        evf = @(CPv,oo) arrayfun(@(y) pick(dfdet_contour_tr(CPv,merge(merge(cc,oo),struct('c0',c0 + y))),which),xf);
        if verified, vv = dfdet_verify(evf,CP,vtol); else, vv = evf(CP,struct); end
    case 'lin'
        specs = cell(1,numel(xf));
        for k = 1:numel(xf), specs{k} = xspec(a,0,Mo.a + xf(k),struct('dA',1)); end
        vv = spec_eval(CP,specs,1 + ~cdf,verified,vtol,co);
    case 'adfmax'
        specs = cell(1,numel(xf));
        for k = 1:numel(xf), specs{k} = adfmax_specs(xf(k),Mo.a); end
        vv = spec_eval(CP,specs,1 + ~cdf,verified,vtol,co);
end
v(fin) = vv;
end

function r = dfdet_contour_tr(CP,o)
[F,f] = dfdet_contour('tratio',CP,0,o); r = [F f];
end

function v = pick(r,k)
v = r(k);
end

function r = cpair(CP,o)
[F,f] = dfdet_contour('tratio',CP,0,o); r = [F f];
end

function r = spair(CP,specs,o)
[F,f] = dfdet_contour('specs',CP,specs,o); r = [F f];
end

function v = spec_eval(CP,specs,which,verified,vtol,co)
% the values (which = 1) or the derivatives (2) of the events specs on the outer contour CP, the line refined to their
% slopes and constants
qm = zeros(1,numel(specs)); am = [];
for i = 1:numel(specs)
    qm(i) = max([0 abs([specs{i}.q])]); am = [am abs([specs{i}.A])]; %#ok<AGROW>
end
evf = @(CPv,oo) pickv(dfdet_contour_s(dfdet_contour('refine',CPv,qm,am),specs,merge(co,oo)),which);
if verified, v = dfdet_verify(evf,CP,vtol); else, v = evf(CP,struct); end
end

function r = dfdet_contour_s(CP,specs,o)
[F,f] = dfdet_contour('specs',CP,specs,o); r = {F,f};
end

function v = pickv(r,k)
v = r{k};
end

function S = xspec(c0,q,A,o)
S = dfdet_xspec(c0,q,A,o);
end

function S = adfmax_specs(x,a0)
% the two parabola terms of {ADFmax <= x}: R <= x v + a0 minus R < -x v + a0 - 1, on sqrt V <= -1/(2x) for x < 0, where
% both parabolas reach the same value, so that their boundary terms cancel exactly and are left out
if x < 0, vm = -0.5/x; else, vm = Inf; end
S = [dfdet_xspec(0,x,a0,struct('vmax',vm,'w',1,'dq',1,'bnd',false)), dfdet_xspec(0,-x,a0 - 1,struct('vmax',vm,'w',-1,'dq',-1,'bnd',false))];
end

function acc = vsum(CP,y,cdf)
% (1/pi) sum' Im{M_V(t_k) e^{-t_k y} t'(u_k)/t_k} (P(V > y); without the 1/t_k the density): log M_V(t) = log P(t) -
% (1/2) sum_k log(+-s_k(t)) at s = 0, each logarithm continued from node to node
Mo = CP.Mo; r = size(CP.R,1); rp = Mo.rp;
Z = CP.R; Z(rp+1:r,:) = -Z(rp+1:r,:);
lg = log(Z);
for k = 2:size(Z,2)
    lg(:,k) = lg(:,k) + 2i*pi*round((imag(lg(:,k-1)) - imag(lg(:,k)))/(2*pi));
end
lM = CP.lP(:) - sum(lg,1).'/2;
if strcmp(Mo.start,'fixed'), lM = lM + Mo.g0^2*CP.fx.cf(:); end
if cdf, d = CP.dt./CP.t; else, d = CP.dt; end
acc = sum(CP.wt.*imag(exp(lM - CP.t*y).*d));
end

function o = merge(o,p)
for f = fieldnames(p)', o.(f{1}) = p.(f{1}); end
end

function v = getf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
