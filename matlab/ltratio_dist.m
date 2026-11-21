function [F,info] = ltratio_dist(x,cas,c,opts,mode)
%LTRATIO_DIST  Engine dispatch of LTRATIO_CDF and LTRATIO_PDF (mode 'cdf' or 'pdf').
%
%   Validates the arguments, chooses the representation (opts.method:
%   'auto', the default, 'single' or 'line', 'double' or 'conditional',
%   'contour'; see LTRATIO_ARGS for a fixed start), maps the case to the
%   engine's conventions and calls DFX_CORE (line integral), TR_LIMDIST
%   (conditional integral) or TR_CONTOUR (contour integral).  With 'auto' the
%   choice is made point by point by LTRATIO_METHOD, from the timings of the
%   representations in the MATLAB programs (single for c0 = 0; with
%   c0 ~= 0 double except for c0 > 0 at x >= 0 in cases DK, DU, DT, SK, and
%   single except for c0 < 0 at x >= 0 in cases SU, ST; the contour integral
%   with a fixed start in cases DU, DT); info.ndouble is
%   the number of points sent to the conditional integral by the rule.  A point
%   at which the conditional integral reports that its error estimate is not
%   met (info.status, info.err) gets the line integral's value instead, when a
%   zero or stationary start makes it available; info.nfallback counts them.
%   With a fixed start in cases DU, DT ('auto', the contour integral at every
%   point, info.ncontour of them) a point at which the contour integral
%   reports that its error estimate is not met gets the conditional
%   integral's value where that is clean (info.nfallback; otherwise the
%   flagged value stays, with info.status), and the contour values below
%   1e-3 are recomputed by the conditional integral where it is clean
%   (info.ntail), as the C and Julia engines do, both only where the
%   conditional integral is checked in these programs (c >= -20, or c >= -40
%   with (c g0)^2/2 <= 100; LTRATIO_METHOD); info.err and info.status
%   are those of the values returned.  With 'contour' the contour values are
%   returned as they are, with a warning where the status is nonzero.
%
%   Tails.  The line integral has an absolute accuracy, not a relative one:
%   its terms are of order one, and the computed density or distribution
%   function carries an error of order 1e-14 to 1e-21 whatever its size,
%   which finer rules do not reduce (2e-14 in the lower tail of cases DT
%   and ST at c0 ~= 0; densities below the floor are noise and can be
%   negative: -6e-21 at x = 8.5 in case SU, c = c0 = -1.2, where the density
%   is 3e-22).  With 'auto', the values of the line integral with a
%   density below 1e-3 in absolute value, or a distribution function below
%   1e-3, are therefore recomputed by the conditional integral, which keeps its
%   relative accuracy in the tails and returns an error estimate (info.ntail
%   of them; info.err holds the estimates at every point computed by the
%   conditional integral and NaN elsewhere), so that the relative error stays
%   near 1e-10 or below (not with a stationary start for c < -50, nor with
%   c0 ~= 0 and a zero start for c < -100, where LTRATIO_METHOD keeps the
%   line integral throughout).  With opts.method = 'single' given explicitly, the
%   line integral is used throughout.  Far out, below the error of either
%   method, a computed density can be a tiny negative number; it is returned
%   as zero.  Non-finite x: F = 0 at -Inf, 1 (cdf) or 0 (pdf) at +Inf, NaN
%   at NaN.
%
%   The bounds of the engine.  Where the line integral would need a Gauss-
%   Legendre rule of more than 4096 nodes (large x above the unit root: the
%   segment to the saddle of the Gaussian factor, from x ~ 660 to 780 at
%   c = 7) or a path of more than two million nodes (c0 < 0 near zero with
%   x > 0, where a root of X(v) = c0 v^2 + x v + 1/2 lies far out), it
%   refuses the point before evaluating it (DFX_CORE raises dfx:bound:rule or
%   dfx:bound:path; the bounds of the C and Julia engines).  With 'auto' the
%   conditional integral serves there, with its checks (at large x above the
%   unit root its one evaluation where the checked value is flagged; the
%   warning ltratio:bound where a value near c0 = 0 is flagged; ROUTE_DOUBLE),
%   info.nrouted of them; with 'single', and in case DG, which has the line
%   integral only, the error ltratio:bound names the arguments.  With
%   opts.verified = true a flagged point whose fallback, the line integral,
%   refuses it keeps its flagged value, with the warning ltratio:bound.
%
%   See also LTRATIO_METHOD, LTRATIO_CDF, LTRATIO_PDF, DFX_CDF, TR_LIMDIST.
t0 = tic;
[cas,c,c0,g0,method] = ltratio_args(cas,c,opts);
if ~(isnumeric(x) && isreal(x)), error('ltratio:args','ltratio_%s: x must be real.',mode); end
F = zeros(size(x)); fin = isfinite(x);
F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
info = struct('method',method,'time',0);
if any(fin(:))
    switch method
        case 'single'
            F(fin) = single_values(x(fin),cas,c,c0,opts,mode,'error');
        case 'double'
            if ~logical(getf_(opts,'verified',false)), ltratio_brittle(cas,c,c0,g0,{'double'}); end
            o = rmfield_safe(opts,{'method'});
            o.c0 = c0; o.g0 = g0;
            [Fv,inf2] = tr_limdist(x(fin),cas,c,o,mode);
            F(fin) = Fv;
            fn = fieldnames(inf2);
            for k = 1:numel(fn), info.(fn{k}) = inf2.(fn{k}); end
        case 'contour'
            if ~logical(getf_(opts,'verified',false)), ltratio_brittle(cas,c,c0,g0,{'contour'}); end
            o = rmfield_safe(opts,{'method'});
            o.c0 = c0; o.g0 = g0;
            [Fv,inf2] = tr_contour(x(fin),cas,c,o,mode);
            F(fin) = Fv;
            fn = fieldnames(inf2);
            for k = 1:numel(fn), info.(fn{k}) = inf2.(fn{k}); end
        case 'auto'
            verified = logical(getf_(opts,'verified',false));
            [dbl,tails,ctr] = ltratio_method(cas,c,c0,g0,x,verified); ctr = fin & ctr;
            if ~verified
                % the default: one evaluation by one method at each point (no second evaluation, no refinement driven by
                % an estimate, no fallback, no recomputation of the tails; ltratio_dist of julia/DFExact.jl)
                o = rmfield_safe(opts,{'method'}); o.c0 = c0; o.g0 = g0; o.verified = false;
                dbl = fin & dbl & ~ctr; sgl = fin & ~dbl & ~ctr;
                used = {}; if any(ctr(:)), used{end+1} = 'contour'; end; if any(dbl(:)), used{end+1} = 'double'; end
                ltratio_brittle(cas,c,c0,g0,used);
                info.ndouble = nnz(dbl); info.ncontour = nnz(ctr);
                if any(ctr(:)), F(ctr) = tr_contour(x(ctr),cas,c,o,mode); end
                if any(dbl(:)), F(dbl) = tr_limdist(x(dbl),cas,c,o,mode); end
                if any(sgl(:)), [F(sgl),info.nrouted] = single_values(x(sgl),cas,c,c0,opts,mode,'auto'); end
                if strcmp(mode,'pdf'), F(F < 0) = 0; end
                info.time = toc(t0);
                return
            end
            if any(ctr(:))                                   % a fixed start in cases DU, DT: the contour integral
                [F,info] = contour_auto(F,info,x,ctr,tails,cas,c,c0,g0,opts,mode);
                if strcmp(mode,'pdf'), F(F < 0) = 0; end
                info.time = toc(t0);
                return
            end
            dbl = fin & dbl;
            sgl = fin & ~dbl; tl = false(size(x));
            info.ndouble = nnz(dbl); info.ntail = 0;
            if any(sgl(:))
                [F(sgl),info.nrouted] = single_values(x(sgl),cas,c,c0,opts,mode,'auto');
                if tails
                    if strcmp(mode,'pdf'), tl = sgl & abs(F) < 1e-3; else, tl = sgl & F < 1e-3; end
                end
                info.ntail = nnz(tl); dbl = dbl | tl;
            end
            Fs = F;                                          % the line integral's values, kept for the fallback at tail points
            if any(dbl(:))
                o2 = rmfield_safe(opts,{'method'}); o2.c0 = c0; o2.g0 = g0;
                ws = warning('off','tr_limdist:status');
                [Fv,inf2] = tr_limdist(x(dbl),cas,c,o2,mode);
                warning(ws);
                F(dbl) = Fv;
                info.err = NaN(size(x)); if isfield(inf2,'err'), info.err(dbl) = inf2.err; end
                info.status = zeros(size(x)); if isfield(inf2,'status'), info.status(dbl) = inf2.status; end
                % a point at which the conditional integral reports that its error estimate is not met gets the line
                % integral's value (the conditional integral is the only route with a fixed start)
                fb = info.status ~= 0 & g0 == 0; info.nfallback = nnz(fb);
                F(fb & tl) = Fs(fb & tl);                    % tail points: the line integral's value was computed above
                fbr = fb & ~tl;
                if any(fbr(:)), F(fbr) = single_values(x(fbr),cas,c,c0,opts,mode,'keep',F(fbr)); end
            end
    end
end
if strcmp(mode,'pdf'), F(F < 0) = 0; end            % a density below the error of the method
info.time = toc(t0);
end

function [F,info] = contour_auto(F,info,x,ctr,tails,cas,c,c0,g0,opts,mode)
% opts.method = 'auto' with a fixed start in cases DU, DT: the contour integral at the points ctr; a point at which it
% reports that its error estimate is not met gets the conditional integral's value where that is clean (these programs
% have no line integral with a fixed start: otherwise the flagged contour value stays, with its status), and the values
% below 1e-3 (a density in absolute value, a distribution function in the lower tail) are recomputed by the conditional
% integral, whose error is relative, where it is clean (the rule of the C and Julia engines); both only where the
% conditional integral is checked in these programs (tails, from LTRATIO_METHOD)
o = rmfield_safe(opts,{'method'}); o.c0 = c0; o.g0 = g0;
ws = warning('off','tr_contour:status');
[Fc,ic] = tr_contour(x(ctr),cas,c,o,mode);
warning(ws);
F(ctr) = Fc;
info.ncontour = nnz(ctr); info.nfallback = 0; info.ntail = 0;
info.err = NaN(size(x)); info.err(ctr) = ic.err;
info.status = zeros(size(x)); info.status(ctr) = ic.status;
fb = ctr & info.status ~= 0;
if strcmp(mode,'pdf'), small = abs(F) < 1e-3; else, small = F < 1e-3; end
tl = ctr & info.status == 0 & small;
re = (fb | tl) & tails;                                % (tails: the conditional integral is checked there, LTRATIO_METHOD)
if any(re(:))
    ws = warning('off','tr_limdist:status');
    [Fd,id] = tr_limdist(x(re),cas,c,o,mode);
    warning(ws);
    ir = find(re); ok = id.status(:) == 0; ig = ir(ok);
    F(ig) = Fd(ok); info.err(ig) = id.err(ok); info.status(ig) = 0;
    info.nfallback = nnz(fb(ig)); info.ntail = nnz(tl(ig));
end
end

function [F,nrouted] = single_values(x,cas,c,c0,opts,mode,route,Fkeep)
% the line integral (DFX_CORE) at the points x; at a point that it refuses (its paths would exceed the bounds of the
% engine: the errors dfx:bound:rule and dfx:bound:path, raised before the point is evaluated), by route: 'auto' the
% conditional integral (ROUTE_DOUBLE, with the warning ltratio:bound where its values near c0 = 0 are flagged; case DG,
% which has the line integral only, an error), 'error' the error ltratio:bound that names the arguments, 'keep' the
% value Fkeep, with the warning ltratio:bound (the fallback of opts.verified = true at a flagged point)
if strcmp(cas,'DG'), model = 4; else, model = find('KUT' == cas(2)); end
o = rmfield_safe(opts,{'method','g0'});
o.c0 = c0;
if cas(1) == 'S', o.start = 'stationary'; else, o.start = 'zero'; end
F = zeros(size(x)); nrouted = 0; nflag = 0; nkeep = 0;
for k = 1:numel(x)
    try
        F(k) = dfx_core(mode,x(k),model,c,o);
    catch ME
        if ~strncmp(ME.identifier,'dfx:bound:',10), rethrow(ME); end
        if strcmp(route,'keep'), F(k) = Fkeep(k); nkeep = nkeep + 1; continue; end
        if strcmp(route,'error') || model == 4, error('ltratio:bound','%s',bound_msg(mode,cas,c,c0,x(k),ME)); end
        [F(k),fl] = route_double(x(k),cas,c,c0,opts,mode,ME.identifier);
        nrouted = nrouted + 1; nflag = nflag + fl;
    end
end
head = sprintf('%s, c = %g, c0 = %g',cas,c,c0);
if nflag > 0
    ltratio_brittle('warn','ltratio:bound',[head ': the line integral''s paths exceed the bounds of the engine at x > 0 (c0 < 0 near ' ...
        'zero), and the conditional integral that takes these points flags some of them; such values may be off by more than ' ...
        '1e-5 (the density of case SK at c = c0 = -1e-3, x = 0.05 by 1.2e-3)']);
end
if nkeep > 0
    ltratio_brittle('warn','ltratio:bound',[head ': at some points the value of the conditional integral is flagged and the line ' ...
        'integral''s paths exceed the bounds of the engine; the flagged values are returned']);
end
end

function [F,fl] = route_double(x,cas,c,c0,opts,mode,id)
% The conditional integral at a point that the line integral refuses, with opts.method = 'auto' (route_double of
% DFExact.jl and of the C engine), by the bound: a far root (dfx:bound:path, c0 < 0 near zero with x > 0) with its checks
% (opts.verified = true), whatever its status, fl where it is flagged (against the contour integral of a newer engine,
% November 2026: within 5e-7 at most of the points of case SK with c = c0 from -1e-4 to -0.1 and x from 0.05 to 1.5,
% flagged at most of them, but the density at c = c0 = -1e-3, x = 0.05 off by 1.2e-3; case DK within 1e-13); large x
% above the unit root (dfx:bound:rule) with its checks where its status is 0, otherwise one evaluation (in case DT the
% checked value is off by up to 3e-2 with the status 2 at c = 7 and c = 10, where one path is right)
o = rmfield_safe(opts,{'method'}); o.c0 = c0; o.g0 = 0; o.verified = true;
ws = warning('off','tr_limdist:status');
[F,inf2] = tr_limdist(x,cas,c,o,mode);
warning(ws);
st = 0; if isfield(inf2,'status'), st = inf2.status(1); end
fl = false;
if strcmp(id,'dfx:bound:path'), fl = st ~= 0; return; end
if st == 0, return; end
o.verified = false;
F = tr_limdist(x,cas,c,o,mode);
end

function m = bound_msg(mode,cas,c,c0,x,ME)
% the error of a point that the line integral refuses, with the arguments (path_error of DFExact.jl, dfl_bound_msg of
% the C engine, whose message the MEX file passes on as it is)
if ~isempty(strfind(ME.message,'the line integral needs more nodes than the engine allows')), m = ME.message; return; end
args = sprintf('case %s, c = %g',cas,c);
if c0 ~= 0, args = [args sprintf(', c0 = %g',c0)]; end
args = [args sprintf(', x = %g',x)];
if strcmp(ME.identifier,'dfx:bound:rule'), why = 'at large x above the unit root';
elseif c0 < 0 && x > 0, why = sprintf('c0 < 0 near zero with x > 0: the root %.3g of X(v) = 1/2 + x v + c0 v^2 lies far out',(x + sqrt(x*x + 2*abs(c0)))/(2*abs(c0)));
else, why = 'a root of X(v) far out';
end
if strcmp(cas,'DG'), alt = 'case DG has the line integral representation only';
else, alt = 'the automatic method takes the conditional integral at such points'; end
m = sprintf('ltratio_%s: %s: the line integral needs more nodes than the engine allows here (%s; %s); %s',mode,args,why,ME.message,alt);
end

function o = rmfield_safe(o,names)
for k = 1:numel(names)
    if isfield(o,names{k}), o = rmfield(o,names{k}); end
end
end

function v = getf_(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
