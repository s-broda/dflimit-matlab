function engine_table_matlab(outfile,jfile,sel,autoonly)
%ENGINE_TABLE_MATLAB  The timings behind the choice of engine of the MATLAB programs (LTRATIO_ENGINE).
%
%   ENGINE_TABLE_MATLAB(outfile,jfile) times the two MATLAB engines, the
%   classical programs (opts.engine = 'classical') and the general engine
%   (opts.engine = 'general'), on the cells of the table of the Julia package
%   (dev/engine_table.txt of DFLimit.jl, jfile, whose times it prints beside):
%   the warm time of one call (after a call that builds the outer paths and
%   fills the caches; the smaller of two timed calls when one takes less than
%   2 s, one call otherwise; a call of more than 60 s is timed once, without
%   the warm-up, marked *), the largest difference of the two engines, the
%   ratio of the times, the ratio of the Julia engines, and the engine of the
%   MATLAB rule (LTRATIO_ENGINE) and of the Julia rule in the cell.  A cell
%   outside the checks of the classical programs (LTRATIO_ARGS) is marked
%   n/a.  Writes outfile (default engine_table_matlab.txt beside this file).
%   sel: the numbers of the rows to run (default all).  About two hours (the
%   row of case DU at c = 10 in the bulk is not timed: the general engine takes
%   26 s per value there, the classical programs 1.7 s).  autoonly (false):
%   only the engine of the MATLAB rule per row, without timings.
%
%   See also LTRATIO_ENGINE.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
if nargin < 1 || isempty(outfile), outfile = fullfile(here,'engine_table_matlab.txt'); end
if nargin < 2 || isempty(jfile), jfile = ''; end
if nargin < 3, sel = []; end
if nargin < 4, autoonly = false; end
J = jtable(jfile);
fid = fopen(outfile,'w'); cl = onCleanup(@() fclose(fid));
pr = @(varargin) fprintf_both(fid,varargin{:});
pr('The choice of engine of the MATLAB programs: warm times in ms of the two MATLAB engines (the classical programs,\n');
pr('engine = ''classical'', and the general engine, engine = ''general''), MATLAB %s with its default threads; diff the largest difference of\n',version);
pr('the two; g/c the ratio of the times; Julia g/c the ratio of the engines of DFLimit.jl (dev/engine_table.txt); auto the engine of\n');
pr('the MATLAB rule (LTRATIO_ENGINE) and of the Julia rule (DFLimit.autoengine) in the cell.  * a single call of more than 60 s;\n');
pr('n/a the cell lies outside the checks of the classical programs (LTRATIO_ARGS).\n\n');
R = rows();
ws = warning('off','all'); wcl = onCleanup(@() warning(ws));
cur = '';
for i = 1:numel(R)
    if ~isempty(sel) && ~any(sel == i), continue, end
    r = R{i};
    if ~strcmp(r.sec,cur), cur = r.sec; pr('\n%s\n',cur); end
    if ~isempty(r.xsf)
        % (the points of the cell: mean + k sd from the moments, the classical ones where they are available)
        X = r.xsf(); fc = @() r.fc(X); fg = @() r.fg(X); fa = @() r.auto(X);
    else
        fc = r.fc; fg = r.fg; fa = r.auto;
    end
    if autoonly
        au = fa(); if au, au = 'general'; else, au = 'classical'; end
        pr('  %3d %-62s auto %s\n',i,r.label,au); continue
    end
    [tc,vc,mc] = timed(fc); [tg,vg,mg] = timed(fg);
    if isnan(tc), d = NaN; else, d = maxdiff(vc,vg,r.rel); end
    jr = NaN; ja = ''; if isKey(J,r.label), jj = J(r.label); jr = jj{1}; ja = jj{2}; end
    au = fa(); if au, au = 'general'; else, au = 'classical'; end
    if isnan(tc), cs = '      n/a'; else, cs = sprintf('%9.1f%s',1e3*tc,mc); end
    pr('  %3d %-62s classical %s  general %9.1f%s  g/c %7.2f  diff %8.1e  Julia g/c %6.2f  auto %-9s Julia %s\n', ...
        i,r.label,cs,1e3*tg,mg,tg/tc,d,jr,au,ja);
end
end

function fprintf_both(fid,varargin)
fprintf(fid,varargin{:}); fprintf(varargin{:});
end

function [t,v,mark] = timed(f)
mark = ' ';
try
    t0 = tic; v = f(); t1 = toc(t0);
catch ME
    if strcmp(ME.identifier,'ltratio:args'), t = NaN; v = []; return, end
    rethrow(ME);
end
if t1 > 60, t = t1; mark = '*'; return, end
t0 = tic; v = f(); t = toc(t0);
if t < 2, t0 = tic; f(); t = min(t,toc(t0)); end
end

function d = maxdiff(a,b,rel)
a = a(:); b = b(:);
if rel, d = max(abs(b./a - 1)); else, d = max(abs(a - b)); end
end

function J = jtable(jfile)
% the Julia rows: label -> {ratio g/c, auto}
J = containers.Map();
if isempty(jfile) || ~exist(jfile,'file'), return, end
L = strsplit(fileread(jfile),newline);
for i = 1:numel(L)
    l = L{i};
    k = strfind(l,' classical ');
    if isempty(k), continue, end
    lab = strtrim(l(1:k(1))); tok = regexp(l,'ratio g/c\s+([\d.]+).*auto (\w+)','tokens','once');
    if ~isempty(tok), J(lab) = {str2double(tok{1}),tok{2}}; end
end
end

% =====================================================================================================================
% the cells (those of dev/engine_table.jl of DFLimit.jl)
% =====================================================================================================================
function R = rows()
R = {};
s1 = '(1) the t-ratio, F and f at x = mean + (-2, -1, 0, 1, 2) sd';
for cas = {'DK','DU','DT'}, for c = [0 -10 -50], R{end+1} = tcell(s1,'zero start, c0 = 0, -100 <= c <= 5',cas{1},c); end, end %#ok<AGROW>
R{end+1} = tcell(s1,'zero start, c0 = 0, -100 <= c <= 5','DG',0,0,0,[],[-3.5 -2.8 -2.0 -1.2 -0.5]);
R{end+1} = tcell(s1,'zero start, c0 = 0, -100 <= c <= 5','DG',-10,0,0,[],[-4.5 -3.5 -2.5 -1.5 -0.5]);
for c = [3 5], R{end+1} = tcell(s1,'zero start, c0 = 0, 0 < c <= 5','DT',c); end %#ok<AGROW>
for cas = {'DK','DU','DT'}, R{end+1} = tcell(s1,'zero start, c0 = 0, c < -100',cas{1},-300); end %#ok<AGROW>
for q = {{'DU',-10,-10},{'DK',-5,-5},{'DT',-5,-5},{'DU',-50,-50},{'DK',-20,-20},{'DT',-100,-100}}
    R{end+1} = tcell(s1,'zero start, c0 < 0, -100 <= c <= 0',q{1}{1},q{1}{2},q{1}{3}); %#ok<AGROW>
end
for q = {{'DT',-20,3},{'DU',-5,2},{'DT',-50,5},{'DK',-10,2},{'DK',-2,1},{'DK',0,2},{'DU',0,2},{'DT',1,3}}
    R{end+1} = tcell(s1,'zero start, c0 > 0, -100 < c <= 1',q{1}{1},q{1}{2},q{1}{3}); %#ok<AGROW>
end
for q = {{'DU',-100,3},{'DU',2,1},{'DK',2,2}}, R{end+1} = tcell(s1,'zero start, c0 > 0, otherwise',q{1}{1},q{1}{2},q{1}{3}); end %#ok<AGROW>
for cas = {'DK','DU','DT'}, R{end+1} = tcell(s1,'zero start, c0 != 0, c < -100',cas{1},-300,-300); end %#ok<AGROW>
R{end+1} = tcell(s1,'zero start, c0 != 0, c < -100','DU',-480,-480);
for q = {{'SU',-10},{'ST',-50}}, R{end+1} = tcell(s1,'stationary start, c0 = 0',q{1}{1},q{1}{2}); end %#ok<AGROW>
for q = {{'SK',-5},{'SU',-10},{'ST',-50},{'SU',-100},{'SK',-300}}, R{end+1} = tcell(s1,'stationary start, c0 = c, c <= -5',q{1}{1},q{1}{2},q{1}{2}); end %#ok<AGROW>
for q = {{'SK',-1},{'SU',-0.5},{'SK',-0.1},{'SU',-0.1},{'SU',-3e-3},{'ST',-3e-3}}, R{end+1} = tcell(s1,'stationary start, c0 = c, -5 < c <= -3e-3',q{1}{1},q{1}{2},q{1}{2}); end %#ok<AGROW>
for q = {{'SU',-1e-3},{'ST',-1e-3},{'SK',-1e-4}}, R{end+1} = tcell(s1,'stationary start, c0 = c, c > -3e-3',q{1}{1},q{1}{2},q{1}{2}); end %#ok<AGROW>
R{end+1} = tcell(s1,'stationary start, c0 > 0','ST',-1,2);
for q = {{'DK',5,-3},{'DU',5,-3},{'DT',5,-3},{'DK',4,-0.5}}, R{end+1} = tcell(s1,'above the unit root, c0 < 0, c >= 4',q{1}{1},q{1}{2},q{1}{3},0,[-1 0 1]); end %#ok<AGROW>
for q = {{'DT',7},{'DK',7},{'DU',10}}, R{end+1} = tcell(s1,'above the unit root, c >= 7, the bulk',q{1}{1},q{1}{2},0,0,[-1 0 1]); end %#ok<AGROW>
for q = {{'DK',0,1},{'DK',-0.5,0.5},{'DK',-5,2},{'DK',-50,0.5},{'DK',1,1}}, R{end+1} = tcell(s1,'fixed start, case DK, c <= 1, (c g0)^2/2 <= 400',q{1}{1},q{1}{2},0,q{1}{3}); end %#ok<AGROW>
R{end+1} = tcell(s1,'fixed start, case DK, c <= 1, (c g0)^2/2 <= 400','DK',-10,-2,1);
for q = {{'DK',-100,0.5},{'DK',-80,0.3}}, R{end+1} = tcell(s1,'fixed start, case DK, c >= -100, g0 <= 0.5',q{1}{1},q{1}{2},0,q{1}{3}); end %#ok<AGROW>
for q = {{'DK',2,1},{'DK',-20,2},{'DK',-100,1}}, R{end+1} = tcell(s1,'fixed start, case DK, otherwise',q{1}{1},q{1}{2},0,q{1}{3}); end %#ok<AGROW>
for q = {{'DU',-10,1},{'DT',0,2},{'DU',-50,1.5}}, R{end+1} = tcell(s1,'fixed start, cases DU, DT',q{1}{1},q{1}{2},0,q{1}{3}); end %#ok<AGROW>
s2 = '(2) the moments E t^k, k = 1..4 (diff: relative)';
M = {{'moments','DK',-5},{'moments','DU',-5},{'moments','DT',-5},{'moments','DT',0},{'moments','SK',-5},{'moments','SU',-5},{'moments','ST',-5}, ...
     {'moments, fixed start','DU',-5,0,1},{'moments','DU',-50},{'moments','DU',-100},{'moments','DK',-480}, ...
     {'moments, stationary, c0 != 0','SK',-300,-300},{'moments, stationary, c0 != 0','ST',-200,-200},{'moments, stationary, c0 != 0','SU',-5,-5}, ...
     {'moments, zero start, c0 != 0','DK',-100,-100},{'moments, zero start, c0 != 0','DU',-480,-480},{'moments, zero start, c0 != 0','DT',-300,-300},{'moments, zero start, c0 != 0','DU',-5,-5}};
for q = M, R{end+1} = mcell(s2,q{1}{:}); end %#ok<AGROW>
s3 = '(3) quantiles and the coefficient statistic';
for q = {{'DU',0,0},{'DU',-10,0},{'DT',-5,0},{'ST',-10,0},{'DK',-5,1},{'DU',-5,1},{'DK',0,0},{'DU',-50,0},{'SU',-5,0},{'DT',3,0},{'DG',0,0}}
    cas = q{1}{1}; c = q{1}{2}; g0 = q{1}{3};
    lab = sprintf('ltratio_icdf at p = 0.01, 0.05, 0.5 (%s, c = %s%s)',cas,jnum(c),gs(g0));
    R{end+1} = mk(s3,lab,@() ltratio_icdf([0.01 0.05 0.5],cas,c,struct('g0',g0,'engine','classical')), ...
        @() ltratio_icdf([0.01 0.05 0.5],cas,c,struct('g0',g0,'engine','general')),@() ltratio_engine('icdf',ltratio_spec(cas,struct('g0',g0)),c,0,0)); %#ok<AGROW>
end
for q = {{'DU',0},{'SU',-5},{'DT',-10}}
    cas = q{1}{1}; c = q{1}{2};
    R{end+1} = mk(s3,sprintf('coef_cdf at y = -20, -10, -5 (%s, c = %s)',cas,jnum(c)),@() coef_cdf([-20 -10 -5],cas,c,struct('engine','classical')), ...
        @() coef_cdf([-20 -10 -5],cas,c,struct('engine','general')),@() ltratio_engine('dist',ltratio_spec(cas),c,-10,0)); %#ok<AGROW>
end
R{end+1} = mk(s3,'coef_icdf at p = 0.05 (DU, c = 0)',@() coef_icdf(0.05,'DU',0,struct('engine','classical')),@() coef_icdf(0.05,'DU',0,struct('engine','general')),@() false);
s4 = '(4) the linear form P(U - aV <= b) and the law of the denominator P(V <= b)';
Lq = {{'DK',-3.5,[-0.6 -0.2639 0.3],0},{'DK',-3.5,[-0.6 -0.2639 0.3],-7},{'DK',-1,[-0.4 0 0.5],0},{'DK',2,[-2 0 0.5],-5}, ...
      {'DK',0,[-0.499 -0.4 0 1],0},{'DK',-0.1,[-0.49 -0.4 0 0.5],0},{'DK',0,[-0.4 0.5 3],3},{'DG',-7.3,[-0.9 -0.5 0],0}, ...
      {'DG',0,[-0.49 -0.3 0.5],0},{'DK',NaN,[],-5},{'DK',NaN,[],2},{'DK',NaN,[],5},{'DK',NaN,[],7},{'DK',-1,[-0.3 0 50],5},{'DK',2.5,[-160 -0.49 0],5}};
for q = Lq
    cas = q{1}{1}; a = q{1}{2}; bs = q{1}{3}; c = q{1}{4};
    if isnan(a), bs = evdk(c)*[0.3 1 2.5]; pv = ' (P(V <= b))'; else, pv = ''; end
    S = lspec(cas,-13.5);
    R{end+1} = mk(s4,sprintf('ltratio_lincdf, a = %s%s (%s, c = %s)',jnum(a),pv,cas,jnum(c)),@() ltratio_lincdf(bs,a,cas,c,struct('engine','classical')), ...
        @() ltratio_lincdf(bs,a,cas,c,struct('engine','general')),@() ltratio_engine('lincdf',S,c,0,a)); %#ok<AGROW>
end
for q = {{-5,0,-2},{-2,-5,-7.3},{-1,0,-2},{-1,2,-2}}
    cb = q{1}{1}; c = q{1}{2}; a = q{1}{3}; bs = [-0.9 -0.5 0 0.5];
    R{end+1} = mk(s4,sprintf('ltratio_lincdf, a = %s (DG, cbar = %s, c = %s)',jnum(a),jnum(cb),jnum(c)),@() ltratio_lincdf(bs,a,'DG',c,struct('cbar',cb,'engine','classical')), ...
        @() ltratio_lincdf(bs,a,'DG',c,struct('cbar',cb,'engine','general')),@() ltratio_engine('lincdf',lspec('DG',cb),c,0,a)); %#ok<AGROW>
end
s5 = '(5) the power envelope, ADFmax, regions';
for q = {{'DK',[-5 -10 -20]},{'DK',[-0.01 -0.5]},{'DG',[-5 -10 -20]},{'DG',[-1 -2]}}
    cas = q{1}{1}; cs = q{1}{2};
    lab = sprintf('power_envelope at c = %s (%s)',strjoin(arrayfun(@jnum,cs,'UniformOutput',false),', '),cas);
    if strcmp(cas,'DK'), S = lspec('DK',-13.5); else, S = struct('D',dfdet_trend(1),'start','zero','g0',0,'gls',[],'case','DG'); end
    R{end+1} = mk(s5,lab,@() power_envelope(cs,cas,struct('engine','classical')),@() power_envelope(cs,cas,struct('engine','general')),@() ltratio_engine('envelope',S,cs(1),0,0)); %#ok<AGROW>
end
for q = {{'DU',0},{'DT',-10},{'SU',-10},{'ST',-5},{'SU',-0.1},{'ST',-50}}
    cas = q{1}{1}; c = q{1}{2}; xa = [-3 -2.4 -1.5];
    R{end+1} = mk(s5,sprintf('adfmax_cdf, adfmax_pdf at x = -3, -2.4, -1.5 (%s, c = %s)',cas,jnum(c)), ...
        @() [adfmax_cdf(xa,cas,c,struct('engine','classical')) adfmax_pdf(xa,cas,c,struct('engine','classical'))], ...
        @() [adfmax_cdf(xa,cas,c,struct('engine','general')) adfmax_pdf(xa,cas,c,struct('engine','general'))],@() ltratio_engine('adfmax',ltratio_spec(cas),c,0,0)); %#ok<AGROW>
end
bl = struct('type','line','a',-3.5,'b',-0.2639); bp = struct('type','parabola','x',-1.94087);
R{end+1} = mk(s5,'ltratio_region, line a = -3.5, b = -0.2639 (DK, c = -5)',@() ltratio_region(bl,'DK',-5,struct('engine','classical')),@() ltratio_region(bl,'DK',-5,struct('engine','general')),@() false);
R{end+1} = mk(s5,'ltratio_region, parabola x = -1.94087 (DK, c = -7)',@() ltratio_region(bp,'DK',-7,struct('engine','classical')),@() ltratio_region(bp,'DK',-7,struct('engine','general')),@() false);
% further cells of the MATLAB programs, where their rule differs from that of the Julia package (its boundaries)
s6 = '(6) further cells for the boundaries of the MATLAB rule';
for c = [-20 -30], R{end+1} = tcell(s6,'zero start, c0 = 0, -100 <= c <= 5','DK',c); end %#ok<AGROW>
for q = {{'DU',1},{'DU',2},{'DT',1},{'DK',1}}, R{end+1} = tcell(s6,'zero start, c0 = 0, 0 < c <= 5',q{1}{1},q{1}{2}); end %#ok<AGROW>
for q = {{'SK',-10},{'SK',-50},{'SU',-1},{'ST',-1}}, R{end+1} = tcell(s6,'stationary start, c0 = 0',q{1}{1},q{1}{2}); end %#ok<AGROW>
for q = {{'DU',-20,-20},{'DT',-20,-20},{'DK',-10,-10},{'DU',-1,-1},{'DT',-50,-50}}
    R{end+1} = tcell(s6,'zero start, c0 < 0, -100 <= c <= 0',q{1}{1},q{1}{2},q{1}{3}); %#ok<AGROW>
end
for q = {{'DU',0,1},{'DT',-20,2}}, R{end+1} = tcell(s6,'fixed start, cases DU, DT',q{1}{1},q{1}{2},0,q{1}{3}); end %#ok<AGROW>
end

function r = mk(sec,label,fc,fg,auto,rel,xsf)
if nargin < 6, rel = false; end
if nargin < 7, xsf = []; end
r = struct('sec',sec,'label',label,'fc',fc,'fg',fg,'auto',auto,'rel',rel,'xsf',xsf);
end

function r = tcell(sec,lab,cas,c,c0,g0,ks,xs)
if nargin < 5, c0 = 0; end
if nargin < 6, g0 = 0; end
if nargin < 7 || isempty(ks), ks = [-2 -1 0 1 2]; end
if nargin < 8 || isempty(xs), xsf = @() grid(cas,c,c0,g0,ks); else, xsf = @() xs; end
label = sprintf('%s (%s, c = %g%s%s)',lab,cas,c,cs0(c0),gs(g0));
oc = struct('c0',c0,'g0',g0,'engine','classical'); og = oc; og.engine = 'general';
S = ltratio_spec(cas,struct('g0',g0));
r = mk(sec,label,@(X) [ltratio_cdf(X,cas,c,oc) ltratio_pdf(X,cas,c,oc)],@(X) [ltratio_cdf(X,cas,c,og) ltratio_pdf(X,cas,c,og)], ...
    @(X) ltratio_engine('dist',S,c,c0,X(ceil(numel(X)/2))),false,xsf);
end

function xs = grid(cas,c,c0,g0,ks)
try
    m = ltratio_mom(2,cas,c,struct('c0',c0,'g0',g0,'engine','classical'));
catch ME
    if ~strcmp(ME.identifier,'ltratio:args'), rethrow(ME); end
    m = ltratio_mom(2,cas,c,struct('c0',c0,'g0',g0,'engine','general'));
end
xs = m(1) + ks*sqrt(max(m(2) - m(1)^2,1e-12));
end

function r = mcell(sec,lab,cas,c,c0,g0)
if nargin < 5, c0 = 0; end
if nargin < 6, g0 = 0; end
label = sprintf('%s (%s, c = %s%s%s)',lab,cas,jnum(c),cs0(c0),gs(g0));
oc = struct('c0',c0,'g0',g0,'engine','classical'); og = oc; og.engine = 'general';
r = mk(sec,label,@() ltratio_mom(4,cas,c,oc),@() ltratio_mom(4,cas,c,og),@() ltratio_engine('mom',ltratio_spec(cas,struct('g0',g0)),c,c0,0),true);
end

function S = lspec(cas,cb)
if strcmp(cas,'DG'), S = struct('D',dfdet_trend(1),'start','zero','g0',0,'gls',cb,'case','DG');
else, S = struct('D',dfdet_trend(-1),'start','zero','g0',0,'gls',[],'case','DK'); end
end

function v = evdk(c)
% E V in case DK
if c == 0, v = 0.5; else, v = (exp(2*c) - 1 - 2*c)/(4*c*c); end
end

function s = jnum(v)
% a number as Julia prints a Float64 in the labels of its table
if isnan(v), s = 'NaN'; return, end
if v == round(v) && abs(v) < 1e15, s = sprintf('%.1f',v); return, end
s = sprintf('%.15g',v);
if abs(v) < 1e-4 || abs(v) >= 1e15, s = regexprep(sprintf('%.15g',v),'e-0*(\d)','e-$1'); end
end

function s = cs0(c0)
if c0 == 0, s = ''; else, s = [', c0 = ' jnum(c0)]; end
end

function s = gs(g0)
if g0 == 0, s = ''; else, s = [', g0 = ' jnum(g0)]; end
end
