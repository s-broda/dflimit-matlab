function varargout = dfdet_rules(op,varargin)
%DFDET_RULES  The rules of one evaluation of the general engine's inversion on paths (the t-ratio).
%
%   One evaluation of DFDET_CDF, DFDET_PDF and DFDET_ICDF (and of LTRATIO_CDF
%   and the other user-level routines where the rule takes the general engine)
%   takes a thinner rule than before November 2026 where it holds: the outer
%   hyperbola with 1.5 times the step and the U-shapes with 16 nodes per panel
%   and 32 on the tails of their rays for 24 and 40 (32 and 48 for 48 and 64
%   with a fixed start), and by regions of the start, the deterministic terms,
%   c and c0 the present rule or a finer one; the rules, their regions and
%   their constants are those of DFLimit.jl (PATH_RULES, path_region,
%   path_rule and rule_path of DFDet.jl) and of the C engine
%   (dfl_det_path_region): measured there against a reference with half the
%   outer step and twice the inner nodes, within max(1e-9, twice the present
%   rule's own error) in F and max(1e-8, twice its own) of the largest density
%   in f.  The rules (hfac the factor of the step of the outer hyperbola, n the
%   nodes per panel of the U-shapes, ntail those on the tails of their rays;
%   n, ntail with a fixed start in brackets):
%     'present'  the rule until November 2026: hfac 1, n 24, ntail 40 (48, 64);
%     'thin'     hfac 1.5, n 16, ntail 32 (32, 48);
%     'inner'    the thin inner rule on the present hyperbola;
%     'outer'    the present inner rule on the thin hyperbola;
%     'fine'     the thin inner rule on the present hyperbola with half its step;
%     'finestep' the present inner rule on the present hyperbola with the step
%                min(1/2, 2.5/(g0 |c0|)) times its own where it is flatter than
%                70 degrees (a fixed start without deterministic terms).
%   The regions of opts.rule = 'auto' (the default): near and above the unit
%   root (c > -1) with c0 < 0 'inner' (zero and stationary starts), 'fine' (a
%   fixed start with deterministic terms) or 'finestep' (a fixed start without
%   them; the present step was off by up to 0.86 in the upper tail at c = 0.5,
%   g0 = 1, c0 = -20), with c0 = 0, a fixed start and no deterministic terms
%   'inner'; the present rule beyond the range measured (c > 2 with a fixed
%   start, c > 5, c < -60 with c0 other than 0 and c); 'thin' elsewhere.  On the
%   vertical line (above the unit root) hfac does not apply and the line keeps
%   its 24 nodes per panel.  opts.verified = true takes 'present' (and the
%   variants of DFDET_VERIFY); opts.hfac, opts.n and opts.ntail override the
%   rule.  Under volatility and persistence regimes (opts.vol) the outer
%   hyperbola starts at 70 degrees, where the engine's designs start at 45,
%   and 'thin' and 'outer' take twice the step of that hyperbola (DFDET_VOL
%   'startphi' and 'hfac'); a given opts.phi or opts.hfac is taken as given.
%
%   r = DFDET_RULES('region',Mo,c,c0): the rule of 'auto' for the model Mo
%   (DFDET_STARTMODEL) at c with the local null value c0 (path_region).
%   [hfac,n,ntail,r] = DFDET_RULES('rule',rule,verified,Mo,c,c0,o): the step
%   factor, the inner nodes and the rule taken (path_rule): those of rule
%   ('auto': 'region'; 'present' with verified), with o.hfac (r 'given'),
%   o.n, o.ntail in their place where given.
%   CP = DFDET_RULES('path',Mo,c,phi,hfac,r,c0): the outer contour of the rule
%   r (rule_path: DFDET_CPATH, except that 'fine' and 'finestep' keep the
%   angle of the present rule's hyperbola).
%   T = DFDET_RULES('table'): the rules, a structure of [hfac n ntail nfixed
%   ntailfixed] per name (PATH_RULES).
%
%   See also DFDET_CDF, DFDET_ICDF, DFDET_CPATH.
switch op
    case 'table'
        varargout = {rules_()};
    case 'region'
        varargout = {path_region(varargin{:})};
    case 'rule'
        [hf,n,nt,r] = path_rule(varargin{:});
        varargout = {hf,n,nt,r};
    case 'path'
        varargout = {rule_path(varargin{:})};
    case 'check'
        check_rule(varargin{1});
    otherwise
        error('dfdet:args','dfdet_rules: unknown operation %s.',op);
end
end

function T = rules_()
% PATH_RULES: (hfac, n, ntail, n with a fixed start, ntail with a fixed start)
T = struct('present',[1.0 24 40 48 64],'thin',[1.5 16 32 32 48],'inner',[1.0 16 32 32 48], ...
           'outer',[1.5 24 40 48 64],'fine',[0.5 16 32 32 48],'finestep',[1.0 24 40 48 64]);
end

function r = check_rule(rule)
r = lower(char(rule));
if ~any(strcmp(r,{'auto','present','thin','inner','outer','fine','finestep'}))
    error('dfdet:args','dfdet: opts.rule must be ''auto'', ''present'', ''thin'', ''inner'', ''outer'', ''fine'' or ''finestep'' (the rules of one evaluation of the inversion on paths, DFDET_RULES).');
end
end

function r = path_region(Mo,c,c0)
% the rule of one evaluation by regions of (the class of the design, the start, c, c0), measured in November 2026 (see the
% help).  Near and above the unit root (c > -1) the outer integrand of the upper tail needs the outer step of the present
% rule or a finer one: with c0 < 0, where the parabola X(v) has its root b ~ x/|c0| far out, the thinner outer step lost
% up to 1e-4 with a fixed start (case DK at c = 0, c0 = -5 by 3.9e-5, case DU at c = 2, c0 = -5 by 3.1e-6) and 1.3e-8 of
% the density's peak with a zero start (case DK at c = 0, c0 = -5); without deterministic terms and with a fixed start the
% present step itself was off by up to 0.86 in the upper tail (case DK at c = 0.5, g0 = 1, c0 = -20; 1.8e-3 at c0 = -5),
% where the steps of 'finestep' are within 1e-11 of the classical engine, and with deterministic terms 'fine' within
% 1e-12; there at c0 = 0 too the thinner step lost 1.4e-8 of the peak at mean + 5.75 sd (c = 0.5)
PATH_UNITROOT_CMIN = -1; PATH_THIN_CMAX_FIXED = 2; PATH_THIN_CMAX = 5; FAR_CMIN = -60;
fixed = strcmp(Mo.start,'fixed');
rank1 = Mo.r == 1 && Mo.rm == 0;
if c > PATH_UNITROOT_CMIN
    % (beyond the range of the sweep above the unit root, with a fixed start beyond c = 2 and beyond c = 5 with a zero
    % start: the present rule; at c = 5 with a fixed start the thin rule differed from it by up to 5.6e-6 of the peak of
    % the density, Perron's model C with two breaks)
    if (c > PATH_THIN_CMAX_FIXED && fixed) || c > PATH_THIN_CMAX, r = 'present'; return, end
    if c0 < 0
        if ~fixed, r = 'inner'; elseif rank1, r = 'finestep'; else, r = 'fine'; end
        return
    end
    if c0 == 0 && fixed && rank1, r = 'inner'; return, end
end
% (far from the unit root, below the range of the sweep, with c0 other than 0 and c: the present rule; in case DK at
% c = -480, c0 = 2 the thin rule differed from it by 7.7e-9)
if c < FAR_CMIN && c0 ~= 0 && c0 ~= c, r = 'present'; return, end
r = 'thin';
end

function [hf,n,nt,r] = path_rule(rule,verified,Mo,c,c0,o)
% the step factor of the outer hyperbola, the inner nodes and the rule of one evaluation
r = check_rule(rule);
if verified, r = 'present'; elseif strcmp(r,'auto'), r = path_region(Mo,c,c0); end
T = rules_(); v = T.(r);
hf = v(1);
if strcmp(Mo.start,'fixed'), n = v(4); nt = v(5); else, n = v(2); nt = v(3); end
if isfield(o,'hfac') && ~isempty(o.hfac), hf = double(o.hfac); r = 'given'; end
if isfield(o,'n') && ~isempty(o.n), n = o.n; end
if isfield(o,'ntail') && ~isempty(o.ntail), nt = o.ntail; end
end

function CP = rule_path(Mo,c,phi,hf,r,c0)
% the outer contour of the rule r with the step factor hf: DFDET_CPATH, except that the rules with a finer step keep the
% angle of the present rule's hyperbola (where the angle is chosen from the largest term, its criterion weighs the step,
% and a finer step took a flatter hyperbola with larger terms: case DK at c = -0.5, g0 = 2, c0 = -20 was off by 1.7e-4 on
% it, against 4.5e-12 on the present one): 'fine' half the present step, 'finestep' (a fixed start without deterministic
% terms) the step min(1/2, 2.5/(g0 |c0|)) times the present one where the present hyperbola is flatter than 70 degrees (at
% c = 0.5, g0 = 1, c0 = -20 the present step was off by 0.86, half of it by 0.31, this one by 7e-12; on the steep
% hyperbolas of g0 = 2, 88 degrees, the present step is within 1e-12 and the thin inner rule is not); on the vertical line
% the line itself
FINESTEP_PHI = 70*pi/180;
if isempty(phi) && any(strcmp(r,{'fine','finestep'}))
    CP1 = dfdet_cpath(Mo,c,struct('phi',[],'hfac',1));
    if isnan(CP1.h), CP = CP1; return, end
    phi = CP1.phi;
    if strcmp(r,'finestep')
        if CP1.phi < FINESTEP_PHI, hf = min(0.5,2.5/(abs(Mo.g0)*abs(c0))); else, hf = 1; end
    end
end
CP = dfdet_cpath(Mo,c,struct('phi',phi,'hfac',hf));
end
