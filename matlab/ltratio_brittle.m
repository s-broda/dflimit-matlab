function msg = ltratio_brittle(cas,c,c0,g0,used)
%LTRATIO_BRITTLE  The warning of the one-evaluation default in a regime known to be brittle.
%
%   msg = LTRATIO_BRITTLE(cas,c,c0,g0,used) gives the warning ltratio:brittle
%   (and returns its text, '' if none) when the representations used (a cell
%   of 'single', 'double', 'contour') with opts.verified = false, the default,
%   lie in a regime where one evaluation is known to be brittle, decided from
%   the parameters alone (no computed diagnostics): the boundaries lie where
%   one evaluation stops being good to about 1e-5 absolute, with a margin
%   (November 2026; brittle_msg of julia/DFExact.jl, dfl_brittle of the C
%   engine):
%     contour, c0 > 0, c < -15      more than 1e-5 (6e-6 at c = -27, c0 = 2.3;
%                                   in the density 2e-2 at c = -60, c0 = 5)
%     contour, c < -20, |c g0| > 60 more than 1e-5 (3.9e-5 at c = -60, g0 = 1.5,
%                                   c0 = -2; the router takes the general engine
%                                   there, LTRATIO_ENGINE)
%     contour, c0 < 0, c > 1        the far upper tail fails
%     contour, case DU, c0 > 3, c > 1.5   fails (outright at c = 2, c0 >= 5)
%     double, zero start, c0 > 0, c < 0    up to 1.5e-2 at x < 0
%     double, case DK, c0 ~= 0, c > -5     up to 1.2e-2 in the upper tail
%     double, zero start, c < -100; stationary start, c < -50
%     double, case SK, c > -0.5            up to 3e-3 near c = 0
%     double, fixed start, case DK, c < -20; cases DU, DT, c < -30
%                                          (3e-5 at c = -27.5 in case DK,
%                                          1.2e-4 at c = -50 in case DU)
%   opts.verified = true checks each value.
%   LTRATIO_ICDF and COEF_ICDF give it once per call.
%
%   LTRATIO_BRITTLE('warn',id,msg) gives the warning id with the text msg under the same rule (once per call of
%   LTRATIO_ICDF, COEF_CDF, COEF_ICDF): the warning ltratio:bound of LTRATIO_DIST, where the line integral refused points
%   whose values from the conditional integral are flagged.
%
%   See also LTRATIO_CDF, LTRATIO_METHOD.
persistent depth seen seenid
if isempty(depth), depth = 0; seen = {}; seenid = {}; end
if nargin == 1
    % ltratio_brittle('begin') ... ltratio_brittle('end'): the messages in between are given once each at the end
    % (LTRATIO_ICDF, COEF_CDF, COEF_ICDF evaluate the distribution many times)
    switch cas
        case 'begin', depth = depth + 1;
        case 'end'
            depth = max(depth - 1, 0);
            if depth == 0
                [u, iu] = unique(seen,'stable'); ids = seenid(iu); seen = {}; seenid = {};
                for k = 1:numel(u), warning(ids{k},'%s',u{k}); end
            end
    end
    msg = ''; return
end
if nargin == 3 && strcmp(cas,'warn')
    msg = c0;
    if depth > 0, seen{end+1} = msg; seenid{end+1} = c; else, warning(c,'%s',msg); end
    return
end
msg = '';
hint = '; opts.verified = true checks each value';
stat = cas(1) == 'S';
% (the contour integral with c0 > 0 below c = -15: 6e-6 at c = -27, c0 = 2.3, (c g0)^2/2 = 146, in the density up to 2e-2
% at c = -60, g0 = 1.5, c0 = 5 in case DU; and with any c0 below c = -20 where |c g0| > 60: in case DU 3.9e-5 at c = -60,
% g0 = 1.5, c0 = -2, and 6.3e-6 at c = -100, g0 = 0.9; brittle_msg of DFExact.jl, November 2026)
if any(strcmp(used,'contour')) && c0 > 0 && c < -15
    msg = sprintf('%s, fixed start, c0 > 0 far from the unit root (c < -15): the contour integral by one evaluation may be off by more than 1e-5%s',cas,hint);
elseif any(strcmp(used,'contour')) && c < -20 && abs(c*g0) > 60
    msg = sprintf('%s, fixed start far from the unit root with |c g0| > 60: the contour integral by one evaluation may be off by more than 1e-5%s',cas,hint);
elseif any(strcmp(used,'contour')) && c > 1 && c0 < 0
    msg = sprintf('%s, fixed start, c0 < 0 at c > 1: the contour integral by one evaluation fails in the far upper tail%s',cas,hint);
elseif any(strcmp(used,'contour')) && strcmp(cas,'DU') && c > 1.5 && c0 > 3
    msg = sprintf('%s, fixed start, c0 > 3 at c > 1.5: the contour integral by one evaluation fails (the automatic method takes the conditional integral there)%s',cas,hint);
elseif any(strcmp(used,'double'))
    reg = '';
    if g0 ~= 0
        if strcmp(cas,'DK') && c < -20, reg = 'fixed start far from the unit root (c < -20)';
        elseif ~strcmp(cas,'DK') && c < -30, reg = 'fixed start far from the unit root (c < -30)';
        end
    elseif ~stat && c0 > 0 && c < 0, reg = 'c0 > 0 below the unit root';
    elseif ~stat && strcmp(cas,'DK') && c0 ~= 0 && c > -5, reg = 'case DK with c0 != 0 near or above the unit root (c > -5)';
    elseif ~stat && c < -100, reg = 'c < -100';
    elseif stat && c < -50, reg = 'a stationary start with c < -50';
    elseif strcmp(cas,'SK') && c > -0.5, reg = 'case SK near the unit root (c > -0.5)';
    end
    if ~isempty(reg), msg = sprintf('%s, %s: the conditional integral by one evaluation may be off by more than 1e-5%s',cas,reg,hint); end
end
if ~isempty(msg)
    if depth > 0, seen{end+1} = msg; seenid{end+1} = 'ltratio:brittle'; else, warning('ltratio:brittle','%s',msg); end
end
end
