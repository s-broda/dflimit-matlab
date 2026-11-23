function varargout = ltratio_rules(what,varargin)
%LTRATIO_RULES  The rules of one evaluation of the classical programs' contour and conditional integrals.
%
%   [n,hfac] = LTRATIO_RULES('contour',c,c0,g0): the rule of one evaluation
%   of the contour integral (TR_CONTOUR with opts.verified false; ctr_rule of
%   DFLimit.jl, dfl_ctr_rule of the C engine): n Gauss-Legendre nodes per panel
%   and per ray of the U-shapes and the factor hfac of the outer step, by
%   regions of c, c0 and c g0 measured in November 2026 against the checked
%   computation on a dense grid of x (cases DU, DT, g0 from 0.5 to 3, c from 2
%   to -60, c0 in {0, 2, -2, -5, c}; within 6.2e-11 in F).  Before, n = 64 and
%   hfac = 1 everywhere (opts.rule = 'present' of TR_CONTOUR).  The outer step
%   1.5 times as long holds for |c g0| <= 20 and c <= 1 (above the unit root it
%   lost 1.5e-2 in case DU at c = 2, g0 = 1.5; at c g0 = -40.5 in case DT
%   3.9e-8, at -90 2.4e-2); 32 nodes hold near the unit root (with c0 > 0 down
%   to c = -2.5: in case DU at c = -5, g0 = 1.5, c0 = 2 they lost 2.1e-9 in F)
%   and down to c = -10 (from c = -13.5 on they lost up to 3e-4 of the peak of
%   the density).  0.42 times the time near the unit root (median), 0.59 above
%   it.
%
%   n = LTRATIO_RULES('conditional',model,c,c0,g0,stat): the rule of one
%   evaluation of the conditional integral (TR_LIMDIST with opts.verified
%   false; cond_rule of DFLimit.jl, dfl_cond_rule of the C engine): the nodes
%   Nv of the path of each inner inversion and Nl of the integral over l, by
%   regions measured in November 2026 against 192 nodes and the outer
%   tolerance 1e-14 on a dense grid of x.  Before, 96 nodes everywhere.  48
%   nodes hold with a stationary start for c <= -2 and with a zero start in
%   cases DU, DT (models 2, 3) for c >= -5 and in case DK for c <= -13.5; 64 with a
%   stationary start for c > -2 (48 lost 7e-8 of the peak of the density in
%   case SK at c = -0.5); 96 elsewhere (with a zero start in cases DU, DT 48
%   nodes lost 1.8e-8 of the peak at c = -13.5 and up to 1e-3 far from the
%   unit root), with a fixed start (measured at one configuration only: 48
%   nodes lost 7e-7 in case DU at c = 2, c0 = -2) and above c = 2.  The options
%   Nv, Nl of TR_LIMDIST override it.
%
%   See also TR_CONTOUR, TR_LIMDIST, DFDET_RULES.
switch what
    case 'contour'
        [c,c0,g0] = varargin{:};
        CTR_RULE_CG0 = 20;
        if c > 1, n = 32; h = 1;
        elseif abs(c*g0) > CTR_RULE_CG0, n = 64; h = 1;
        else
            if c >= -2.5 || (c >= -10 && c0 <= 0), n = 32; else, n = 64; end
            h = 1.5;
        end
        varargout = {n,h};
    case 'conditional'
        [model,c,~,g0,stat] = varargin{:};
        COND_RULE_CMAX = 2;
        if stat
            if c <= -2, n = 48; else, n = 64; end
        elseif g0 ~= 0
            n = 96;
        elseif model == 1
            if c <= -13.5, n = 48; else, n = 96; end
        else
            if c >= -5, n = 48; else, n = 96; end
        end
        if c > COND_RULE_CMAX, n = 96; end                 % (beyond the range measured)
        varargout = {n};
    otherwise
        error('ltratio:args','ltratio_rules: what must be ''contour'' or ''conditional''.');
end
end
