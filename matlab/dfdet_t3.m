function varargout = dfdet_t3(op,varargin)
%DFDET_T3  The first seasonal t-ratio: its outer contour, the residues, the doubled transform (HEGY_T3_CDF).
%
%   CP = DFDET_T3('path',fn,D,c,opts): the outer contour of t3, the engine's
%   (the present rule) for the design D at c with a zero start, or the
%   volatility path of opts.vol (t3_path of DFDet_seas.jl); refused for c > 0,
%   c < -480 and on a vertical line.  [F,f] = DFDET_T3('eval',CP,x,c0,n,ntail):
%   F and f at the finite points x: the residues at the roots of the engine's
%   determinant (DFDET_CONTOUR 't3'), or where two roots of a group come
%   within 1e-4 (relative) of each other along the contour the engine's
%   inversion on paths of the doubled transform (every root twice, log P
%   twice, the shift a = 1; n, ntail > 0: its inner nodes, which also take it).
%   g = DFDET_T3('gap',CP): the smallest relative distance of two roots of a
%   group along the contour.
%
%   See also HEGY_T3_CDF, HEGY_T3_PDF, HEGY_T3_QUANTILE.
switch op
    case 'path', varargout = {t3_path(varargin{:})};
    case 'eval', [F,f] = t3_eval(varargin{:}); varargout = {F,f};
    case 'gap', varargout = {t3_gap(varargin{:})};
    otherwise, error('dfdet:args','dfdet_t3: unknown operation %s.',op);
end
end

function CP = t3_path(fn,D,c,opts)
if ~(isstruct(D) && isfield(D,'br') && isfield(D,'C') && isfield(D,'detrend')), error('dfdet:args','dfdet: the design must be a structure built by dfdet_determ, dfdet_trend, dfdet_perron, ...'); end
vol = []; if isfield(opts,'vol'), vol = opts.vol; end
if ~isempty(vol)
    [VM,c1] = dfdet_vol('args',fn,D,c,struct('vol',vol),'zero',0,[]);
    if ~isempty(VM), CP = dfdet_cpath(VM,min(VM.cr)); return, end
    c = c1;
end
if ~isscalar(c), error('dfdet:args','%s: c is a vector (one value per regime): that needs a volatility path, opts.vol = volpath(breaks, omega)',fn); end
if ~(isnumeric(c) && isreal(c) && isfinite(c)), error('dfdet:args','%s: c must be a finite real scalar',fn); end
if ~(c <= 0), error('dfdet:args','%s: c must be at most 0 (the seasonal statistics are implemented below and at the unit root)',fn); end
if ~(c >= -480), error('dfdet:args','%s: c must be at least -480 (the range in which the computation is checked)',fn); end
Mo = dfdet_startmodel(D,c,'zero',0,[]);
CP = dfdet_cpath(Mo,c);
if isnan(CP.h), error('dfdet:args','%s: the outer contour of this design is the vertical line (no hyperbola keeps its positive roots right of the origin); the seasonal statistics are not implemented there',fn); end
end

function g = t3_gap(CP)
r = size(CP.R,1); rp = CP.Mo.rp; g = Inf;
for grp = 1:2
    if grp == 1, ix = 1:rp; else, ix = rp+1:r; end
    for i = ix
        for j = ix(ix > i)
            g = min([g abs(CP.R(i,:) - CP.R(j,:))./max(abs(CP.R(i,:)),abs(CP.R(j,:)))]);
        end
    end
end
end

function [F,f] = t3_eval(CP,x,c0,n,ntail)
if nargin < 4, n = []; end
if nargin < 5, ntail = []; end
if isempty(n) && t3_gap(CP) >= 1e-4, [F,f] = dfdet_contour('t3',CP,x,struct('c0',c0)); return, end
% the transform e^{-s}/D(s,t) in the form the engine's inversion on paths integrates (t3_doubled): every root twice
r = size(CP.R,1); rp = CP.Mo.rp;
idx = [reshape([1:rp; 1:rp],1,[]) reshape([rp+1:r; rp+1:r],1,[])];
C2 = CP; C2.R = CP.R(idx,:); C2.lP = 2*CP.lP;
Mo = CP.Mo; Mo.r = 2*r; Mo.rp = 2*rp; Mo.rm = 2*Mo.rm; Mo.a = 1; Mo.kappa = 0.5;
C2.Mo = Mo; C2.fx = [];
if isfield(CP,'key'), C2.key = [CP.key '|t3doubled']; end
o = struct('c0',c0);
if ~isempty(n), o.n = n; end
if ~isempty(ntail), o.ntail = ntail; end
[F,f] = dfdet_contour('tratio',C2,x,o);
end
