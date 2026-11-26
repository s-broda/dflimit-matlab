function v = dfdet_verify(evalf,CP,vtol)
%DFDET_VERIFY  The checked computation of the general engine (opts.verified = true of the user-level routines).
%
%   v = DFDET_VERIFY(evalf,CP,vtol): evalf(CP',o) evaluates the values on the
%   outer contour CP' with the inner options o (alpha, n, ntail).  Two
%   structurally different evaluations: the default rule, and a variant in
%   which every node moves (the hyperbola 7.5 degrees steeper, at most 88,
%   with 0.8 times the step, or on the vertical line 32 nodes per panel for
%   24; the U-shapes at alpha = 0.3 for 0.5 with 16 and 24 more nodes per
%   panel and tail).  When they agree to vtol (1e-9) the default values are
%   returned; otherwise a third, finer evaluation (half the step, twice the
%   nodes) is returned, with the warning dfdet:verified if it and the variant
%   still differ by more than vtol.  About twice the cost of the default.
%   Where the variant fails with an error (any error but an exhausted
%   memory, which reaches the caller) the variant at the default angle
%   takes its place.
%
%   See also DFDET_CDF, DFDET_CONTOUR.
if nargin < 3 || isempty(vtol), vtol = 1e-9; end
fixed = strcmp(CP.Mo.start,'fixed');
if fixed, n0 = 48; nt0 = 64; else, n0 = 24; nt0 = 40; end
v1 = evalf(CP,struct('alpha',0.5,'n',n0,'ntail',nt0));
try
    v2 = evalf(variant_path(CP,1),struct('alpha',0.3,'n',n0 + 16,'ntail',nt0 + 24));
catch ME
    % (any error, as in the C and Julia engines, except an exhausted memory, which reaches the caller: the groups not
    % separated on the steeper hyperbola, the guard on the nodes of its U-shapes, an error of MATLAB in its node loops)
    if nofallback(ME), rethrow(ME); end
    v2 = evalf(variant_path(CP,0),struct('alpha',0.3,'n',n0 + 16,'ntail',nt0 + 24));
end
d = max([0; abs(v1(:) - v2(:))]);
if d <= vtol || ~all(isfinite(v1(:))), v = v1; return, end
v3 = evalf(variant_path(CP,2),struct('alpha',0.5,'n',2*n0,'ntail',2*nt0));
d3 = max([0; abs(v3(:) - v2(:))]);
if ~(d3 <= vtol)
    warning('dfdet:verified','dfdet, verified: two evaluations by different rules differ by %.2g (above %g) (%s, start %s, c = %g).',d3,vtol,CP.Mo.P.D.name,CP.Mo.start,CP.c);
end
v = v3;
end

function tf = nofallback(ME)
% whether the error ME must reach the caller instead of the fallback: an exhausted memory (MATLAB's out of memory, an array
% beyond the size limit or the preference of the maximum array size), also as the cause of another error
id = ME.identifier;
tf = any(strcmp(id,{'MATLAB:nomem','MATLAB:array:SizeLimitExceeded','MATLAB:pmaxsize','MATLAB:array:sizeLimitExceeded'})) || ...
    ~isempty(regexpi(id,'nomem|outofmemory','once'));
for k = 1:numel(ME.cause), tf = tf || nofallback(ME.cause{k}); end
end

function CPv = variant_path(CP,level)
% the outer contour of a variant rule: level 1 steeper and with a shorter step, level 0 the same with the default
% angle, level 2 the default angle with half the step; on the vertical line more nodes per panel
Mo = CP.Mo; c = CP.c;
if isfield(CP,'TM'), Mo = CP.TM; end      % (the variants of a volatility path: vol_variant of DFDet_vol.jl)
if CP.phi == pi/2 && isnan(CP.h)
    if level == 2, nq = 48; else, nq = 32; end
    CPv = dfdet_cpath(Mo,c,struct('nq',nq,'ts',CP.tstar,'tau',CP.tau));
    return
end
switch level
    case 1, CPv = dfdet_cpath(Mo,c,struct('phi',min(CP.phi + 7.5*pi/180,88*pi/180),'hfac',0.8));
    case 0, CPv = dfdet_cpath(Mo,c,struct('phi',CP.phi,'hfac',0.8));
    otherwise, CPv = dfdet_cpath(Mo,c,struct('phi',CP.phi,'hfac',0.5));
end
end
