function [F,info] = adfmax_dist(x,cas,c,opts,mode)
%ADFMAX_DIST  Argument checks and non-finite x of ADFMAX_CDF and ADFMAX_PDF (mode 'cdf' or 'pdf').
%
%   opts.engine: 'classical' TR_ADFMAX; 'general' the general engine
%   (DFDET_ADFMAX_CDF, DFDET_ADFMAX_PDF); 'auto' (the default) the engine of
%   LTRATIO_ENGINE, whose values below 1e-3 go back to the conditional integral
%   (relative errors).  cas may be a design of the general engine (LTRATIO_SPEC):
%   cases DU, DT, SU, ST by the classical programs, any other design by the
%   general engine.
%
%   See also ADFMAX_CDF, TR_ADFMAX, DFDET_ADFMAX_CDF.
if nargin < 4 || isempty(opts), opts = struct; end
[engine,oc] = ltratio_route('engine',opts);
if ~(ischar(cas) || isstring(cas))
    S = ltratio_spec(cas,opts);
    if ~strcmp(engine,'general') && any(strcmp(S.case,{'DU','DT','SU','ST'}))
        oc.engine = engine; [F,info] = adfmax_dist(x,S.case,c,oc,mode); return
    end
    if strcmp(engine,'classical'), error('ltratio:args','adfmax: engine = ''classical'': ADFmax of the classical programs has the cases DU, DT, SU, ST only.'); end
    F = general_adfmax(x,S,c,oc,mode); info = struct('engine','general'); return
end
if strcmp(engine,'general')
    F = general_adfmax(x,ltratio_spec(cas,oc),c,oc,mode); info = struct('engine','general'); return
end
opts = oc;
if ~(ischar(cas) || isstring(cas)) || ~any(strcmpi(cas,{'DU','DT','SU','ST'}))
    error('ltratio:args','adfmax: the case must be one of DU, DT, SU, ST.');
end
cas = upper(char(cas));
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c)), error('ltratio:args','adfmax: c must be a finite real scalar.'); end
if cas(1) == 'S' && ~(c <= -1e-4), error('ltratio:args','ltratio: a stationary start (cases SK, SU, ST) requires c <= -1e-4: closer to zero the variance -1/(2c) of the start is too large for the quadrature.'); end   % (LTRATIO_ARGS)
if c < -550, error('ltratio:args','adfmax: c must be at least -550 (the range in which the programs are checked).'); end
if ~(isnumeric(x) && isreal(x)), error('ltratio:args','adfmax: x must be real.'); end
F = zeros(size(x)); F(isnan(x)) = NaN; if strcmp(mode,'cdf'), F(x == Inf) = 1; end
fin = isfinite(x); info = struct('err',NaN(size(x)),'status',zeros(size(x)));
if strcmp(engine,'auto') && isempty(fieldnames(opts)) && any(fin(:)) && ltratio_engine('adfmax',ltratio_spec(cas,opts),c,0,0)
    % the general engine; its values below 1e-3 by the conditional integral, whose error is relative
    Fg = general_adfmax(x(fin),ltratio_spec(cas,opts),c,opts,mode);
    F(fin) = Fg;
    ig = find(fin); tl = ig(Fg < 1e-3);
    if ~isempty(tl), [Fv,inf1] = tr_adfmax(x(tl),cas,c,opts,mode); F(tl) = Fv; info.err(tl) = inf1.err; info.status(tl) = inf1.status; end
    if strcmp(mode,'pdf'), F(F < 0) = 0; end
    info.engine = 'general';
    return
end
if any(fin(:))
    [Fv,inf1] = tr_adfmax(x(fin),cas,c,opts,mode);
    F(fin) = Fv; info.err(fin) = inf1.err; info.status(fin) = inf1.status;
end
if strcmp(mode,'pdf'), F(F < 0) = 0; end
end

function F = general_adfmax(x,S,c,oc,mode)
o = ltratio_route('gopts',S,oc); o.gls = [];
if strcmp(mode,'cdf'), F = dfdet_adfmax_cdf(x,S.D,c,o); else, F = dfdet_adfmax_pdf(x,S.D,c,o); F(F < 0) = 0; end
end
