function v = ltratio_switch(name,v)
%LTRATIO_SWITCH  The switches of the quantile searches of the t-ratio.
%
%   v = LTRATIO_SWITCH(name): the value of the switch name (true by default):
%     'mom_start'  LTRATIO_ICDF (the classical Newton search, and the general
%                  engine's search where the rule gives it the quantiles)
%                  starts at the normal approximation with the first two
%                  moments from jets (the device of Yong Bao and Raymond Kan,
%                  KAN_CORE), where they hold, in place of the engines'
%                  circles; false: the start of before November 2026;
%     'dist2'      the classical Newton search takes F and f from one pass of
%                  the general engine's paths where the rule (LTRATIO_ENGINE)
%                  gives a Newton step's point to it, the values of LTRATIO_CDF
%                  and LTRATIO_PDF bit for bit (the quantiles are unchanged bit
%                  for bit); false: the two evaluations of before.
%   old = LTRATIO_SWITCH(name,v) sets it and returns the value before.  The
%   switches are for checking (MOM_START and DIST2 of DFLimit.jl, DFL_MOM_START
%   and DFL_DIST2 of the C engine); with the folder of the MEX file of the C
%   engine ahead on the path (LTRATIO_CROUTE('on')) the C engine's switch is set
%   as well.
%
%   See also LTRATIO_ICDF.
persistent S
if isempty(S), S = struct('mom_start',true,'dist2',true); end
name = lower(char(name));
if ~isfield(S,name), error('ltratio:args','ltratio_switch: the switches are ''mom_start'' and ''dist2''.'); end
old = S.(name);
if nargin > 1
    S.(name) = logical(v);
    if ltratio_croute('on'), ltratio_mex('switch',name,double(logical(v))); end
end
v = old;
end
