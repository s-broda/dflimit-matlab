function m = dfdet_msg(id)
%DFDET_MSG  The texts of the argument errors of the general engine (the messages of DFDet.jl).
switch id
    case 'first_start', m = 'the statistic that detrends first (dfdet_detrend_first) is implemented for the zero start only; start = ''fixed'' (g0 ~= 0) and start = ''stationary'' are not supported with it';
    case 'first_gls', m = 'GLS detrending (opts.gls = cbar) detrends first already: use the design itself, not dfdet_detrend_first(D), with gls';
    case 'first_mom', m = 'the moments of the statistic that detrends first (dfdet_detrend_first) are not implemented (not validated)';
    case 'gls_start', m = 'GLS detrending (opts.gls = cbar) is implemented for the zero start only: the derivation of the GLS limit assumes u_0 = 0, and start = ''fixed'' (g0 ~= 0) and start = ''stationary'' are not supported with it';
    case 'stat_c', m = 'the stationary start needs c < 0';
    otherwise, m = id;
end
end
