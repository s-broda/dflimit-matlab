function m = dfdet_msg(id)
%DFDET_MSG  The texts of the argument errors of the general engine (the messages of DFDet.jl).
switch id
    case 'first_start', m = 'the statistic that detrends first (dfdet_detrend_first) is implemented for the zero start only; start = ''fixed'' (g0 ~= 0) and start = ''stationary'' are not supported with it';
    case 'first_gls', m = 'GLS detrending (opts.gls = cbar) detrends first already: use the design itself, not dfdet_detrend_first(D), with gls';
    case 'first_mom', m = 'the moments of the statistic that detrends first (dfdet_detrend_first) are not implemented (not validated)';
    case 'gls_start', m = 'GLS detrending (opts.gls = cbar) is implemented for the zero start only: the derivation of the GLS limit assumes u_0 = 0, and start = ''fixed'' (g0 ~= 0) and start = ''stationary'' are not supported with it';
    case 'stat_c', m = 'the stationary start needs c < 0';
    case 'vol_stat', m = 'with a volatility path (vol) the stationary start is that of the first regime, J(0) ~ N(0, omega_1^2/(-2 c_1)), which needs c < 0 on the first regime';
    case 'vol_fixed', m = 'with a volatility path (vol) the distribution with a fixed start (g0 ~= 0) is not computed: the paths of the engine''s fixed start are laid for a constant variance and persistence, and under regimes they fail or lose accuracy for starts of moderate size; the transform of the fixed start under regimes is DFDet.vol_qf of DFLimit.jl, and a zero or a stationary start is computed';
    otherwise, m = id;
end
end
