function D = dfdet_dval(Mo,c,s,t)
%DFDET_DVAL  The determinant of the joint moment generating function of the general engine.
%   D = DFDET_DVAL(Mo,c,s,t) = e^{c + theta} det(I - 2 (s A + t B + gamma E11) Sigma(theta)) for the model Mo
%   (DFDET_MODEL, or a prepared design, DFDET_PREP, for the zero start), so that M(s,t) = e^{-s/2} D^{-1/2}; with a
%   stationary start the extended matrices; with a fixed start that of the zero start (M_{g0} = e^{-s/2} D^{-1/2}
%   e^{g0^2 Qf}).  In cases DU and DT this is D(s,t) of the classical engine, in case DK e^c times it.  A check of the
%   engine against the classical one (DFDET_TEST).
if ~isfield(Mo,'start'), Mo = dfdet_model(Mo); end
P = Mo.P; theta = sqrt(complex(c*c - 2*t)); gam = (c + theta)/2;
if strcmp(Mo.start,'stationary')
    [S,h] = dfdet_sigma(P,theta);
    v = -1/(2*c); n = P.N; N = n + 1;
    Sx = complex(zeros(N)); Sx(1:n,1:n) = S + v*(h*h.'); Sx(N,1:n) = v*h.'; Sx(1:n,N) = v*h; Sx(N,N) = v;
    Q = s*Mo.A + t*Mo.B; Q(1,1) = Q(1,1) + gam; Q(N,N) = Q(N,N) - gam;
    D = exp(c + theta)*det(eye(N) - 2*(Q*Sx));
    return
end
S = dfdet_sigma(P,theta); N = P.N;
M = eye(N) - 2*((s*P.A + t*P.B)*S);
M(1,:) = M(1,:) - 2*gam*S(1,:);
if abs(theta) < 0.5
    z = 2*theta;
    if abs(z) < 1e-3, e = 1 - z/2 + z*z/6 - z^3/24; else, e = (1 - exp(-z))/z; end
    m11 = 1 - (c + theta)*e;
else
    if c > 0, d = -2*t/(theta + c); else, d = theta - c; end
    m11 = (d + (c + theta)*exp(-2*theta))/(2*theta);
end
M(1,1) = M(1,1) + m11 - (1 - 2*gam*S(1,1));
D = exp(c + theta)*det(M);
end
