function [x, w] = glrule(n)
%GLRULE  Gauss-Legendre nodes and weights on [0, 1], cached; shared by DFX_CORE and TR_LIMDIST.
%   Golub-Welsch up to 1024 nodes; above, where its matrix of eigenvectors (8 n^2 bytes, 134 MB at n = 4096) would
%   dominate the memory and its time grows like n^3, Newton's method on the three-term recurrence of the Legendre
%   polynomial from the Chebyshev initial guesses, in memory proportional to n (the rule of the C engine, dfl_gl_rule;
%   gl_newton of DFExact.jl; against Golub-Welsch the nodes agree to 5e-16 and the weights to 4e-10 relative at
%   n = 4096, the small weights at the ends, which the eigenvectors give less accurately).  At most 4096 nodes, the bound
%   of the C and Julia engines (DFL_RULEMAX, RULEMAX): a larger rule raises the error dfx:bound:rule, which the user-level
%   routines turn into the conditional integral or into an error that names the arguments (LTRATIO_DIST).
persistent R
if isempty(R), R = {}; end
if n > 4096
    error('dfx:bound:rule', 'a Gauss-Legendre rule of the inner integrals needs %.3g nodes, at most 4096', n);
end
if n <= numel(R) && ~isempty(R{n}), x = R{n}{1}; w = R{n}{2}; return; end
if n <= 1024
    k = (1:n-1)';
    b = k ./ sqrt(4 * k.^2 - 1);
    [V, D] = eig(diag(b, 1) + diag(b, -1));
    [x, ix] = sort(diag(D)); V = V(:, ix);
    x = (x + 1) / 2; w = V(1, :)'.^2;
else
    t = cos(pi * ((0:n-1)' + 0.75) / (n + 0.5));
    act = (1:n)';                                              % the nodes still iterated (each stops at |dt| < 1e-16)
    for it = 1:100
        [p1, p0] = legendre_pn(n, t(act));
        dp = n * (t(act) .* p1 - p0) ./ (t(act) .* t(act) - 1);
        dt = p1 ./ dp; t(act) = t(act) - dt;
        act = act(~(abs(dt) < 1e-16));
        if isempty(act), break; end
    end
    [p1, p0] = legendre_pn(n, t);
    dp = n * (t .* p1 - p0) ./ (t .* t - 1);
    x = flipud((1 + t) / 2);                                   % ascending, as Golub-Welsch
    w = flipud(1 ./ ((1 - t .* t) .* dp .* dp));             % 2/((1 - t^2) P_n'(t)^2) on [-1, 1], halved for [0, 1]
end
R{n} = {x, w};
end

function [p1, p0] = legendre_pn(n, t)
% P_n(t) and P_{n-1}(t) by the three-term recurrence
p0 = ones(size(t)); p1 = t;
for k = 2:n
    p2 = ((2 * k - 1) * t .* p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2;
end
end
