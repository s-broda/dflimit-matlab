function [x, w] = glrule(n)
%GLRULE  Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch), cached; shared by DFX_CORE and TR_LIMDIST.
persistent R
if isempty(R), R = {}; end
if n <= numel(R) && ~isempty(R{n}), x = R{n}{1}; w = R{n}{2}; return; end
k = (1:n-1)';
b = k ./ sqrt(4 * k.^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, ix] = sort(diag(D)); V = V(:, ix);
x = (x + 1) / 2; w = V(1, :)'.^2;
R{n} = {x, w};
end
