function e = erfz(z)
%ERFZ  erf for complex z through the Faddeeva function (FADDEEVA): erf(z) = 1 - exp(-z^2) w(iz) for Re z >= 0,
%   erf(z) = -erf(-z) for Re z < 0, and MATLAB's erf on the real axis.  Shared by DFX_CORE and TR_LIMDIST.
e = complex(zeros(size(z)));
r = imag(z) == 0; e(r) = erf(real(z(r)));
pos = ~r & real(z) >= 0; neg = ~r & real(z) < 0;
if any(pos(:)), zp = z(pos); e(pos) = 1 - exp(-zp.^2) .* faddeeva(1i * zp); end
if any(neg(:)), zn = -z(neg); e(neg) = -(1 - exp(-zn.^2) .* faddeeva(1i * zn)); end
end
