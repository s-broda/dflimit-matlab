# The C engine in MATLAB

`ltratio_mex` is a MEX file built from the C engine in `c/` (a port of the Julia engine `DFExact.jl` of the package
DFLimit.jl, checked against its values by the tests in `c/test/`).  It needs nothing else: no runtime, no library, 2.4 MB.  Two ways to use it:

- **Self-contained**: `ltratio.m` offers the functions of the MATLAB programs in `matlab/` with the same arguments.
- **Under the MATLAB programs**: `dfx_core.m`, `tr_limdist.m`, `tr_contour.m` and `tr_limmom.m` have the names and
  signatures of the engines in `matlab/`; with this folder ahead of `matlab/` on the path, `ltratio_cdf`, `ltratio_pdf`, `ltratio_icdf`,
  `ltratio_mom`, `ltratio_test` and the paper's replication script `ltratio_paper` (in the paper's replication package) run unchanged on the C engine (`ltratio_test` passes all 268
  checks in about 30 s instead of about five minutes, when started from a folder other than `matlab/`: the current folder
  comes first on MATLAB's path).  In the repository this folder also has `tr_adfmax.m`, `tr_region.m` and
  `tr_gilpelaez.m`, the replacements of the engines of the applications in `matlab/`: `adfmax_cdf`, `adfmax_pdf`,
  `adfmax_icdf` (the ADFmax test), `ltratio_region` (the regions {U <= g(V)} of case DK, for the boundary types
  `'line'`, `'parabola'` and `'wap'`; a function handle `g` needs the MATLAB engine), `ltratio_lincdf` (the linear
  form U - aV) and `power_envelope` (the envelope of Elliott, Rothenberg and Stock) then run on the C engine
  (`dfl_adfmax`, `dfl_region` in `c/dfl_double.c`, `dfl_lincdf`, `dfl_power_envelope` in `c/dfl_gilpelaez.c`; the
  replacement `tr_gilpelaez.m` computes the characteristic function anew at every call of the bisection of
  `power_envelope`, which `ltratio('envelope', ...)` does not; near c = 0, where the envelope's critical value comes
  from the conditional integral, `power_envelope` calls `ltratio_mex('envelope', ...)` itself).

```matlab
addpath ltratio-macos                                 % the unpacked zip (see below)
ltratio('icdf', [0.01 0.05 0.10], 'DU')               % asymptotic critical values of tau_mu
ltratio('cdf', -2.86157, 'DU', -5)                    % limiting power of the 5% test at c = -5
ltratio('cdf', -2, 'DU', -5, struct('g0', 1.5))       % fixed start
ltratio('pdf', [-1 0 1], 'DK', 2, struct('c0', 2))    % t-ratio centered at the true value
ltratio('mom', 4, 'SU', -5)                           % moments
ltratio('lincdf', -0.2639, -3.5, 'DK', 0)             % P(U - aV <= b): the size of the point-optimal test at cbar = -7
ltratio('envelope', [-5 -10], 'DG')                   % the power envelope of the trend case
[F, e] = ltratio('coefcdf', -14.0936, 'DU', -10)      % the coefficient statistic (the single inversion of Bao and Kan) and its error estimate
ltratio('coefpdf', -8.039, 'DK')                      % its density
ltratio('coefmom', 2, 'DK', -5)                       % its moments
ltratio('threads')                                    % the number of threads (all processors; set with a second argument)
```

The general engine (`c/dfl_det_*.c`, any piecewise smooth deterministic design) is reached through designs:
`D = ltratio('design', kind, ...)` makes one (`'trend'` p, `'perron'` model lam with several breaks, `'ltrend'` p,
`'matrix'` X [TB] from the regressor matrix of a sample, `'first'` D for the statistic that detrends first, `'fun'` f K
[breaks] from a function handle), a structure with the break points `br` and the coefficient matrices `C` of the pieces
side by side; `ltratio('cdf', x, D, c, opts)` (and `'pdf'`, `'icdf'`, `'mom'`) then takes the router over the two
engines (`opts.engine` `'auto'`, the default, `'classical'` or `'general'`; `opts.start` `'zero'`, `'fixed'` with
`opts.g0`, `'stationary'`; `opts.gls` = cbar for GLS detrending), and `ltratio('det', fn, D, args..., opts)` the general
engine's other functions: the coefficient statistic, the linear form and the law of V, ADFmax, local power, the
envelope, the point-optimal statistic, the LM statistic, MZ_alpha and MSB (`ltratio_mex('det', ...)`, the list in the
header of `ltratio_mex.c`).  A case string takes the router too, with `opts.engine` `'auto'` by default (also for
`'lincdf'` and `'envelope'`, without `opts.N` and `opts.umax`); `opts.engine = 'classical'`, a method named or
`verified` compute with the classical engine as before, and where the rule takes the general engine for a case string
the two engines agree to about 1e-9.  `det_check.m` checks a built MEX file on the general engine.  With this folder ahead of
`matlab/` on the path the user-level routines of the MATLAB programs (`ltratio_cdf`, `ltratio_pdf`, `ltratio_icdf`,
`ltratio_mom`, `coef_cdf`, `coef_icdf`, `adfmax_cdf`, `adfmax_pdf`, `ltratio_lincdf`, `power_envelope`,
`ltratio_region`) hand a call with a design or a case string, without options of the MATLAB engines, to the C router
(`ltratio_mex('route', ...)`, `matlab/ltratio_croute.m`; `opts.engine = 'auto'` is the default, for a case string as
for a design, and `'classical'` keeps the classical programs, here on the C engine through the drop-ins): the C rule
chooses the engine, both engines compute in C, and the warning `ltratio:brittle`
names a brittle regime of the classical engine's one evaluation where the router takes it.  The route operation also
has ADFmax, the linear form, the envelope and the regions (`'adfcdf'`, `'adfpdf'`, `'lincdf'`, `'envelope'`,
`'region'`).

```matlab
D = ltratio('design', 'perron', 'C', [0.3 0.7]);      % Perron's model C with two breaks
ltratio('icdf', [0.01 0.05 0.10], D)                  % its critical values
ltratio('cdf', -4, D, -10, struct('start', 'stationary'))
ltratio('det', 'lincdf', ltratio('design', 'trend', 1), -0.5, -7.3, 0, struct('gls', -13.5))
ltratio('lincdf', -0.2639, -3.5, ltratio('design', 'trend', 0), -5, struct('start', 'stationary'))   % the linear form through a design
```

The case `'DG'` is the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996) at cbar = -13.5 (the DF-GLS
tests of the paper): zero start, the single integral only, no moments, and no tail recomputation (far in the tails the values carry the single integral's absolute error).  `opts.method` is `'auto'`
(the default), `'single'` or `'line'` (the line integral representation), `'double'` or `'conditional'` (the conditional integral
representation) or `'contour'` (the contour integral representation: the inversion formula on paths around the
singularities of the moment generating function, with the integral over v in closed form and the paths in s and t deformed, a double integral with an explicit integrand, for a
fixed start g0 not 0 in cases DU, DT only; `dfl_contour` in `c/dfl_gen.c`, a port of `tr_contour` of `DFExact.jl`).
By default (`opts.verified = false`) every value is one evaluation by one method: no second evaluation, no refinement
driven by an error estimate, no fallback to another representation and no recomputation of the tails. `'auto'` takes the
fastest representation that one evaluation gives to well below 1e-5 (`dfl_automethod`): with a fixed start the contour
integral in cases DU, DT, except where it fails (c > 1 with c0 < 0; in case DU c > 1.5 with c0 > 3), where it takes the
conditional integral, and the line integral (the noncentral kernel) in case DK; with a zero or stationary start the
conditional integral with c0 < 0 and a zero start (cases DU, DT; case DK for c <= -5) and with a stationary start in case
SK for -50 <= c <= -0.5 and with c0 > 0 at x < 0 for c >= 0 (cases DU, DT), and the line integral elsewhere (c0 = 0, a
zero start with c0 > 0 otherwise, a stationary start in cases
SU, ST, and far from the unit root: c < -100 with a zero start, c < -50 with a stationary one). Where the parameters lie
in a regime in which one evaluation is known to be brittle the warning `ltratio:brittle` names it, once per call
(`ltratio_mex('brittle', ...)`, `dfl_brittle`). `opts.verified = true` gives the checked computation at 1 to 4.5 times
the cost of the default (about 1.8 with the contour integral; more in the tails below 1e-3, which it recomputes): error estimates and a status, the agreement of two evaluations on different paths, the fallbacks to another
representation and the recomputation of the values below 1e-3 by the conditional integral, with the rule before
November 2026. A computed density below zero is returned as zero.  The MEX file itself takes numeric arguments only (see the header of
`ltratio_mex.c`): `'single'`, `'double'`, `'contour'`, `'dist'`, `'icdf'`, `'mom'`, `'tstar'`, `'adfmax'`, `'region'`, `'lincdf'`, `'envelope'`, `'threads'`, `'brittle'`
(`'dist'` and `'icdf'` take the method 0 automatic, 1 single, 2 double, 3 contour, and an optional last argument verified; `[F, err, status, nodes, level] =
ltratio_mex('contour', kind, x, model, c, c0, g0)` returns the error estimate, the status, the number of outer nodes and
the level of the inner rules at each point).

Seconds per value of the distribution function by the default (one evaluation), MATLAB R2026b on an M4 laptop (ten
cores), `timeit` medians of warm calls, November 2026 (the script `timing_fixed.m` produces the rows of the second table):

| case c c0 x | MATLAB programs (1 thread) single / double | C MEX 1 thread single / double | C MEX 10 threads single / double |
|---|---|---|---|
| DK 0 0 -1.94 | 0.034 / 0.019 | 0.0090 / 0.0018 | 0.0066 / 0.0022 |
| DU 0 0 -2.86 | 0.065 / 0.082 | 0.012 / 0.020 | 0.0027 / 0.0095 |
| DT 0 0 -3.41 | 0.059 / 0.109 | 0.011 / 0.025 | 0.0025 / 0.011 |
| DU -5 0 -2.86 | 0.067 / 0.051 | 0.013 / 0.014 | 0.0028 / 0.0087 |
| DU 10 0 20 | 0.106 / 0.107 | 0.027 / 0.039 | 0.0060 / 0.012 |
| DU 0 1 -1 | 0.191 / 0.059 | 0.040 / 0.021 | 0.0083 / 0.0093 |
| DT -2.5 -5 -1.1 | 0.182 / 0.099 | 0.037 / 0.017 | 0.0079 / 0.0077 |
| DK 2 2 0 | 0.0056 / 0.044 | 0.0005 / 0.0052 | 0.0002 / 0.0014 |
| DU 5 -5 126.85 | 1.12 / 0.102 | 0.394 / 0.026 | 0.074 / 0.0091 |
| SU -5 0 -2 | 0.084 / 0.410 | 0.016 / 0.039 | 0.0036 / 0.016 |
| ST -10 0 -3 | 0.066 / 0.513 | 0.013 / 0.041 | 0.0029 / 0.017 |

With a fixed start the MATLAB programs have the conditional integral and, in cases DU, DT, the contour integral (on one
thread, as in the first table); the line integral uses the noncentral kernel, written in C to be vectorized, and
`'auto'` takes it in case DK and the contour integral in cases DU, DT:

| case c c0 x, g0 | MATLAB programs (1 thread) double / contour | C MEX 1 thread single / double / contour | C MEX 10 threads single / double / contour |
|---|---|---|---|
| DK -5 0 -2, g0 = 1 | 0.066 / -- | 0.042 / 0.0099 / -- | 0.0082 / 0.0044 / -- |
| DK 2 0 1, g0 = 1 | 0.101 / -- | 0.0039 / 0.022 / -- | 0.0009 / 0.010 / -- |
| DU -5 0 -2.86, g0 = 1 | 1.52 / 0.137 | 0.327 / 0.741 / 0.043 | 0.078 / 0.247 / 0.011 |
| DU -5 0 0.5, g0 = 1 | 1.34 / 0.047 | 0.032 / 0.674 / 0.019 | 0.0080 / 0.251 / 0.0055 |
| DU 2 0 -3, g0 = 1 | 1.35 / 0.136 | 0.334 / 0.876 / 0.051 | 0.082 / 0.275 / 0.012 |
| DT -5 0 -3, g0 = 2 | 1.67 / 0.109 | 0.302 / 0.774 / 0.045 | 0.073 / 0.258 / 0.011 |
| DT 2 0 -2, g0 = 1 | 2.13 / 0.109 | 0.366 / 0.814 / 0.047 | 0.089 / 0.264 / 0.011 |
| DU -5 -5 -2, g0 = 1 | 1.52 / 0.086 | 0.321 / 0.721 / 0.037 | 0.077 / 0.240 / 0.0095 |
| DT -2.5 -5 -1.1, g0 = 1 | 1.59 / 0.073 | 0.351 / 0.757 / 0.031 | 0.086 / 0.263 / 0.0079 |

## Downloading a built MEX file

The workflow "MATLAB MEX" (`.github/workflows/mex.yml`) builds, when started by hand (Actions tab, "Run workflow") and on every tag `v*`, one zip
per platform (`ltratio-linux-x86_64`, `ltratio-macos` with the Apple silicon and the Intel MEX files,
`ltratio-windows-x86_64`) and attaches them as artifacts to the run (Actions tab; kept for 30 days); a tag `v*`
attaches them to a GitHub release.  Each zip unpacks to a folder with the MEX file, `ltratio.m`, `dfx_core.m`, `tr_limdist.m`, `tr_contour.m`, `tr_limmom.m`, `tr_adfmax.m`, `tr_region.m`, `tr_gilpelaez.m` and this README.  Add the
folder to the MATLAB path.  No MATLAB is installed on the runners: the MEX file is compiled against the stub
`mexstub.h`, which declares the few functions of the MEX API it uses by their versioned symbol names (the "800"
API of R2018a and later), and on Windows is linked against import libraries generated from `win/libmx.def` and
`win/libmex.def`.

## Building locally

From the repository root, on macOS (Linux: `c++ -shared`, extension `mexa64`; add `-fPIC`):

```sh
mkdir -p build
for f in dfl_special dfl_mgf dfl_single dfl_moments dfl_double dfl_gilpelaez dfl_api dfl_threads tables \
         dfl_det_la dfl_det_design dfl_det_mgf dfl_det_path dfl_det_inner dfl_det_api dfl_det_paper2 dfl_det_route dfl_det_line dfl_det_tables dfl_kan dfl_det_kan; do
  cc -std=c11 -O2 -c c/$f.c -o build/$f.o; done
cc -std=c11 -O3 -c c/dfl_gen.c -o build/dfl_gen.o    # -O3: the l-rule of the fixed start is written to be vectorized
cc -std=c11 -O3 -c c/dfl_batch.c -o build/dfl_batch.o    # (and the sums over the nodes of the fixed layouts)
c++ -std=c++11 -O2 -c c/faddeeva/Faddeeva.cc -o build/Faddeeva.o
c++ -std=c++11 -O2 -c c/faddeeva/faddeeva_shim.cpp -o build/faddeeva_shim.o
cc -std=c11 -O2 -DDFLIMIT_STUB -c c/mex/ltratio_mex.c -o build/ltratio_mex.o
c++ -bundle -undefined dynamic_lookup -o ltratio_mex.mexmaca64 build/*.o
```

On Windows, with the Visual Studio tools: `lib /def:c\mex\win\libmx.def /machine:x64 /out:libmx.lib` (and
`libmex.def`), then `cl /LD /O2 /std:c11 /DDFLIMIT_STUB /Fe:ltratio_mex.mexw64 c\mex\ltratio_mex.c c\*.c
c\faddeeva\Faddeeva.cc c\faddeeva\faddeeva_shim.cpp libmx.lib libmex.lib` (leave out the `dfl_test_*.c` files).
With a MATLAB installation, compiling against its headers instead (drop `-DDFLIMIT_STUB`, add
`-I$MATLAB/extern/include`) checks that the stub still matches.  The engine's own tests (`c/dfl_test_*.c` against
the Julia reference values in `c/test/`, written by `c/test/*_ref.jl`) are run by the workflow "C engine"
(`c.yml`), with `dfl_test_app` (the applications, the linear form U - aV and the power envelope, against the functions of `DFExact.jl`); `c/build.sh` builds
them all; `-DDFL_STRICT_COMPLEX` compiles the engine with complex numbers as a struct, which catches any use of
the complex arithmetic operators that MSVC lacks.  `c/tables.c` (the K_0 cells and the series) is generated by
`c/gen_tables.jl` from the Julia engine.

## The devices of Bao and Kan

The coefficient statistic (`ltratio('coefcdf' | 'coefpdf' | 'coeficdf' | 'coefmom', ...)`), the moments, the linear form
and the envelope take by default the devices of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments
of the Least Squares Estimator in First-Order Autoregressions: O(1) Algorithms" (`c/dfl_kan.c`, the router's rules of the
Julia package; see the top-level README): one inversion of the closed-form characteristic function along the ray of the
half-plane, the moments from one jet per node.  `opts.method` `'inversion'` (`'jets'`) forces them, `'engine'` gives the
computation of before, `opts.tol` and `opts.verified` as in the MATLAB programs, and `[v, err] = ltratio(...)` returns the
error estimates (NaN where the engines computed).  `ltratio_mex('route', ...)` takes the method as the C codes of
`dfl_inv_opts` (0, 5, 6, and 1 to 4 for the representations of the coefficient statistic), and `ltratio_mex('det', ...)` the
point-optimal functions and `coefgp` (`coef_cdf_gp`) with method 0, 5 or 7 (the fixed grid), tol, N and umax.
