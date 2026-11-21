# dflimit-matlab

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24404382.svg)](https://doi.org/10.5281/zenodo.24404382)

The exact limiting distribution of the Dickey-Fuller t-ratio in MATLAB, and the C engine of the R and Python packages.

The programs compute the distribution function, density, quantiles and moments of the limits of the Dickey-Fuller
t-ratios under local alternatives, phi = 1 + c/n, by the exact representations of

> Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio and Some Related
> Statistics", working paper, 2026, available at <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619>.

for the regression without deterministic terms, with an intercept and with an intercept and a trend; with a zero, a
fixed or a stationary initial observation; for the t-ratio centred at the unit root or at a local null value c0; and for
the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996). They also have the ADFmax test of Leybourne (1995),
the coefficient statistic, the probabilities of regions {U <= g(V)} of the limit experiment (point-optimal and
weighted-average-power tests), the distribution of the linear form U - aV and the asymptotic power envelope.

```matlab
addpath /path/to/dflimit-matlab/matlab
ltratio_icdf([0.01 0.05 0.10], 'DU')         % asymptotic critical values of tau_mu
ltratio_cdf(-2.86157, 'DU', -5)               % limiting power of the 5% test at c = -5
ltratio_cdf(-2.86157, 'SU', -5)               % the same with a stationary start
ltratio_cdf(0, 'DK', 2, struct('c0', 2))      % t-ratio centred at the true value
ltratio_cdf(-2, 'DU', -5, struct('g0', 1.5))  % fixed start
ltratio_mom(4, 'DT')                          % first four moments of tau_tau
```

## Two engines

The repository holds two implementations of the same functions with the same arguments.

The MATLAB programs in `matlab/` are plain MATLAB (R2019b or later; they also run, more slowly, on GNU Octave). They
take 0.05 to 0.5 seconds per value on one thread and are the readable form of the method: `Contents.m` lists them and
the manual `matlab/doc/ltratio_manual.pdf` describes the user-level routines, the engines (the line integral `dfx_*`,
the conditional integral `tr_limdist`, the contour integral `tr_contour`), accuracy and cost.

The C engine in `c/` is a port of the Julia engine `DFExact.jl` of the package DFLimit.jl, multithreaded and 20 to 50
times faster. In MATLAB it is reached through the MEX file `ltratio_mex`, built from `c/` with the gateway
`c/mex/ltratio_mex.c`. The folder of the MEX file has `ltratio.m`, which offers the functions of the MATLAB programs
with the same arguments (`ltratio('cdf', -2.86157, 'DU', -5)`), and drop-in replacements of the engines of the MATLAB
programs (`dfx_core.m`, `tr_limdist.m`, `tr_contour.m`, `tr_limmom.m`, `tr_adfmax.m`, `tr_region.m`,
`tr_gilpelaez.m`): with that folder ahead of `matlab/` on the path, `ltratio_cdf` and the other user-level routines run
unchanged on the C engine. The C engine also has the line integral with a fixed start, which the MATLAB programs lack
(they take the conditional integral there).

Use the MATLAB programs when nothing but MATLAB may be installed or when the code is to be read, and the MEX file for
speed: power curves, confidence belts, anything with many values. Both give the same values to well below 1e-5.

## Installation

Clone the repository (or download it as a zip) and add the folder `matlab` to the MATLAB path:
`addpath /path/to/dflimit-matlab/matlab` (and `savepath` to keep it).

The MEX file needs nothing else: no compiler, no runtime, one file of about 2.5 MB. There are two ways to get it.

1. Download a built one. The workflow "MATLAB MEX" (`.github/workflows/mex.yml`), which runs on a push that changes `c/`, on a tag `v*`
   and by hand from the Actions tab ("Run workflow"), builds one zip per platform (`ltratio-linux-x86_64.zip`, `ltratio-macos.zip`
   with the Apple silicon and the Intel MEX files, `ltratio-windows-x86_64.zip`) and attaches them to the run as
   artifacts; on a tag it also attaches them to a GitHub release of the repository
   (https://github.com/s-broda/dflimit-matlab/releases). Unzip, and add the folder to the
   path ahead of `matlab/`: `addpath ltratio-macos`.
2. Build it yourself, on macOS or Linux from the repository root with a C and a C++ compiler, or on Windows with the
   Visual Studio tools; `c/mex/README.md` has the commands, the timings and the details of the drop-ins.

## Tests

The C engine: `sh c/build.sh` compiles the engine and its test programs in `c/`; each test compares the engine with the
reference values of the Julia engine in `c/test/`:

```sh
cd c
for t in prim single gen fixed mom dbl app; do ./dfl_test_$t test/${t}_ref.txt; done
./dfl_test_vmath
./dfl_test_bounds
```

Every program ends with a line that counts the checks and the failures (`dfl_test_fixed` also reads
`test/contour_ref.txt`, the reference values of the contour integral; `dfl_test_bounds` has its reference values in
the source: the points of the repairs of 0.5.0, where the line integral exceeds the bounds of the engine or, with a
fixed start, climbed). The workflow "C engine" (`.github/workflows/c.yml`)
runs the same tests on Linux, macOS and Windows, on a push that changes `c/`, on a tag `v*` and by hand.

The MATLAB programs: `ltratio_test` runs 268 checks in about five minutes and returns the number of failures
(`ltratio_test(true)` is the quick version, `ltratio_test(false, true)` raises an error on a failure); `dfx_test` checks
the line integral against stored values of the Julia engine (`matlab/validation/python_exact_reference.json`). Start
MATLAB in a folder other than `matlab/`. With the folder of the MEX file ahead of `matlab/` on the path the same
`ltratio_test` runs on the C engine, in about 30 seconds, which is the check of a built MEX file.

`matlab/validation/` has the script behind the accuracy tables of earlier versions of the Supplementary Material
(`acc_table.m`), the timings behind the automatic choice of the representation (`method_sweep.m` with its logs) and a
record of older computations (`tratio_legacy/`).

The reference values in `c/test/*_ref.txt` and the tables `c/tables.c` are written from the Julia engine by
`c/test/*_ref.jl` and `c/gen_tables.jl`, in the Julia environment `c/test/Project.toml`, which loads DFLimit.jl from
its GitHub repository: once `julia --project=c/test -e 'using Pkg; Pkg.instantiate()'`, then for example
`julia --project=c/test -t auto c/test/single_ref.jl c/test/single_ref.txt` from the repository root.

## The C engine upstream of R and Python

The C sources in `c/` are the upstream of the R package [`dflimit`](https://github.com/s-broda/dflimit) and of the
Python package dflimit-py: each copies the engine sources of `c/` and `c/faddeeva/` (not the test programs, `c/test/`
or `c/mex/`) into its own tree with its `tools/sync-engine.sh /path/to/dflimit-matlab` and records the commit it copied from. A change to the
engine is made here, tested here, and then synced into the two packages. The Julia engine, of which the C engine is a
port, lives in the Julia package [`DFLimit.jl`](https://github.com/s-broda/DFLimit.jl), which also has the same
functions with the same names and arguments.

The scripts that reproduce the tables and figures of the paper are in the paper's replication package, which runs
them on these MATLAB programs and on DFLimit.jl.

## Licence and citation

MIT (see `LICENSE`); the Faddeeva package in `c/faddeeva/` is by Steven G. Johnson, MIT licence (see the header of
`c/faddeeva/Faddeeva.cc`). Please cite the paper when using the programs.

The code was developed in the paper's repository up to commit b484d4f (see `NEWS.md`).
