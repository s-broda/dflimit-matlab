# News

## 0.5.0

Three errors are repaired, in the C engine (and with it the MEX file and the R and Python packages) and in the MATLAB
programs, as in DFLimit.jl 0.5.0. Values outside the regions named below are unchanged.

**Large x above the unit root.** In the bulk of the distribution at c = 7 and beyond (x in the hundreds) the inner
integrals of the line integral need a Gauss-Legendre rule whose node count grows like x^2. Beyond 4096 nodes the C
engine used a rule that it does not hold and crashed: in case DT at c = 7, x = 751 a segmentation fault after 20
seconds, which ends MATLAB with the MEX file, and R or Python with their packages. The MATLAB programs built every rule
by Golub-Welsch, with its full matrix, and passed 2 GB at the same point. The MATLAB programs now compute rules above
1024 nodes by Newton's method, as the C engine does, and both refuse a rule above 4096 nodes before anything is
evaluated: the automatic method takes the conditional integral at such a point (0.990823749760806 at the point above,
in 0.2 seconds with the C engine).

**A negative local null value near zero.** With c0 < 0 close to zero and x > 0, in every case and with every start, a
root of X(v) = c0 v^2 + x v + 1/2 lies far out, near x/|c0|, and the paths of the line integral's inner integrals grow
like (x/c0)^2; the C engine stored every node, and a call took more than 2 GB within two seconds (case SK at c = c0 =
-1e-3, x = 0.05 and at c = c0 = -1e-2, x = 1; the same in the MATLAB programs). The C engine now sums the paths as it
builds them, and both refuse a path of more than two million nodes before it is evaluated. The automatic method takes
the conditional integral at such a point, with its checks, and warns where such a value is flagged (`dfl_warning`; the
warning `ltratio:bound` in MATLAB); with the line integral requested, and in case DG, which has the line integral only,
an error names the arguments (the return code 10 + the bound in C, the identifier `ltratio:bound` in MATLAB). Below the
bound the footprint of such a point fell from 1 GB to 5 MB (case SK, c = c0 = -0.1, x = 1.5).

**A piece of the path that climbed.** Where the exponent of the line integral's integrand rises into a segment at one
of its ends, the piece of steepest descent from that end could climb by many orders of magnitude and cancel against the
next piece. With a fixed start (the line integral of the C engine: case DK by default, cases DU, DT with the line
integral requested) and c0 of -0.01 or -0.001 at x > 0 and c >= 0, the values were wrong, from errors of 1e-7 and
2.4e-3 that look plausible to 6e74: in case DK at c = 0, c0 = -0.01, g0 = 1, x = 1.5 the engine returned 6.03e74 for
0.960992, at c = 2 0.08440 for 0.08676. The same fault gave -3.5e86 and 1.2e8 for 0.48 with c0 = 1 at c = -20 in cases
DU, DT with the line integral requested, and, with a zero start, 1.1e13 (the MATLAB programs 1.7e13) for 0.68267 in
case DK at c = 0, c0 = -1e-4, x = 0. Such a piece is now cut back where it climbs above the larger end of the segment.

Of 6945 values of the C engine across its functions (the six cases and DG, the three starts, c from -50 to 10, c0 of
both signs, the three methods with and without `verified`, quantiles, moments, the coefficient statistic, ADFmax, the
regions, the linear form and the power envelope) 6932 are unchanged bit for bit; the others are four of the wrong values
above, seven that move by at most 3e-13 where the descent piece is now cut back, and two in case DK at c = 5,
c0 = -0.1, x = 14.17, whose paths now exceed the bound and which the conditional integral gives within 3e-15 of
before. At the repaired points the C engine agrees with DFLimit.jl 0.5.0 within 1e-13. Of 1398 values of the MATLAB
programs across the same functions all are unchanged bit for bit, and `ltratio_test` passes, also on the MEX file. At
the flagged points near c0 = 0 the conditional integral of the MATLAB programs differs from that of the C engine by up
to 4e-9 (the density of case SK at c = c0 = -1e-3, x = 0.05, which is off by 1.2e-3 in both). The C test program
`dfl_test_bounds` checks the repaired points. The MEX file gives the warning of `dfl_warning` and passes the refusals
of the line integral to the MATLAB programs as `dfx:bound:rule` and `dfx:bound:path`, so that `ltratio_dist` routes
them as it does those of the MATLAB engine.

## 0.4.1

Documentation only; the programs and the values are those of 0.4.0. The README gives the paper's address on SSRN;
`CITATION.cff` names the paper as the suggested citation, and `.zenodo.json` describes the package for its archive on
Zenodo; the pointers to the reproduction of the paper name the paper's replication package; the workflows run again on a
push that changes `c/`.

## 0.4.0

The first release as a repository of its own. The MATLAB programs, their manual and validation material, the C
engine with its tests and build script, the MEX gateway and the workflows that test the C engine and build the MEX
file were developed in the paper's repository (DF-Limit, the directories `public/` and `c/`); the code here is that of
DF-Limit commit b484d4f, with the MATLAB programs in `matlab/` (formerly `public/`). The values are unchanged.

The changes are those of the move. The paths that named `public/` name `matlab/`; `dfx_test` reads its reference
values from `matlab/validation/`; the reference scripts `c/test/*_ref.jl` and `c/gen_tables.jl` load the Julia engine
from the package DFLimit.jl (the environment `c/test/Project.toml`); the reproduction of the paper (`ltratio_paper`,
`mackinnon94_p`) stays in the paper's repository; the logs of earlier validation runs stay in its history. The workflows
run when started by hand and on a tag `v*`, no longer on every push. `Contents.m` and the manual state the version.
The R package dflimit and the Python package dflimit-py now copy the C engine from this repository.
