# News

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
