acc_table.m    the tables "The two representations against each other" of earlier versions of the Supplementary Material
               (tab:acc: zero and stationary start; tab:accfixed: fixed start): the line
               integral against the conditional integral at seven points per
               configuration, with the median times; writes acc_table.log (every point)
               and acc_rows.tex (the LaTeX rows; the fixed-start rows carry gamma_0 in
               the label).  About 25 minutes on the MATLAB programs, 5 on the C MEX file
               (c/mex, used when the MEX file is there).
method_sweep.m   per-point timings of the two representations in the MATLAB programs (one
               thread; c in [-5,5], c0 in {-5,-1,1,3}, x in [-4,32], the six cases, plus
               c0 = 0), the basis of the automatic choice in ltratio_method.m (November
               2026); writes method_sweep_<case>.log; method_sweep_summary.py aggregates
               the logs and compares the rule with the C engine's
method_sweep_*.log   the run behind the automatic choice (November 2026), read by
               method_sweep_summary.py
ext_ref.txt    the reference values of the extensions of the limit (volatility regimes, a normal component, the
               seasonal statistics) of the Julia engine, written by dev/ext_ref.jl of DFLimit.jl (the same file as
               c/test/ext_ref.txt), read by ../ext_test.m
python_exact_reference.json   the twelve stored reference values of the Julia engine (computed
               on 31 October 2026), read by ../dfx_test.m
tratio_legacy/ the limit computations of the earlier, combined version of the finite-sample
               paper, kept as a record (see its README.txt)

The logs of earlier runs of acc_table, ltratio_test and ltratio_paper are in the history of the
paper's repository (DF-Limit, public/validation/, up to commit b484d4f).
