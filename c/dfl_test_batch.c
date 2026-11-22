/* dfl_test_batch.c: the fixed node layouts of the line integral (dfl_batch.c) against the adaptive rules, on the points of
 * the single integral's reference grid (test/single_ref.txt: cases DK, DU, DT, DG, zero and stationary start, c from 10 to
 * -480): F and f with opts.batched on and off must agree to 1e-12 (the density relative to max(1, f)) for c > -20, and to
 * 1e-11 for c <= -20, the agreement of the layouts with the adaptive rules and with layouts of twice the nodes on
 * [-480, -20] (k0first of DFExact.jl; the Julia engine's layouts differ from its adaptive rules by up to 3.2e-12 at the
 * same points, cases DU, DT at c = -20 to -480, as these); and the timing of both (one pair of values per point) */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/single_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    dfl_single_opts on, off; dfl_single_defaults(&on); dfl_single_defaults(&off); off.batched = false;
    int model, st, n = 0, nb = 0, fails = 0; double c, x, F, fd, wF = 0, wF20 = 0, wf = 0, tb = 0, ta = 0;
    while (fscanf(f, "%d %lf %d %lf %lf %lf", &model, &c, &st, &x, &F, &fd) == 6) {
        double tau = c <= 0 && !st ? 1.0 : -1;   /* (for the count of the points the layouts take) */
        int used = dfl_usebatch(1, c, x, tau > 0 ? tau : 1.0);
        double t0 = dfl_wtime();
        double Fb = dfl_single_value_c0(2, x, model, c, st, &on), fb = dfl_single_value_c0(1, x, model, c, st, &on);
        double t1 = dfl_wtime();
        double Fa = dfl_single_value_c0(2, x, model, c, st, &off), fa = dfl_single_value_c0(1, x, model, c, st, &off);
        double t2 = dfl_wtime();
        tb += t1 - t0; ta += t2 - t1; n++; nb += used;
        double eF = fabs(Fb - Fa), ef = fabs(fb - fa) / fmax(1.0, fabs(fa));
        if (eF > wF) wF = eF; if (ef > wf) wf = ef; if (c > -20 && eF > wF20) wF20 = eF;
        double tol = c > -20 ? 1e-12 : 1e-11;
        if (!(eF <= tol) || !(ef <= tol)) { fails++; printf("FAIL m=%d c=%g stat=%d x=%g: F %.15f vs %.15f (%.1e)  f %.15e vs %.15e (%.1e)\n", model, c, st, x, Fb, Fa, eF, fb, fa, ef); }
    }
    printf("%d points (%d in the range of the layouts), %d failed, largest |F_batched - F_adaptive| %.1e (c > -20: %.1e), |f_batched - f_adaptive| %.1e; %.2f ms per pair batched, %.2f adaptive (%d threads)\n",
           n, nb, fails, wF, wF20, wf, 1e3 * tb / n, 1e3 * ta / n, dfl_nworkers());
    return fails ? 1 : 0;
}
