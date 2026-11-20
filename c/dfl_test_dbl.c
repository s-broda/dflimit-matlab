/* dfl_test_dbl.c: the double integral against the Julia values (test/dbl_ref.txt) */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/dbl_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    /* verified = 1: the stored values at the tolerances of the present behaviour; and the one-evaluation default (o1) against
     * them at the tolerance it meets (5e-7 in F, 1e-5 absolute or 1e-6 relative in f), except in the regimes named by
     * dfl_brittle, where it is not expected to */
    dfl_double_opts o, o1; dfl_double_defaults(&o); o.verified = 1; dfl_double_defaults(&o1); char err[256];
    int model, st, n = 0, fails = 0, n1 = 0; double c, c0, g0, x, F, fd, wF = 0, wf = 0, w1 = 0, tt = 0, tt1 = 0;
    while (fscanf(f, "%d %lf %lf %lf %d %lf %lf %lf", &model, &c, &c0, &g0, &st, &x, &F, &fd) == 8) {
        double Fc, fc; dfl_double_info iF, ifd;
        double t0 = dfl_wtime();
        int rc = dfl_double(2, 1, &x, model, c, c0, g0, st, &o, &Fc, &iF, err, sizeof err);
        rc |= dfl_double(1, 1, &x, model, c, c0, g0, st, &o, &fc, &ifd, err, sizeof err);
        tt += dfl_wtime() - t0;
        n++;
        if (rc) { fails++; printf("FAIL m=%d c=%g c0=%g g0=%g stat=%d x=%g: %s\n", model, c, c0, g0, st, x, err); continue; }
        double eF = fabs(Fc - F), ef = fabs(fc - fd) / fmax(fabs(fd), 1e-300), ea = fabs(fc - fd);
        if (eF > wF) wF = eF; if (ea > wf) wf = ea;
        if (!(eF < 1e-10) || !(ef < 1e-8 || ea < 1e-10) || iF.status || ifd.status) { fails++; printf("FAIL m=%d c=%g c0=%g g0=%g stat=%d x=%g: F %.15f vs %.15f (%.1e)  f %.15e vs %.15e (%.1e)  status %d %d err %.1e %.1e\n", model, c, c0, g0, st, x, Fc, F, eF, fc, fd, ef, iF.status, ifd.status, iF.err, ifd.err); }
        if (dfl_brittle(2, model, c, c0, g0, st, 2) == NULL) {
            double F1, f1; t0 = dfl_wtime();
            rc = dfl_double(2, 1, &x, model, c, c0, g0, st, &o1, &F1, &iF, err, sizeof err) | dfl_double(1, 1, &x, model, c, c0, g0, st, &o1, &f1, &ifd, err, sizeof err);
            tt1 += dfl_wtime() - t0; n1++;
            double e1 = fmax(fabs(F1 - F), fmin(fabs(f1 - fd), 10 * fabs(f1 - fd) / fmax(fabs(fd), 1e-300))); if (e1 > w1) w1 = e1;
            if (rc || !(fabs(F1 - F) < 5e-7) || !(fabs(f1 - fd) < 1e-5 || fabs(f1 - fd) < 1e-6 * fabs(fd))) { fails++; printf("FAIL (one evaluation) m=%d c=%g c0=%g g0=%g stat=%d x=%g: F %.15f vs %.15f, f %.15e vs %.15e\n", model, c, c0, g0, st, x, F1, F, f1, fd); }
        }
    }
    printf("%d points, %d failed, worst |dF| %.1e, worst |df| %.1e (verified), %.0f ms per pair of values (%d threads); one evaluation at %d of them, worst %.1e, %.0f ms per pair\n",
           n, fails, wF, wf, 1e3 * tt / n, dfl_nworkers(), n1, w1, n1 ? 1e3 * tt1 / n1 : 0.0);
    /* the double zero x = -sqrt(2 c0) (c = -2.5, c0 = 2, x = -2, November 2026): the conditional integral returned NaN with
     * status 0 with a stationary start and an infinite density in cases DU, DT; the default (method 0) must be finite in the
     * six cases (the density of case DK is infinite at x = -2), its distribution function within 1e-9 of the line integral's,
     * its density within 1e-4 (the density has a cusp there: just below x = -2 the line integral's density is off by up to
     * 5e-5, the conditional integral's stable to 1e-13; at x = -2 both are good to about 1e-6 only) */
    {
        const double xs[4] = {-2.0, -2.0 + 1e-9, -2.000001, -1.999}; double wdF = 0, wdf = 0; int bad = 0;
        for (int m = 1; m <= 3; m++) for (int stt = 0; stt <= 1; stt++) for (int k = 0; k < 4; k++) {
            double Fa, fa, Fs, fs; int ns;
            int rc = dfl_dist(2, 1, &xs[k], m, -2.5, 2.0, 0.0, stt, 0, NULL, &o, &Fa, &ns, err, sizeof err) | dfl_dist(1, 1, &xs[k], m, -2.5, 2.0, 0.0, stt, 0, NULL, &o, &fa, &ns, err, sizeof err)
                   | dfl_dist(2, 1, &xs[k], m, -2.5, 2.0, 0.0, stt, 1, NULL, NULL, &Fs, &ns, err, sizeof err) | dfl_dist(1, 1, &xs[k], m, -2.5, 2.0, 0.0, stt, 1, NULL, NULL, &fs, &ns, err, sizeof err);
            bool dkinf = m == 1 && !stt && k == 0;
            bool ok = rc == 0 && isfinite(Fa) && (dkinf ? isinf(fa) : isfinite(fa)) && fabs(Fa - Fs) < 1e-9 && (dkinf || fabs(fa - fs) < 1e-4 * fs);
            if (fabs(Fa - Fs) > wdF) wdF = fabs(Fa - Fs); if (!dkinf && fabs(fa - fs) / fs > wdf) wdf = fabs(fa - fs) / fs;
            if (!ok) { bad++; printf("FAIL double zero m=%d stat=%d x=%.10g: auto F %.15g f %.15g, line integral F %.15g f %.15g\n", m, stt, xs[k], Fa, fa, Fs, fs); }
        }
        printf("the double zero x = -2 (c = -2.5, c0 = 2), six cases: %d failed, auto vs line integral |dF| %.1e, |df/f| %.1e\n", bad, wdF, wdf);
        fails += bad;
    }
    /* the range of c (November 2026): below DFL_CMIN = -550 dfl_dist, dfl_icdf and dfl_moments return an error with DFL_CMIN_MSG */
    {
        double x = -20, v, p = 0.05, m[2], me[2]; int ns, bad = 0;
        if (!(dfl_dist(2, 1, &x, 2, -551, 0, 0, false, 0, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_CMIN_MSG) == 0)) bad++;
        if (!(dfl_dist(1, 1, &x, 1, -700, 0, 0, true, 1, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_CMIN_MSG) == 0)) bad++;
        if (!(dfl_moments(3, -551, 0, 0, false, 2, m, me, err, sizeof err) == 1 && strcmp(err, DFL_CMIN_MSG) == 0)) bad++;
        if (!(dfl_icdf(1, &p, 4, -551, 0, 0, false, 0, 1e-10, NULL, NULL, &v, NULL, NULL, NULL, err, sizeof err) == 1 && strcmp(err, DFL_CMIN_MSG) == 0)) bad++;
        x = -17; if (dfl_dist(2, 1, &x, 2, -550, 0, 0, false, 0, NULL, NULL, &v, &ns, err, sizeof err) != 0) bad++;
        /* the region of a fixed start (task 8): an error outside it in each case, inside the values of test/fixed_ref.txt */
        double xs3[3] = {-2, -2.86, -3}, g3[3] = {1, 1, 2}, r3[3] = {0.71763193275341519, 0.22106021106163276, 0.2889970388707237};
        for (int m = 1; m <= 3; m++) { if (dfl_dist(2, 1, &xs3[m - 1], m, -5, 0, g3[m - 1], false, 0, NULL, NULL, &v, &ns, err, sizeof err) != 0 || !(fabs(v - r3[m - 1]) < 1e-11)) bad++; }
        x = -2;
        if (!(dfl_dist(2, 1, &x, 1, 0, 0, 2.5, false, 0, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_FIXED_MSG) == 0)) bad++;
        if (!(dfl_dist(1, 1, &x, 2, -10, 0, 2, false, 0, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_FIXED_MSG) == 0)) bad++;
        if (!(dfl_moments(3, -110, 0, 0.5, false, 2, m, me, err, sizeof err) == 1 && strcmp(err, DFL_FIXED_MSG) == 0)) bad++;
        printf("the range of c (at least %g) and the region of a fixed start: %d failed\n", DFL_CMIN, bad);
        if (bad) printf("FAIL the range check of c\n");
        fails += bad;
    }
    /* the stationary start near the unit root (November 2026): for c > DFL_CSTAT = -1e-4 (and c > 0) dfl_dist, dfl_single,
     * dfl_double, dfl_icdf, dfl_moments and dfl_adfmax return an error with DFL_STAT_MSG; at c = -1e-4 the values of the
     * Julia engine (quickcheck.jl, section 18b), and the coefficient statistic of case SK at c = -0.001, which the lower end
     * E V e^{-8} of the double integral had put at 6e-6 (it is 0.0010438) */
    {
        double x = -1.94087, x0 = 0, v, p = 0.05, m[2], me[2]; int ns, bad = 0; dfl_double_info inf; dfl_double_opts dd; dfl_double_defaults(&dd);
        if (!(dfl_dist(2, 1, &x, 1, -9e-5, 0, 0, true, 0, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_single(1, 1, &x, 2, -9e-5, 0, 0, true, NULL, &v, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_double(2, 1, &x, 1, -9e-5, 0, 0, true, &dd, &v, &inf, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_icdf(1, &p, 3, -9e-5, 0, 0, true, 0, 1e-10, NULL, NULL, &v, NULL, NULL, NULL, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_moments(2, -9e-5, 0, 0, true, 2, m, me, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_adfmax(2, 1, &x, 2, -9e-5, true, NULL, &v, &inf, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        if (!(dfl_dist(2, 1, &x, 2, 1, 0, 0, true, 0, NULL, NULL, &v, &ns, err, sizeof err) == 1 && strcmp(err, DFL_STAT_MSG) == 0)) bad++;
        double r[3]; const char *lab[3] = {"SK F(-1.94087)", "SK coef F(-8.039)", "SU ADFmax F(-2.43542)"};
        const double ref[3] = {0.028490044498374834, 0.00032997244558452543, 0.050000718603976146}, xm = -2.43542;
        int rc = dfl_dist(2, 1, &x, 1, -1e-4, 0, 0, true, 0, NULL, NULL, &r[0], &ns, err, sizeof err)
               | dfl_dist(2, 1, &x0, 1, -1e-4, -8.039, 0, true, 0, NULL, NULL, &r[1], &ns, err, sizeof err)
               | dfl_adfmax(2, 1, &xm, 2, -1e-4, true, NULL, &r[2], &inf, err, sizeof err);
        for (int k = 0; k < 3; k++) if (rc || !(fabs(r[k] - ref[k]) < 1e-11)) { bad++; printf("FAIL c = -1e-4, %s: %.17g vs %.17g %s\n", lab[k], r[k], ref[k], rc ? err : ""); }
        double cf;
        if (dfl_dist(2, 1, &x0, 1, -1e-3, -8.039, 0, true, 0, NULL, NULL, &cf, &ns, err, sizeof err) != 0 || !(fabs(cf - 0.0010438037922007037) < 1e-12)) { bad++; printf("FAIL SK c = -0.001 coef F(-8.039) = %.17g\n", cf); }
        printf("the stationary start near the unit root (c <= %g): %d failed\n", DFL_CSTAT, bad);
        if (bad) printf("FAIL the range check of the stationary start\n");
        fails += bad;
    }
    return fails ? 1 : 0;
}
