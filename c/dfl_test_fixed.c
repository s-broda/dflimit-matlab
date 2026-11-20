/* dfl_test_fixed.c: the single integral with a fixed start g0 != 0 (the noncentral kernel) against the Julia values
 * (first argument, default test/fixed_ref.txt), and the contour integral of a fixed start in cases DU, DT (dfl_contour)
 * against the Julia values of tr_contour (second argument, default test/contour_ref.txt; columns model c c0 g0 x F f
 * err_F err_f): absolute 1e-12 in F and f (or the reference's error estimate where it is larger), status 0; method 0 of dfl_dist at some of the points with F, f >= 1e-3 (no tail
 * recomputation there) must return the contour integral's values, at one tail point the double integral's; the scope of
 * method 3. */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/fixed_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    dfl_single_opts o; dfl_single_defaults(&o);
    if (getenv("DFL_NLNC")) DFL_NLNC = atoi(getenv("DFL_NLNC"));          /* the l-rule's node count, for checking */
    int model, st, n = 0, fails = 0; double c, c0, g0, x, F, fd, wF = 0, wf = 0, tt = 0;
    while (fscanf(f, "%d %lf %lf %lf %d %lf %lf %lf", &model, &c, &c0, &g0, &st, &x, &F, &fd) == 8) {
        double t0 = dfl_wtime();
        double Fc = dfl_single_value_gen(2, x, model, c, c0, g0, st, &o), fc = dfl_single_value_gen(1, x, model, c, c0, g0, st, &o);
        tt += dfl_wtime() - t0;
        double eF = fabs(Fc - F), ef = fabs(fc - fd) / fmax(fabs(fd), 1e-300), ea = fabs(fc - fd); n++;
        bool inf_ok = isinf(fd) && isinf(fc);
        if (eF > wF) wF = eF; if (ef > wf && ea > 1e-11 && !inf_ok) wf = ef;
        if (!(eF < 1e-11) || !(ef < 1e-9 || ea < 1e-11 || inf_ok)) { fails++; printf("FAIL m=%d c=%g c0=%g g0=%g stat=%d x=%g: F %.15f vs %.15f (%.1e)  f %.15e vs %.15e (%.1e)\n", model, c, c0, g0, st, x, Fc, F, eF, fc, fd, ef); }
        else if (getenv("DFL_VERBOSE")) printf("ok   m=%d c=%g c0=%g g0=%g x=%g: F %.15f (%.1e) f %.15e (%.1e)\n", model, c, c0, g0, x, Fc, eF, fc, ef);
    }
    fclose(f);
    printf("%d points, %d failed, worst |dF| %.1e, worst |df/f| %.1e, %.1f ms per pair of values (%d threads)\n", n, fails, wF, wf, 1e3 * tt / n, dfl_nworkers());

    /* the contour integral */
    f = fopen(argc > 2 ? argv[2] : "test/contour_ref.txt", "r"); if (!f) { printf("no reference file for the contour integral\n"); return 1; }
    int nc = 0, cfails = 0, nauto = 0, ntail = 0; double cF = 0, cf = 0, c1 = 0, ct = 0, eFr, efr; char msg[512];
    while (fscanf(f, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &model, &c, &c0, &g0, &x, &F, &fd, &eFr, &efr) == 9) {
        double Fc, fc, eF, ef; int sF, sf, nF, nf;
        double t0 = dfl_wtime();
        /* verified: the stored values (the present behaviour); and the one-evaluation default against them */
        int rc = dfl_contour(2, 1, &x, model, c, c0, g0, 1, &Fc, &eF, &sF, &nF, msg, sizeof msg) | dfl_contour(1, 1, &x, model, c, c0, g0, 1, &fc, &ef, &sf, &nf, msg, sizeof msg);
        double F1 = NAN, f1 = NAN;
        rc |= dfl_contour(2, 1, &x, model, c, c0, g0, 0, &F1, NULL, NULL, NULL, msg, sizeof msg) | dfl_contour(1, 1, &x, model, c, c0, g0, 0, &f1, NULL, NULL, NULL, msg, sizeof msg);
        ct += dfl_wtime() - t0; nc++;
        if (rc) { cfails++; printf("FAIL contour m=%d c=%g c0=%g g0=%g x=%g: %s\n", model, c, c0, g0, x, msg); continue; }
        double dF = fabs(Fc - F), df = fabs(fc - fd);
        if (!(dF <= cF)) cF = dF; if (!(df <= cf)) cf = df;
        /* (the tolerance 1e-12, or the reference's own error estimate where that is larger, plus 1e-15 times the
         * noncentrality (c g0)^2/2: at the corner (c g0)^2/2 = 7200, c = -40 of case DT the roots of D(., t) at |t| ~ 600
         * carry relative errors of ~1e-11 in both engines, from the quadratic through D(-1, t), D(0, t), D(1, t), which
         * the two evaluations share and the error estimate does not see, and the values differ by 5e-12) */
        double tnc = 1e-15 * (c * g0) * (c * g0) / 2;
        if (!(dF < fmax(1e-12, eFr) + tnc) || !(df < fmax(1e-12, efr) + tnc) || sF || sf) { cfails++; printf("FAIL contour m=%d c=%g c0=%g g0=%g x=%g: F %.15f vs %.15f (%.1e, status %d)  f %.15e vs %.15e (%.1e, status %d)\n", model, c, c0, g0, x, Fc, F, dF, sF, fc, fd, df, sf); }
        else if (getenv("DFL_VERBOSE")) printf("ok   contour m=%d c=%g c0=%g g0=%g x=%g: F %.15f (%.1e, err %.1e, %d nodes) f %.15e (%.1e, err %.1e)\n", model, c, c0, g0, x, Fc, dF, eF, nF, fc, df, ef);
        /* the one evaluation (A at its first level): within 1e-9 of the stored values at these points (5e-9 at the corner
         * (c g0)^2/2 = 7200, where A needs three levels for 1e-10) */
        double d1 = fmax(fabs(F1 - F), fabs(f1 - fd)); if (!(d1 <= c1)) c1 = d1;
        if (!(d1 < (tnc > 1e-12 ? 5e-9 : 1e-9))) { cfails++; printf("FAIL contour (one evaluation) m=%d c=%g c0=%g g0=%g x=%g: F %.15f vs %.15f, f %.15e vs %.15e\n", model, c, c0, g0, x, F1, F, f1, fd); }
        /* method 0 (dfl_automethod: 3 here) away from the tails returns the contour integral's values (at every fourth
         * point inside the region of a fixed start), in the tails the double integral's (one point) */
        if (nc % 4 == 1 && Fc >= 1e-3 && fc >= 1e-3 && dfl_fixed_ok(model, c, g0)) {
            double Fa = NAN, fa = NAN, Fb = NAN, fb = NAN; int ns = -1; dfl_double_opts dv; dfl_double_defaults(&dv); dv.verified = 1;
            if (dfl_automethod(model, c, c0, g0, false, x) != 3 || dfl_dist(2, 1, &x, model, c, c0, g0, false, 0, NULL, &dv, &Fa, &ns, msg, sizeof msg)
                || dfl_dist(1, 1, &x, model, c, c0, g0, false, 0, NULL, &dv, &fa, NULL, msg, sizeof msg) || Fa != Fc || fa != fc || ns != 0
                || dfl_dist(2, 1, &x, model, c, c0, g0, false, 0, NULL, NULL, &Fb, NULL, msg, sizeof msg) || dfl_dist(1, 1, &x, model, c, c0, g0, false, 0, NULL, NULL, &fb, NULL, msg, sizeof msg)
                || Fb != F1 || fb != f1) {
                cfails++; printf("FAIL dist method 0 m=%d c=%g c0=%g g0=%g x=%g: F %.17g vs %.17g, f %.17g vs %.17g (default %.17g vs %.17g, %.17g vs %.17g)\n", model, c, c0, g0, x, Fa, Fc, fa, fc, Fb, F1, fb, f1);
            }
            nauto++;
        }
        if (ntail == 0 && Fc < 1e-3 && fc < 1e-3 && dfl_fixed_ok(model, c, g0)) {
            double Fa = NAN, fa = NAN, Fd = NAN, fdd = NAN; dfl_double_info iF, i_f; dfl_double_opts dop; dfl_double_defaults(&dop); dop.verified = 1;
            int rc2 = dfl_dist(2, 1, &x, model, c, c0, g0, false, 0, NULL, &dop, &Fa, NULL, msg, sizeof msg) | dfl_dist(1, 1, &x, model, c, c0, g0, false, 0, NULL, &dop, &fa, NULL, msg, sizeof msg)
                    | dfl_double(2, 1, &x, model, c, c0, g0, false, &dop, &Fd, &iF, msg, sizeof msg) | dfl_double(1, 1, &x, model, c, c0, g0, false, &dop, &fdd, &i_f, msg, sizeof msg);
            if (rc2 || Fa != (iF.status ? Fc : Fd) || fa != (i_f.status ? fc : fdd) || !(fabs(Fa - F) < 1e-12) || !(fabs(fa - fd) < 1e-12)) {
                cfails++; printf("FAIL dist method 0 (tail) m=%d c=%g c0=%g g0=%g x=%g: F %.17g (double %.17g, contour %.17g), f %.17g (double %.17g, contour %.17g)\n", model, c, c0, g0, x, Fa, Fd, Fc, fa, fdd, fc);
            }
            ntail++;
        }
    }
    fclose(f);
    /* the scope of method 3 */
    { double xx = -2, y; const char *want = DFL_CONTOUR_MSG;
      if (!dfl_contour(2, 1, &xx, 1, -5, 0, 1, 1, &y, NULL, NULL, NULL, msg, sizeof msg) || strcmp(msg, want)) { cfails++; printf("FAIL scope: model 1 (%s)\n", msg); }
      if (!dfl_contour(2, 1, &xx, 2, -5, 0, 0, 1, &y, NULL, NULL, NULL, msg, sizeof msg) || strcmp(msg, want)) { cfails++; printf("FAIL scope: g0 = 0 (%s)\n", msg); }
      if (!dfl_dist(2, 1, &xx, 3, -5, 0, 0, true, 3, NULL, NULL, &y, NULL, msg, sizeof msg) || strcmp(msg, want)) { cfails++; printf("FAIL scope: stationary start (%s)\n", msg); }
      if (dfl_automethod(2, -5, 0, 0, false, xx) == 3 || dfl_automethod(1, -5, 0, 1, false, xx) == 3) { cfails++; printf("FAIL automethod outside the scope\n"); } }
    /* case DU with a fixed start, c0 > 3 at c > 1.5, where the contour integral fails (1e44 at c = 2, c0 = 5, g0 = 1.5):
     * method 0 takes the conditional integral, the line integral's value to 1e-9; the warning for method 3 there */
    { double xx = -8.4075, Fa, Fl;
      int rc = dfl_dist(2, 1, &xx, 2, 2, 5, 1.5, false, 0, NULL, NULL, &Fa, NULL, msg, sizeof msg) | dfl_dist(2, 1, &xx, 2, 2, 5, 1.5, false, 1, NULL, NULL, &Fl, NULL, msg, sizeof msg);
      if (rc || dfl_automethod(2, 2, 5, 1.5, false, xx) != 2 || fabs(Fa - Fl) > 1e-9 || dfl_brittle(2, 2, 2, 5, 1.5, false, 0) != NULL || dfl_brittle(2, 2, 2, 5, 1.5, false, 3) == NULL)
          { cfails++; printf("FAIL DU c=2 c0=5 g0=1.5: method 0 %.17g, line %.17g\n", Fa, Fl); } }
    printf("contour integral: %d points, %d failed, worst |dF| %.1e, worst |df| %.1e (verified), one evaluation %.1e, %d + %d checks of method 0, %.1f ms per value with the first calls (%d threads)\n",
           nc, cfails, cF, cf, c1, nauto, ntail, 1e3 * ct / (4 * nc), dfl_nworkers());
    return (fails || cfails) ? 1 : 0;
}
