/* dfl_test_app.c: the applications (dfl_adfmax, dfl_region, dfl_lincdf, dfl_power_envelope) against the Julia values of
 * julia/DFExact.jl (test/app_ref.txt, written by test/app_ref.jl): F and f of the ADFmax statistic to 1e-10, the probability
 * of a region to 1e-10, the distribution function of U - aV and of V and the power envelope to 1e-10 (all the points of a
 * line in one call, on the threads over the points, and the first point alone, on the threads within the point; the two
 * must agree exactly), and the limits of the distribution function at -Inf, Inf and NaN */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/app_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    char tag[16], err[256]; int n = 0, fails = 0; double w = 0, tt = 0;
    while (fscanf(f, "%15s", tag) == 1) {
        double t0 = dfl_wtime(), d = NAN;
        if (!strcmp(tag, "adfmax")) {
            int model, st; double c, x, F, fd, Fc, fc; dfl_double_info i1, i2;
            if (fscanf(f, "%d %d %lf %lf %lf %lf", &model, &st, &c, &x, &F, &fd) != 6) break;
            int rc = dfl_adfmax(2, 1, &x, model, c, st, NULL, &Fc, &i1, err, sizeof err) | dfl_adfmax(1, 1, &x, model, c, st, NULL, &fc, &i2, err, sizeof err);
            d = fmax(fabs(Fc - F), fabs(fc - fd));
            if (rc || !(d < 1e-10) || i1.status || i2.status) { fails++; printf("FAIL adfmax m=%d stat=%d c=%g x=%g: F %.15f vs %.15f, f %.15f vs %.15f (status %d %d) %s\n", model, st, c, x, Fc, F, fc, fd, i1.status, i2.status, rc ? err : ""); }
        } else if (!strcmp(tag, "region")) {
            int type; double p[2], v0, c, P, Pc; dfl_double_info i1;
            if (fscanf(f, "%d %lf %lf %lf %lf %lf", &type, &p[0], &p[1], &v0, &c, &P) != 6) break;
            int rc = dfl_region(type, p, v0, 0, NULL, c, NULL, &Pc, &i1, err, sizeof err);
            d = fabs(Pc - P);
            if (rc || !(d < 1e-10) || i1.status) { fails++; printf("FAIL region type=%d p=(%g, %g) v0=%g c=%g: %.15f vs %.15f (status %d) %s\n", type, p[0], p[1], v0, c, Pc, P, i1.status, rc ? err : ""); }
        } else if (!strcmp(tag, "lincdf") || !strcmp(tag, "vcdf") || !strcmp(tag, "envelope")) {
            int model = 1, k = 0, lin = !strcmp(tag, "lincdf"), env = !strcmp(tag, "envelope"); double a = NAN, c = 0, cbar = DFL_CBAR_GLS, alpha = 0.05;
            if (lin && fscanf(f, "%d %lf %lf %lf %d", &model, &a, &c, &cbar, &k) != 5) break;
            if (!lin && !env && fscanf(f, "%lf %d", &c, &k) != 2) break;
            if (env && fscanf(f, "%d %lf %d", &model, &alpha, &k) != 3) break;
            if (k < 1 || k > 64) break;
            double x[67], P[64], Pc[67], P1; int ok = 1;
            for (int j = 0; j < k; j++) if (fscanf(f, "%lf %lf", &x[j], &P[j]) != 2) ok = 0;
            if (!ok) break;
            int rc;
            if (env) rc = dfl_power_envelope(k, x, model, alpha, NULL, Pc, err, sizeof err) | dfl_power_envelope(1, x, model, alpha, NULL, &P1, err, sizeof err);
            else { x[k] = -INFINITY; x[k + 1] = NAN; x[k + 2] = INFINITY;
                   rc = dfl_lincdf(k + 3, x, a, model, c, cbar, NULL, Pc, err, sizeof err) | dfl_lincdf(1, x, a, model, c, cbar, NULL, &P1, err, sizeof err);
                   if (!rc && !(Pc[k] == 0 && isnan(Pc[k + 1]) && Pc[k + 2] == 1)) { fails++; printf("FAIL %s: the limits at -Inf, NaN, Inf are %g %g %g\n", tag, Pc[k], Pc[k + 1], Pc[k + 2]); } }
            d = 0; for (int j = 0; j < k; j++) d = fmax(d, fabs(Pc[j] - P[j]));
            if (rc || !(d < 1e-10) || P1 != Pc[0]) {
                fails++; printf("FAIL %s model=%d a=%g c=%g cbar=%g alpha=%g: %s\n", tag, model, a, c, cbar, alpha, rc ? err : "");
                for (int j = 0; j < k; j++) printf("    at %g: %.15f vs %.15f\n", x[j], Pc[j], P[j]);
                if (!rc) printf("    the first point alone: %.17g, in the array %.17g\n", P1, Pc[0]);
            }
            n += k - 1;
        } else { printf("unknown tag %s\n", tag); return 1; }
        tt += dfl_wtime() - t0; n++; if (d > w) w = d;
    }
    printf("%d points, %d failed, worst difference %.1e, %.0f ms per point (%d threads)\n", n, fails, w, 1e3 * tt / n, dfl_nworkers());
    return fails ? 1 : 0;
}
