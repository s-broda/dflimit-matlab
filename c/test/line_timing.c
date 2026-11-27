/* line_timing.c: the measurements behind the rule of method auto of the general engine (det_line_auto, dfl_det_line.c): the
 * line integral with the explicit kernel (method line) against the inversion on paths around the roots (method paths) for
 * designs whose numerator has rank one or two, warm calls (the minimum of 7 after one call; the outer path, the provider and
 * the layouts cached): per value (F and f at x = mean + (-2, -1, 0, 1, 2) sd, ten values) and per quantile (p = 0.05 from the
 * normal approximation; each repetition at a new p, 0.05 + 0.001 k, so that the points of the search are new), 1 and 4
 * threads; and the largest difference of the two methods' values.
 *   cc -std=c11 -O2 -I c -c c/test/line_timing.c ... (link with the objects of c/build.sh); ./line_timing [designs] */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static char E[512];
static void trigf(double r, double *out, void *ctx) { (void) ctx; out[0] = 1; out[1] = sin(2 * M_PI * r); out[2] = cos(2 * M_PI * r); }
static int design(const char *nm, dfl_design *D)
{
    if (!strcmp(nm, "trig")) return dfl_design_fdesign(trigf, NULL, 3, 0, NULL, 1e-14, 96, 0, D, E, sizeof E);
    return dfl_design_trend(atoi(nm + 1), D, E, sizeof E);
}
static int vals(const dfl_design *D, double c, dfl_det_opts *o, const double *xs, double *v)
{
    return dfl_det_dist(2, D, 5, xs, c, o, v, E, sizeof E) || dfl_det_dist(1, D, 5, xs, c, o, v + 5, E, sizeof E);
}
int main(int argc, char **argv)
{
    const char *dd[] = {"t0", "t1", "t2", "t3", "t4", "trig"}; int nd = 6;
    const char **ds = argc > 1 ? (const char **) argv + 1 : dd; if (argc > 1) nd = argc - 1;
    double cs[] = {2, 0, -5, -13.5, -20, -30, -50};
    printf("the line integral (method line) against the inversion on paths (method paths): ms per value (F and f at five points) and per quantile\n");
    for (int id = 0; id < nd; id++) for (int st = 0; st <= 1; st++) for (int ic = 0; ic < 7; ic++) {
        double c = cs[ic]; if (st && c >= 0) continue; if (st && strcmp(ds[id], "t2")) continue;
        dfl_design D; if (design(ds[id], &D)) { printf("%s\n", E); return 1; }
        dfl_det_opts ol, op; dfl_det_defaults(&ol); dfl_det_defaults(&op);
        ol.method = DFL_DET_LINE; op.method = DFL_DET_PATHS; ol.start = op.start = st ? DFL_START_STAT : DFL_START_ZERO;
        double m[2]; if (dfl_det_moments(&D, c, 2, &op, m, NULL, E, sizeof E)) { printf("%s\n", E); continue; }
        double sd = sqrt(m[1] - m[0] * m[0]), xs[5], vl[10], vp[10], x0 = m[0] - 1.6448536269514722 * sd;
        for (int k = 0; k < 5; k++) xs[k] = m[0] + (k - 2) * sd;
        if (vals(&D, c, &ol, xs, vl)) { printf("%-5s %s c = %6.1f: line: %s\n", ds[id], st ? "stat" : "zero", c, E); continue; }
        if (vals(&D, c, &op, xs, vp)) { printf("paths: %s\n", E); continue; }
        double d = 0; for (int k = 0; k < 10; k++) d = fmax(d, fabs(vl[k] - vp[k]));
        printf("%-5s %s c = %6.1f  diff %.1e", ds[id], st ? "stat" : "zero", c, d);
        for (int nth = 1; nth <= 4; nth += 3) {
            dfl_set_nworkers(nth);
            double tl = INFINITY, tp = INFINITY, ql = INFINITY, qp = INFINITY, q;
            for (int r = 0; r < 7; r++) {
                double t = dfl_wtime(); vals(&D, c, &ol, xs, vl); t = dfl_wtime() - t; if (t < tl) tl = t;
                t = dfl_wtime(); vals(&D, c, &op, xs, vp); t = dfl_wtime() - t; if (t < tp) tp = t;
            }
            double p0 = 0.049; dfl_det_quantile(&D, 1, &p0, c, x0, 1e-10, &ol, &q, E, sizeof E); dfl_det_quantile(&D, 1, &p0, c, x0, 1e-10, &op, &q, E, sizeof E);
            for (int r = 0; r < 7; r++) {
                double p = 0.05 + 0.001 * r + 0.0001 * nth;
                double t = dfl_wtime(); dfl_det_quantile(&D, 1, &p, c, x0, 1e-10, &ol, &q, E, sizeof E); t = dfl_wtime() - t; if (t < ql) ql = t;
                t = dfl_wtime(); dfl_det_quantile(&D, 1, &p, c, x0, 1e-10, &op, &q, E, sizeof E); t = dfl_wtime() - t; if (t < qp) qp = t;
            }
            printf("  | %d thr: value line %6.2f paths %6.2f (x%4.2f)  quantile line %6.1f paths %6.1f (x%4.2f)", nth, 1e2 * tl, 1e2 * tp, tp / tl, 1e3 * ql, 1e3 * qp, qp / ql);
        }
        printf("\n"); fflush(stdout);
        dfl_design_free(&D);
        dfl_set_nworkers(4);
    }
    return 0;
}
