/* dfl_test_bounds.c: the bounds of the line integral's quadrature and the repairs of November 2026 (version 0.5.0), against
 * the values of the Julia engine DFLimit.jl 0.5.0, which has the same repairs, and of the contour integral of a newer engine
 * (the values of the Julia tests):
 *   1. large x above the unit root (case DT, c = 7, x = 751), where the segment to the saddle of the Gaussian factor of the
 *      inner integrals needs a Gauss-Legendre rule of more than DFL_RULEMAX nodes: the engine dereferenced the missing rule
 *      and crashed; the automatic method now takes the conditional integral there, the line integral requested explicitly
 *      returns an error that names the arguments;
 *   2. c0 < 0 near zero with x > 0, where a root of X(v) lies far out and the paths of the general kernel grew without bound
 *      (more than 2 GB within two seconds at c = c0 = -1e-3, x = 0.05 with a stationary start): refused beyond DFL_NODEMAX
 *      nodes, the conditional integral with the automatic method (with the warning of dfl_warning where its value is
 *      flagged), an error otherwise (also in case DG, which has the line integral only);
 *   3. a fixed start with c0 in {-0.01, -0.001}, x > 0, c >= 0, where a descent piece of the general kernel's segments
 *      climbed and cancelled (F = 6e74 for 0.96 in case DK at c = 0, c0 = -0.01, g0 = 1, x = 1.5), and the same fault
 *      elsewhere (a fixed start with c0 = 1 at c = -20 by the line integral, a zero start at c0 = -1e-4, x = 0);
 *   4. the upper tail above the unit root (November 23, 2026): case DK with c0 != c beyond the line integral's bound, and for
 *      c > 4 at x > 300, the conditional integral with its checks (before, its flagged points got one evaluation, off by up
 *      to 2.4: values outside [0, 1]); above c = 14 such points are refused; the general engine's line cut to the rate of the
 *      stationary value with c0 (before, off by 1.1e-3 at c = 10, c0 = 9); the final check of a distribution function.
 *      The references: an independent computation, the law of V given J(1) by Girsanov's theorem inverted and integrated
 *      over J(1), within 1e-12 of both engines where they agree. */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int nfail = 0, nchk = 0;
static double worst = 0;
/* one value of dfl_dist against the reference r (tolerance tol absolute, or relative if rel); warn: whether dfl_warning is
 * expected (-1: either) */
static void val(const char *lbl, int kind, int model, double c, double c0, double g0, bool stat, int method, double x, double r, double tol, int rel, int warn)
{
    double v; int ns; char err[1024] = "";
    double t0 = dfl_wtime();
    int rc = dfl_dist(kind, 1, &x, model, c, c0, g0, stat, method, NULL, NULL, &v, &ns, err, sizeof err);
    double t = dfl_wtime() - t0, e = rc ? INFINITY : rel ? fabs(v - r) / fabs(r) : fabs(v - r);
    const char *w = dfl_warning();
    nchk++;
    if (!rc && e > worst && !rel) worst = e;
    if (rc || !(e <= tol) || (warn >= 0 && (w != NULL) != (warn != 0))) {
        nfail++; printf("FAIL %s: %s = %.15g (rc %d%s%s), reference %.15g, difference %.1e, warning %s\n", lbl, kind == 2 ? "F" : "f", v, rc, rc ? ": " : "", rc ? err : "", r, e, w ? w : "none");
    } else printf("ok   %s: %.15g (difference %.1e%s, %.2f s)%s\n", lbl, v, e, rel ? " relative" : "", t, w ? " (warning)" : "");
}
/* dfl_dist must return the error code 10 + bound with a message that contains want1 and want2 */
static void refused(const char *lbl, int kind, int model, double c, double c0, double g0, bool stat, int method, double x, int bound, const char *want1, const char *want2)
{
    double v; int ns; char err[1024] = "";
    int rc = dfl_dist(kind, 1, &x, model, c, c0, g0, stat, method, NULL, NULL, &v, &ns, err, sizeof err);
    nchk++;
    if (rc != 10 + bound || !strstr(err, want1) || (want2 && !strstr(err, want2))) { nfail++; printf("FAIL %s: rc %d, message \"%s\"\n", lbl, rc, err); }
    else printf("ok   %s: \"%s\"\n", lbl, err);
}

int main(void)
{
    double t0 = dfl_wtime();
    /* 1. large x above the unit root (the former crash) */
    val("DT c = 7 x = 751 F", 2, 3, 7, 0, 0, false, 0, 751, 0.990823749760809, 1e-12, 0, 0);
    val("DT c = 7 x = 751 f", 1, 3, 7, 0, 0, false, 0, 751, 9.2924954e-5, 1e-6, 1, 0);
    val("DT c = 7 x = 403.7 F (below the bound; far_tail: the conditional integral)", 2, 3, 7, 0, 0, false, 0, 403.7, 0.838656678280211, 1e-12, 0, 0);
    refused("DT c = 7 x = 751, the line integral requested", 2, 3, 7, 0, 0, false, 1, 751, DFL_BOUND_RULE, "the line integral needs more nodes than the engine allows", "at large x above the unit root");
    /* 2. c0 < 0 near zero with x > 0 (the unbounded paths); the values of the conditional integral, which the automatic
     * method takes there, against the Julia engine, and at x = 0.05 to 1.5 at c = c0 = -1e-4 against the contour integral of
     * the newer engine to 1e-6 (the conditional integral flags some of them: the warning) */
    val("SK c = c0 = -1e-3 x = 0.05 F", 2, 1, -1e-3, -1e-3, 0, true, 0, 0.05, 0.54749158678338294, 1e-10, 0, -1);
    val("SK c = c0 = -1e-2 x = 1 F", 2, 1, -1e-2, -1e-2, 0, true, 0, 1, 0.87075930276798486, 1e-10, 0, -1);
    val("SK c = c0 = -1e-2 x = 1 f", 1, 1, -1e-2, -1e-2, 0, true, 0, 1, 0.20596599779677932, 1e-10, 0, -1);
    {   const double xs[5] = {-1, 0.05, 0.5, 1, 1.5}, Fr[5] = {0.1666135587511895, 0.531233538152, 0.700742367308, 0.847290600876, 0.936188673511};
        double F[5]; int ns; char err[1024] = "";
        int rc = dfl_dist(2, 5, xs, 1, -1e-4, -1e-4, 0, true, 0, NULL, NULL, F, &ns, err, sizeof err);
        const char *w = dfl_warning();
        double e = 0; for (int k = 0; k < 5; k++) e = fmax(e, fabs(F[k] - Fr[k]));
        nchk++;
        if (rc || !(e < 1e-6) || !w || !strstr(w, "flags some of them")) { nfail++; printf("FAIL SK c = c0 = -1e-4, five points: rc %d (%s), worst difference %.1e, warning %s\n", rc, err, e, w ? w : "none"); }
        else printf("ok   SK c = c0 = -1e-4, five points: worst difference %.1e against the newer engine, warning \"%s\"\n", e, w);
    }
    val("DK c = 0 c0 = -1e-3 x = 1.5 F", 2, 1, 0, -1e-3, 0, false, 0, 1.5, 0.9674482411522245, 1e-10, 0, -1);
    refused("SK c = c0 = -1e-4 x = 1, the line integral requested", 2, 1, -1e-4, -1e-4, 0, true, 1, 1, DFL_BOUND_PATH, "the line integral needs more nodes than the engine allows", "c0 = -0.0001, x = 1");
    refused("DG c = -1 c0 = -1e-4 x = 1 (the line integral only)", 2, 4, -1, -1e-4, 0, false, 0, 1, DFL_BOUND_PATH, "the line integral needs more nodes than the engine allows", "case DG has the line integral representation only");
    {   double x = 1, v; char err[1024] = ""; nchk++;
        int rc = dfl_single(2, 1, &x, 1, -1e-4, -1e-4, 0, true, NULL, &v, err, sizeof err);
        if (rc != 10 + DFL_BOUND_PATH || !strstr(err, "the line integral needs more nodes")) { nfail++; printf("FAIL dfl_single at a refused point: rc %d, \"%s\"\n", rc, err); }
        else printf("ok   dfl_single at a refused point: an error\n");
    }
    val("SK c = c0 = -0.1 x = 1.5 F (the line integral, below the bound)", 2, 1, -0.1, -0.1, 0, true, 0, 1.5, 0.9577065510615822, 1e-12, 0, 0);
    /* 3. a fixed start, c0 = -0.01, x > 0: the descent piece that climbed (the values of the contour integral of the newer
     * engine; F was 6.0e74, 0.7983343610, 0.0843986044, 0.7557695044 and -6.2e29, f was 0) */
    val("DK c = 0 c0 = -0.01 g0 = 1 x = 1.5 F", 2, 1, 0, -0.01, 1, false, 0, 1.5, 0.960991862318, 1e-10, 0, 0);
    val("DK c = 0 c0 = -0.01 g0 = 1 x = 1.5 f", 1, 1, 0, -0.01, 1, false, 0, 1.5, 0.079242719297, 1e-10, 0, 0);
    val("DK c = 0 c0 = -0.01 g0 = 1 x = 0.5 F", 2, 1, 0, -0.01, 1, false, 0, 0.5, 0.798334456036, 1e-10, 0, 0);
    val("DK c = 2 c0 = -0.01 g0 = 1 x = 1.5 F", 2, 1, 2, -0.01, 1, false, 0, 1.5, 0.086764388544, 1e-10, 0, 0);
    val("DK c = 0 c0 = -0.01 g0 = 2 x = 0.5 F", 2, 1, 0, -0.01, 2, false, 0, 0.5, 0.755768817096, 1e-10, 0, 0);
    val("DU c = 2 c0 = -0.01 g0 = 1 x = 1.5 F (the line integral)", 2, 2, 2, -0.01, 1, false, 1, 1.5, 0.207535597435, 1e-10, 0, 0);
    /* the same cut of a climbing descent piece elsewhere (the values of the Julia engine; F was -3.5e86, 1.2e8, 1.1e13 and
     * f 5.1e17): a fixed start with c0 > 0 far below the unit root by the line integral, and a zero start at c0 = -1e-4,
     * x = 0 in case DK */
    val("DU c = -20 c0 = 1 g0 = 1 x = -4.7413 F (the line integral)", 2, 2, -20, 1, 1, false, 1, -4.741316631609604, 0.48079941924171321, 1e-10, 0, 0);
    val("DT c = -20 c0 = 1 g0 = 1 x = -4.5724 F (the line integral)", 2, 3, -20, 1, 1, false, 1, -4.5723644715124907, 0.48434943999136959, 1e-10, 0, 0);
    val("DK c = 0 c0 = -1e-4 x = 0 F", 2, 1, 0, -1e-4, 0, false, 0, 0, 0.68266529536039322, 1e-10, 0, 0);
    val("DK c = 0 c0 = -1e-4 x = 0 f", 1, 1, 0, -1e-4, 0, false, 0, 0, 0.32556973831816555, 1e-10, 0, 0);
    /* 4. the upper tail above the unit root (the values of the independent computation; before: 0.8415, -1.3674, 1.1422,
     * 0.6076, 1.0479, 0.8686) */
    val("DK c = 10 c0 = 2 x = 1e4 F", 2, 1, 10, 2, 0, false, 0, 1e4, 0.743623945123, 1e-10, 0, 0);
    val("DK c = 10 c0 = 2 x = 3e4 F", 2, 1, 10, 2, 0, false, 0, 3e4, 0.999338346257, 1e-10, 0, 0);
    val("DK c = 7 c0 = 2 x = 2e3 F", 2, 1, 7, 2, 0, false, 0, 2e3, 0.999999671714, 1e-10, 0, 0);
    val("DK c = 10 c0 = 9 x = 1e3 F", 2, 1, 10, 9, 0, false, 0, 1e3, 0.636120922790, 1e-10, 0, 0);
    val("DK c = 10 c0 = 9 x = 3e3 F", 2, 1, 10, 9, 0, false, 0, 3e3, 0.993550212928, 1e-10, 0, 0);
    val("DK c = 10 c0 = 2 x = 1e5 F", 2, 1, 10, 2, 0, false, 0, 1e5, 1.0, 1e-12, 0, 0);
    {   double x = 1e7, v; char err[1024] = ""; nchk++;
        int rc = dfl_dist(2, 1, &x, 1, 16, 2, 0, false, 0, NULL, NULL, &v, NULL, err, sizeof err);
        if (rc != 1 || !strstr(err, "computes for c <= 14 only") || !strstr(err, "c = 16, c0 = 2, x = 1e+07")) { nfail++; printf("FAIL DK c = 16 c0 = 2 x = 1e7: rc %d, \"%s\"\n", rc, err); }
        else printf("ok   DK c = 16 c0 = 2 x = 1e7: refused, \"%s\"\n", err);
    }
    {   /* the same with verified, where the rule names the conditional integral (c0 < 0) */
        double x = 1e7, v; char err[1024] = ""; dfl_double_opts dv; dfl_double_defaults(&dv); dv.verified = 1; nchk++;
        int rc = dfl_dist(2, 1, &x, 1, 16, -2, 0, false, 0, NULL, &dv, &v, NULL, err, sizeof err);
        if (rc != 1 || !strstr(err, "computes for c <= 14 only")) { nfail++; printf("FAIL DK c = 16 c0 = -2 x = 1e7 verified: rc %d, \"%s\"\n", rc, err); }
        else printf("ok   DK c = 16 c0 = -2 x = 1e7 verified: refused\n");
    }
    {   /* cases DU, DT far above the unit root: the conditional integral at x > 300 (the line integral exhausted the memory at
         * c = 16, c0 = 0, x = 4.4e6), against the large-c law; above c = 14 refused */
        double x1 = 4.07e5, x2 = 4.766e5, x3 = 4.4e6, v; char err[1024] = "";
        nchk++;
        int rc1 = dfl_dist(2, 1, &x1, 2, 13, 0, 0, false, 0, NULL, NULL, &v, NULL, err, sizeof err);
        if (rc1 || !(fabs(v - 0.9545202704) < 1e-6)) { nfail++; printf("FAIL DU c = 13 x = 4.07e5: rc %d, %.12g\n", rc1, v); } else printf("ok   DU c = 13 x = 4.07e5: %.12g\n", v);
        nchk++;
        rc1 = dfl_dist(2, 1, &x2, 3, 13, 0, 0, false, 0, NULL, NULL, &v, NULL, err, sizeof err);
        if (rc1 || !(fabs(v - 0.9973023556) < 1e-6)) { nfail++; printf("FAIL DT c = 13 x = 4.766e5: rc %d, %.12g\n", rc1, v); } else printf("ok   DT c = 13 x = 4.766e5: %.12g\n", v);
        nchk++;
        rc1 = dfl_dist(2, 1, &x3, 3, 16, 0, 0, false, 0, NULL, NULL, &v, NULL, err, sizeof err);
        if (rc1 != 1 || !strstr(err, "computes for c <= 14 only")) { nfail++; printf("FAIL DT c = 16 x = 4.4e6: rc %d, \"%s\"\n", rc1, err); } else printf("ok   DT c = 16 x = 4.4e6: refused\n");
    }
    {   /* the general engine (cases DK, c = 10, c0 = 9: the line cut to the rate with c0; was 0.6372, 0.99735, 1.0000000112) */
        dfl_design D; char err[1024] = ""; double xs[3] = {1e3, 3e3, 1e4}, v[3], r[3] = {0.636120922790, 0.993550212928, 1.0};
        dfl_det_opts o; dfl_det_defaults(&o); o.c0 = 9;
        int rc = dfl_design_trend(-1, &D, err, sizeof err) || dfl_det_dist(2, &D, 3, xs, 10, &o, v, err, sizeof err);
        double e = 0; for (int k = 0; k < 3; k++) e = fmax(e, fabs(v[k] - r[k]));
        nchk++;
        if (rc || !(e < 1e-10)) { nfail++; printf("FAIL the general engine, DK c = 10 c0 = 9: rc %d (%s), worst difference %.1e\n", rc, err, e); }
        else printf("ok   the general engine, DK c = 10 c0 = 9, x = 1e3, 3e3, 1e4: worst difference %.1e\n", e);
        dfl_design_free(&D);
    }
    {   /* the final check of a distribution function */
        double x[2] = {0.5, 1.0}, F1[2] = {0.3, 1 + 5e-9}, F2[2] = {0.3, 1 + 2e-8}; char err[1024] = "";
        nchk++;
        if (dfl_cdf_range_check("ltratio_cdf", "case DK, c = 1", "x", 2, x, F1, err, sizeof err) != 0 ||
            dfl_cdf_range_check("ltratio_cdf", "case DK, c = 1", "x", 2, x, F2, err, sizeof err) != 1 ||
            !strstr(err, "ltratio_cdf: case DK, c = 1, x = 1: the computed distribution function is 1.00000002")) { nfail++; printf("FAIL dfl_cdf_range_check: \"%s\"\n", err); }
        else printf("ok   dfl_cdf_range_check: \"%s\"\n", err);
    }
    printf("%d checks, %d failed, worst absolute difference %.1e, %.1f s (%d threads)\n", nchk, nfail, worst, dfl_wtime() - t0, dfl_nworkers());
    return nfail != 0;
}
