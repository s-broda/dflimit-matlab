/* dfl_test_det.c: the general engine (dfl_det_*.c) against the Julia engine DFLimit.DFDet (reference values of
 * test/det_ref.jl in test/det_ref.txt, first argument): the designs and their preparation, the determinant at complex
 * points, the reciprocal roots, t_*, the outer paths, and the user-level functions (distribution functions, densities,
 * quantiles, the coefficient statistic, moments, the linear form, the law of V, ADFmax, the power envelope, local power,
 * the companions, the router).  Each value line carries its tolerance; the program reports the number of comparisons
 * and the largest deviation of each group, and fails if any value is outside its tolerance. */
#include "dflimit.h"
#include "dfl_det_int.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif
#ifdef __APPLE__
#include <libproc.h>
#include <unistd.h>
#endif

#define MAXD 64
static dfl_design DS[MAXD];
static int HAVE[MAXD];

static void f_sincos(double r, double *o, void *c) { (void) c; o[0] = sin(2 * M_PI * r); o[1] = cos(2 * M_PI * r); }
static void f_perronC04(double r, double *o, void *c) { (void) c; o[0] = 1; o[1] = r; o[2] = r > 0.4 ? 1.0 : 0.0; o[3] = r > 0.4 ? r - 0.4 : 0.0; }
static void f_expr(double r, double *o, void *c) { (void) c; o[0] = exp(r); }

/* the groups of comparisons: name, count, largest deviation (absolute), failures */
typedef struct { const char *name; int n, fail; double worst, worst_rel_tol; } group_t;
static group_t G[64]; static int NG = 0;
static char GNAMES[64][48];
static group_t *grp(const char *name)
{
    for (int i = 0; i < NG; i++) if (!strcmp(G[i].name, name)) return &G[i];
    snprintf(GNAMES[NG], sizeof GNAMES[NG], "%s", name); G[NG].name = GNAMES[NG];
    return &G[NG++];
}
static int VERB = 0;
static void check(const char *gname, double got, double want, double tol, const char *what)
{
    group_t *g = grp(gname); g->n++;
    double d = fabs(got - want);
    if (isnan(want) && isnan(got)) d = 0;
    if (isinf(want) && got == want) d = 0;
    if (!(d <= tol)) { g->fail++; printf("FAIL %-10s %s: %.17g vs %.17g (diff %.2e, tol %.1e)\n", gname, what, got, want, d, tol); }
    else if (VERB) printf("ok   %-10s %s: %.17g (diff %.2e)\n", gname, what, got, d);
    if (d > g->worst) g->worst = d;
    if (tol > 0 && d / tol > g->worst_rel_tol) g->worst_rel_tol = d / tol;
}

static int read_design(char *rest, int id)
{
    char err[512]; char kind[32]; int nc = 0;
    if (sscanf(rest, "%31s%n", kind, &nc) != 1) return 1;
    char *p = rest + nc;
    dfl_design *D = &DS[id]; int rc = 1;
    if (!strcmp(kind, "trend")) { int q = (int) strtol(p, &p, 10); rc = dfl_design_trend(q, D, err, sizeof err); }
    else if (!strcmp(kind, "ltrend")) { int q = (int) strtol(p, &p, 10); rc = dfl_design_ltrend(q, D, err, sizeof err); }
    else if (!strcmp(kind, "perron")) {
        char m[4]; int k; sscanf(p, "%3s%n", m, &k); p += k;
        int q = (int) strtol(p, &p, 10); double lam[16];
        for (int i = 0; i < q; i++) lam[i] = strtod(p, &p);
        rc = dfl_design_perron(m[0], q, lam, D, err, sizeof err);
    } else if (!strcmp(kind, "fdesign")) {
        char nm[32]; sscanf(p, "%31s", nm);
        if (!strcmp(nm, "sincos")) rc = dfl_design_fdesign(f_sincos, NULL, 2, 0, NULL, 0, 0, 0, D, err, sizeof err);
        else if (!strcmp(nm, "perronC04")) { double b = 0.4; rc = dfl_design_fdesign(f_perronC04, NULL, 4, 1, &b, 0, 0, 0, D, err, sizeof err); }
        else if (!strcmp(nm, "expr")) rc = dfl_design_fdesign(f_expr, NULL, 1, 0, NULL, 0, 0, 0, D, err, sizeof err);
    } else if (!strcmp(kind, "matrix")) {
        int base = (int) strtol(p, &p, 10), n = (int) strtol(p, &p, 10);
        double *X = malloc((size_t) n * DS[base].K * sizeof(double));
        dfl_design_regmat(&DS[base], n, X);
        rc = dfl_design_from_matrix(n, DS[base].K, X, -1, NULL, 0, 0, 0, D, err, sizeof err);
        free(X);
    } else if (!strcmp(kind, "first")) { int base = (int) strtol(p, &p, 10); rc = dfl_design_detrend_first(&DS[base], D, err, sizeof err); }
    else if (!strcmp(kind, "raw")) {
        int np = (int) strtol(p, &p, 10), deg = (int) strtol(p, &p, 10), K = (int) strtol(p, &p, 10), loc = (int) strtol(p, &p, 10);
        D->K = K; D->np = np; D->deg = deg; D->loc = loc; D->first = 0;
        D->br = malloc((np + 1) * sizeof(double)); D->C = malloc((size_t) np * (deg + 1) * K * sizeof(double));
        for (int i = 0; i <= np; i++) D->br[i] = strtod(p, &p);
        for (int i = 0; i < np * (deg + 1) * K; i++) D->C[i] = strtod(p, &p);
        rc = 0;
    }
    if (rc) { printf("FAIL design %d (%s): %s\n", id, rest, err); return 1; }
    HAVE[id] = 1;
    return 0;
}
/* the next number of the line (NaN, Inf as Julia prints them) */
static double num(char **p)
{
    while (**p == ' ') (*p)++;
    if (!strncmp(*p, "NaN", 3)) { *p += 3; return NAN; }
    if (!strncmp(*p, "Inf", 3)) { *p += 3; return INFINITY; }
    if (!strncmp(*p, "-Inf", 4)) { *p += 4; return -INFINITY; }
    return strtod(*p, p);
}
static int inum(char **p) { return (int) num(p); }

/* the cells in which the C rule was moved from the Julia rule by the timings of the C engines (dfl_autoengine; the
 * measurements in test/engine_table_c.txt): there the two rules may differ, elsewhere they must agree */
static int moved(int f, int model, int stat, double g0, double gls, double c, double c0, double x)
{
    int fixed = g0 != 0 && !stat, DU = model == 2 && !stat, DT = model == 3 && !stat, DG = model == 4;
    (void) gls;
    switch (f) {
    case DFL_F_DIST: return (!fixed && c0 == 0 && !DG && c < -20) || (!fixed && !stat && c0 < 0 && ((DT && c >= -30 && c <= 0) || (DU && c >= -10 && c <= 0)));
    case DFL_F_ICDF: return c0 == 0 && !fixed && !DG && ((!stat && c <= 2) || (stat && c <= -1));
    case DFL_F_MOM: return c0 != 0;
    case DFL_F_LINCDF: return DG && x != 0;
    case DFL_F_ENVELOPE: return 1;
    default: return 0;
    }
}
/* the user-level value functions of the 'val' and 'err' lines: rc and the values in got */
static int val_call(const char *fn, int id, double cb, int st, double g0, double c, double c0, int ver, double par, int n, const double *x, double *got, char *err)
{
    dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0; o.start = st; o.gamma0 = g0; o.gls = cb; o.verified = ver;
    int rc = 1;
    if (!strcmp(fn, "cdf")) rc = dfl_det_dist(2, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "pdf")) rc = dfl_det_dist(1, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "quant")) rc = dfl_det_quantile(&DS[id], n, x, c, par, 1e-9, &o, got, err, 512);
    else if (!strcmp(fn, "coefcdf")) rc = dfl_det_coef(2, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "coefpdf")) rc = dfl_det_coef(1, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "coefq")) rc = dfl_det_coef_quantile(&DS[id], n, x, c, par, 1e-9, &o, got, err, 512);
    else if (!strcmp(fn, "coefgp")) rc = dfl_det_coef_cdf_gp(&DS[id], n, x, c, &o, 60000, 4000, got, err, 512);
    else if (!strcmp(fn, "mom")) rc = dfl_det_moments(&DS[id], c, n, &o, got, NULL, err, 512);
    else if (!strcmp(fn, "lpow")) { o.c0 = 0; rc = dfl_det_local_power(&DS[id], n, x, 0.05, (int) par, &o, got, err, 512); }
    else if (!strcmp(fn, "lincdf")) rc = dfl_det_lin(2, &DS[id], n, x, par, c, &o, got, err, 512);
    else if (!strcmp(fn, "linpdf")) rc = dfl_det_lin(1, &DS[id], n, x, par, c, &o, got, err, 512);
    else if (!strcmp(fn, "linq")) rc = dfl_det_lin_quantile(&DS[id], x[0], par, c, 0.0, 1e-10, &o, got, err, 512);
    else if (!strcmp(fn, "adfcdf")) rc = dfl_det_adfmax(2, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "adfpdf")) rc = dfl_det_adfmax(1, &DS[id], n, x, c, &o, got, err, 512);
    else if (!strcmp(fn, "adfq")) rc = dfl_det_adfmax_quantile(&DS[id], x[0], c, par, 1e-10, &o, got, err, 512);
    else if (!strcmp(fn, "env")) rc = dfl_det_power_envelope(&DS[id], n, x, par, got, err, 512);
    else if (!strcmp(fn, "ptcdf")) rc = dfl_det_pt(2, &DS[id], n, x, par, c, 0, 0, got, err, 512);
    else if (!strcmp(fn, "ptcdfN")) rc = dfl_det_pt(2, &DS[id], n, x, par, c, (int) c0, 0, got, err, 512);   /* (c0 carries N) */
    else if (!strcmp(fn, "ptpdf")) rc = dfl_det_pt(1, &DS[id], n, x, par, c, 0, 0, got, err, 512);
    else if (!strcmp(fn, "ptq")) rc = dfl_det_pt_quantile(&DS[id], x[0], par, 0, 0, got, err, 512);
    else if (!strcmp(fn, "ptpow")) rc = dfl_det_pt_power(&DS[id], par, c, x[0], 0, 0, &got[0], &got[1], err, 512);
    else if (!strcmp(fn, "ptenv")) rc = dfl_det_pt_envelope(&DS[id], n, x, par, 0, 0, got, err, 512);
    else if (!strcmp(fn, "lmcdf")) rc = dfl_det_lm_cdf(&DS[id], n, x, c, got, err, 512);
    else if (!strcmp(fn, "lmq")) rc = dfl_det_lm_quantile(&DS[id], x[0], c, got, err, 512);
    else if (!strcmp(fn, "lmpow")) rc = dfl_det_lm_power(&DS[id], n, x, par, &got[0], got + 1, err, 512);
    else if (!strcmp(fn, "mzacdf")) rc = dfl_det_mza_cdf(&DS[id], n, x, c, par, NULL, got, err, 512);
    else if (!strcmp(fn, "mzaq")) rc = dfl_det_mza_quantile(&DS[id], x[0], c, par, NAN, 1e-9, NULL, got, err, 512);
    else if (!strcmp(fn, "msbcdf")) rc = dfl_det_msb_cdf(&DS[id], n, x, c, par, got, err, 512);
    else if (!strcmp(fn, "msbq")) rc = dfl_det_msb_quantile(&DS[id], x[0], c, par, got, err, 512);
    else { printf("unknown value function %s\n", fn); snprintf(err, 512, "unknown value function %s", fn); }
    return rc;
}
/* the router's functions of the 'route' and 'rerr' lines */
static int route_call(const char *fn, int e, const char *sp, int a, int st, double g0, double gl, double c, double c0, double par, int n, const double *x, double *got, char *err)
{
    dfl_spec S; S.D = strcmp(sp, "design") ? NULL : &DS[a]; S.model = S.D ? 0 : a; S.start = st; S.g0 = g0; S.gls = gl;   /* (NaN on the lines of a case, except where a line tests the refusal of gls) */
    int rc = 1;
    if (!strcmp(fn, "cdf")) rc = dfl_engine_dist(2, n, x, &S, c, c0, e, 0, 0, got, err, 512);
    else if (!strcmp(fn, "pdf")) rc = dfl_engine_dist(1, n, x, &S, c, c0, e, 0, 0, got, err, 512);
    else if (!strcmp(fn, "icdf")) rc = dfl_engine_icdf(n, x, &S, c, c0, e, 0, 1e-10, got, err, 512);
    else if (!strcmp(fn, "mom")) rc = dfl_engine_mom((int) x[0], &S, c, c0, e, got, err, 512);
    else if (!strcmp(fn, "coef")) rc = dfl_engine_coef(n, x, &S, c, c0, e, got, err, 512);
    else if (!strcmp(fn, "coefq")) rc = dfl_engine_coef_icdf(n, x, &S, c, c0, e, 1e-9, got, err, 512);
    else if (!strcmp(fn, "lin")) rc = dfl_engine_lincdf(n, x, par, a, c, gl, e, got, err, 512);
    else if (!strcmp(fn, "env")) rc = dfl_engine_envelope(n, x, a, par, e, got, err, 512);
    else if (!strcmp(fn, "adf")) rc = dfl_engine_adfmax(2, n, x, &S, c, e, got, err, 512);
    else if (!strcmp(fn, "adfpdf")) rc = dfl_engine_adfmax(1, n, x, &S, c, e, got, err, 512);
    else if (!strcmp(fn, "region1")) { double pr[2] = {par, x[0]}; rc = dfl_engine_region(1, pr, &S, c, e, got, err, 512); }
    else if (!strcmp(fn, "region2")) { double pr[2] = {x[0], 0.0}; rc = dfl_engine_region(2, pr, &S, c, e, got, err, 512); }
    else snprintf(err, 512, "unknown route function %s", fn);
    return rc;
}
/* the warnings of the router and the general engine: cleared at the start of every router call and by
 * dfl_design_from_matrix, the public clear functions, the hold released where a search fails (no reference values) */
static void selftest_warnings(void)
{
    char err[512]; dfl_design T0, M; double v[4], x = -2.0;
    dfl_design_trend(0, &T0, err, sizeof err);
    int n = 50; double X[150];
    for (int t = 0; t < n; t++) { X[t] = 1; X[n + t] = (t + 1.0) / n; X[2 * n + t] = 2 * X[n + t]; }
#define WCHK(cond, what) check("warnings", (cond) ? 0 : 1, 0, 0, what)
    dfl_design_from_matrix(n, 3, X, 0, NULL, 0, 0, 0, &M, err, sizeof err); dfl_design_free(&M);
    WCHK(dfl_det_warning() && strstr(dfl_det_warning(), "linearly dependent"), "design_from_matrix warns of a dependent column");
    dfl_spec S = {NULL, 2, 0, 0.0, NAN};
    dfl_engine_dist(2, 1, &x, &S, -5, 0, DFL_ENGINE_CLASSICAL, 0, 0, v, err, sizeof err);
    WCHK(dfl_det_warning() == NULL, "a classical router call clears dfl_det_warning");
    dfl_warning_join("a stale classical warning");
    dfl_engine_dist(2, 1, &x, &S, -5, 0, DFL_ENGINE_GENERAL, 0, 0, v, err, sizeof err);
    WCHK(dfl_warning() == NULL, "a general router call clears dfl_warning");
    dfl_warning_join("a stale classical warning");
    dfl_design_from_matrix(n, 3, X, 0, NULL, 0, 0, 0, &M, err, sizeof err); dfl_design_free(&M);
    WCHK(dfl_warning() == NULL, "design_from_matrix clears dfl_warning");
    dfl_warning_join("x"); dfl_warning_clear(); WCHK(dfl_warning() == NULL, "dfl_warning_clear");
    det_warn("x"); dfl_det_warning_clear(); WCHK(dfl_det_warning() == NULL, "dfl_det_warning_clear");
    /* the hold of the routed search, released at p outside (0, 1): the next general call clears the warning */
    double p = 1.5;
    int rc = dfl_engine_icdf(1, &p, &S, -5, -1.0, DFL_ENGINE_AUTO, 0, 1e-10, v, err, sizeof err);
    WCHK(rc != 0, "the routed quantile search refuses p = 1.5");
    dfl_design_from_matrix(n, 3, X, 0, NULL, 0, 0, 0, &M, err, sizeof err); dfl_design_free(&M);
    dfl_det_opts o; dfl_det_defaults(&o);
    dfl_det_dist(2, &T0, 1, &x, -5, &o, v, err, sizeof err);
    WCHK(dfl_det_warning() == NULL, "after a refused routed search the general engine's warning is cleared again (the hold released)");
    /* a hold left by a call is dropped by dfl_det_warning_clear */
    det_warn_hold(1); det_warn("held"); dfl_det_warning_clear();
    dfl_det_dist(2, &T0, 1, &x, -5, &o, v, err, sizeof err);
    det_warn("y"); dfl_det_dist(2, &T0, 1, &x, -5, &o, v, err, sizeof err);
    WCHK(dfl_det_warning() == NULL, "dfl_det_warning_clear drops a hold");
#undef WCHK
    dfl_design_free(&T0);
}
/* the memory of the process in MB: the footprint on macOS (what Activity Monitor shows), the peak resident size on other
 * unix systems, 0 on Windows (not measured there) */
static double mem_mb(void)
{
#ifdef _WIN32
    return 0;
#else
#ifdef __APPLE__
    struct rusage_info_v0 ri;
    if (proc_pid_rusage(getpid(), RUSAGE_INFO_V0, (rusage_info_t *) &ri) == 0) return ri.ri_phys_footprint / 1048576.0;
#endif
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss / 1024.0;          /* (kilobytes on Linux) */
#endif
}
/* the guard on the nodes of the U-shapes, before the outer sum (ushape_guard): the design (1 + r, 1{r <= 0.4}) on the
 * vertical line, whose positive root runs off, and d(r) = 1 + r with a fixed start at c = 0 are refused with the error of
 * the guard within half a second and with a few MB (the guard in the U-shapes refused after 1.4 s and 1.2 GB); run first,
 * so that the memory measured is that of these calls */
static void selftest_guard(void)
{
    char err[512]; double v, x = -2.0;
    double br2[3] = {0, 0.4, 1}, C2[8] = {1, 1, 1, 0, 1, 1, 0, 0}, br1[2] = {0, 1}, C1[2] = {1, 1};
    dfl_design D2 = {2, 2, 1, br2, C2, 0, 0}, D1 = {1, 1, 1, br1, C1, 0, 0};
    dfl_det_opts o; dfl_det_defaults(&o);
    double m0 = mem_mb(), t0 = dfl_wtime();
    int rc = dfl_det_dist(2, &D2, 1, &x, 1.0, &o, &v, err, sizeof err);
    double dt = dfl_wtime() - t0, dm = mem_mb() - m0;
    if (VERB) printf("     guard, zero start: %.3f s, %+.0f MB: %s\n", dt, dm, err);
    check("guard", rc && strstr(err, "the U-shapes of the zero start") && strstr(err, "a root runs off") ? 0 : 1, 0, 0, "(1 + r, 1{r <= 0.4}) at c = 1 refused by the guard");
    check("guard", dt < 0.5 ? 0 : dt, 0, 0, "(1 + r, 1{r <= 0.4}) at c = 1 refused within 0.5 s (the deviation: the time)");
    check("guard", dm < 100 ? 0 : dm, 0, 0, "(1 + r, 1{r <= 0.4}) at c = 1 refused with less than 100 MB (the deviation: MB)");
    o.start = DFL_START_FIXED; o.gamma0 = 1.0;
    m0 = mem_mb(); t0 = dfl_wtime();
    rc = dfl_det_dist(2, &D1, 1, &x, 0.0, &o, &v, err, sizeof err);
    dt = dfl_wtime() - t0; dm = mem_mb() - m0;
    if (VERB) printf("     guard, fixed start: %.3f s, %+.0f MB: %s\n", dt, dm, err);
    check("guard", rc && strstr(err, "the U-shapes of the fixed start") && strstr(err, "the roots run off") ? 0 : 1, 0, 0, "d(r) = 1 + r, g0 = 1 at c = 0 refused by the guard");
    check("guard", dt < 0.5 ? 0 : dt, 0, 0, "d(r) = 1 + r, g0 = 1 at c = 0 refused within 0.5 s (the deviation: the time)");
    check("guard", dm < 100 ? 0 : dm, 0, 0, "d(r) = 1 + r, g0 = 1 at c = 0 refused with less than 100 MB (the deviation: MB)");
    /* the caps of the outer hyperbola (a factor hfac of its step far below one) and of the grid of the moments */
    dfl_design T1; dfl_design_trend(0, &T1, err, sizeof err);
    dfl_det_defaults(&o); o.hfac = 1e-4;
    rc = dfl_det_dist(2, &T1, 1, &x, -5.0, &o, &v, err, sizeof err);
    check("guard", rc && strstr(err, "the outer contour would need") ? 0 : 1, 0, 0, "hfac = 1e-4 refused (the outer contour beyond 1e5 nodes)");
    dfl_det_defaults(&o); double m[2];
    rc = dfl_det_moments(&T1, 1e5, 2, &o, m, NULL, err, sizeof err);
    check("guard", rc && strstr(err, "beyond the range of the moments") ? 0 : 1, 0, 0, "the moments at c = 1e5 refused (the grid beyond 1e6 points)");
    dfl_design_free(&T1);
    dfl_det_empty_caches();
}
/* the public functions at non-finite arguments (November 2026): a limit where there is one, else an error that names the
 * argument; none of them may hang (coef_cdf_gp at y = +-Inf, the searches from x0 = +-Inf, the point-optimal statistic at
 * c = +-Inf and cbar = -Inf did not return) */
static void selftest_nonfinite(void)
{
    char err[512]; double v[4], w[4]; int rc;
    double INF = INFINITY, a2[2];
    dfl_design DU, DT, PA; dfl_design_trend(0, &DU, err, sizeof err); dfl_design_trend(1, &DT, err, sizeof err);
    double lam = 0.5; dfl_design_perron('A', 1, &lam, &PA, err, sizeof err);
    dfl_det_opts o; dfl_det_defaults(&o);
#define NFV(cond, what) check("nonfinite", (cond) ? 0 : 1, 0, 0, what)
#define NFE(call, txt, what) do { err[0] = 0; rc = (call); NFV(rc && strstr(err, txt), what); } while (0)
    a2[0] = INF; a2[1] = -INF;
    rc = dfl_det_coef_cdf_gp(&DU, 2, a2, -5, &o, 0, 0, v, err, sizeof err);
    NFV(!rc && v[0] == 1 && v[1] == 0, "coef_cdf_gp at y = Inf, -Inf: 1, 0");
    a2[0] = NAN; rc = dfl_det_coef_cdf_gp(&DU, 1, a2, -5, &o, 0, 0, v, err, sizeof err);
    NFV(!rc && isnan(v[0]), "coef_cdf_gp at y = NaN: NaN");
    double p = 0.05, x = -2, b = -0.26;
    NFE(dfl_det_coef_cdf_gp(&DU, 1, &x, -5, &o, 0, INF, v, err, sizeof err), "umax must be", "coef_cdf_gp refuses umax = Inf");
    NFE(dfl_det_quantile(&DU, 1, &p, -5, INF, 1e-9, &o, v, err, sizeof err), "x0 must be a finite", "dist_quantile refuses x0 = Inf");
    NFE(dfl_det_quantile(&DU, 1, &p, -5, -INF, 1e-9, &o, v, err, sizeof err), "x0 must be a finite", "dist_quantile refuses x0 = -Inf");
    NFE(dfl_det_quantile(&DU, 1, &p, -5, -2, INF, &o, v, err, sizeof err), "tol must be", "dist_quantile refuses tol = Inf");
    NFE(dfl_det_coef_quantile(&DU, 1, &p, -5, INF, 1e-9, &o, v, err, sizeof err), "y0 must be a finite", "coef_quantile refuses y0 = Inf");
    NFE(dfl_det_lin_quantile(&DU, p, -3.5, -5, -INF, 1e-10, &o, v, err, sizeof err), "b0 must be a finite", "lin_quantile refuses b0 = -Inf");
    NFE(dfl_det_adfmax_quantile(&DU, p, -5, INF, 1e-10, &o, v, err, sizeof err), "x0 must be a finite", "adfmax_quantile refuses x0 = Inf");
    NFE(dfl_det_mza_quantile(&DT, p, -5, -13.5, -INF, 1e-9, NULL, v, err, sizeof err), "y0 must be a finite", "mza_quantile refuses y0 = -Inf");
    NFE(dfl_det_lin(2, &DU, 1, &b, INF, -5, &o, v, err, sizeof err), "a must be a finite", "lin_cdf refuses a = Inf");
    NFE(dfl_det_lin(1, &DU, 1, &b, -INF, -5, &o, v, err, sizeof err), "a must be a finite", "lin_pdf refuses a = -Inf");
    double pl[2] = {NAN, -0.26};
    NFE(dfl_det_region(&DU, 1, pl, -5, &o, v, err, sizeof err), "slope a of the line", "region_prob refuses the slope a = NaN");
    a2[0] = -INF; a2[1] = NAN; double cs3[3] = {-INF, NAN, INF};
    rc = dfl_det_power_envelope(&DU, 3, cs3, 0.05, v, err, sizeof err);
    NFV(!rc && v[0] == 1 && isnan(v[1]) && v[2] == 0.05, "power_envelope at c = -Inf, NaN, Inf: 1, NaN, the level");
    rc = dfl_det_pt_envelope(&DT, 3, cs3, 0.05, 0, 0, v, err, sizeof err);
    NFV(!rc && v[0] == 1 && isnan(v[1]) && v[2] == 0.05, "pt_envelope at c = -Inf, NaN, Inf: 1, NaN, the level");
    double q = 4;
    NFE(dfl_det_pt(2, &DT, 1, &q, -13.5, INF, 0, 0, v, err, sizeof err), "c must be a finite", "pt_cdf refuses c = Inf");
    NFE(dfl_det_pt(2, &DT, 1, &q, -INF, -5, 0, 0, v, err, sizeof err), "needs cbar < 0", "pt_cdf refuses cbar = -Inf");
    NFE(dfl_det_pt_power(&DT, -13.5, -INF, 0.05, 0, 0, v, w, err, sizeof err), "c must be a finite", "pt_power refuses c = -Inf");
    NFE(dfl_det_pt(2, &DT, 1, &q, -13.5, -5, 0, INF, v, err, sizeof err), "umax must be", "pt_cdf refuses umax = Inf");
    double xm = 0.2;
    NFE(dfl_det_msb_cdf(&DT, 1, &xm, NAN, -13.5, v, err, sizeof err), "c must be a finite", "msb_cdf refuses c = NaN");
    NFE(dfl_det_msb_quantile(&DT, p, INF, -13.5, v, err, sizeof err), "c must be a finite", "msb_quantile refuses c = Inf");
    o.gls = -INF;
    NFE(dfl_det_dist(2, &DT, 1, &x, -5, &o, v, err, sizeof err), "needs cbar < 0", "dist_cdf refuses gls = -Inf");
    dfl_det_defaults(&o);
    double cs2[2] = {-5, -5};
    NFE(dfl_det_local_power(&DU, 2, cs2, NAN, 0, &o, v, err, sizeof err), "local_power: the level must lie in (0, 1)", "local_power checks its level");
    o.c0 = INF;
    NFE(dfl_det_local_power(&DU, 2, cs2, 0.05, 0, &o, v, err, sizeof err), "c0 must be a finite", "local_power refuses c0 = Inf");
    dfl_det_defaults(&o);
    /* the router */
    NFE(dfl_engine_lincdf(1, &b, INF, 1, -5, NAN, DFL_ENGINE_AUTO, v, err, sizeof err), "a must be a finite", "router: the linear form refuses a = Inf");
    NFE(dfl_engine_lincdf(1, &b, -3.5, 4, -5, -INF, DFL_ENGINE_AUTO, v, err, sizeof err), "cbar must be a finite", "router: the linear form of case DG refuses cbar = -Inf");
    rc = dfl_engine_envelope(3, cs3, 1, 0.05, DFL_ENGINE_AUTO, v, err, sizeof err);
    NFV(!rc && v[0] == 1 && isnan(v[1]) && v[2] == 0.05, "router: the envelope at c = -Inf, NaN, Inf: 1, NaN, the level");
    NFE(dfl_engine_envelope(1, cs2, 1, NAN, DFL_ENGINE_AUTO, v, err, sizeof err), "power_envelope: the level must lie in (0, 1)", "router: the envelope's level message");
    dfl_spec SK = {NULL, 1, 0, 0.0, NAN}, SA = {&PA, 0, 0, 0.0, NAN};
    for (int e = 0; e < 2; e++) {
        const dfl_spec *S = e ? &SA : &SK; const char *nm = e ? "Perron A (general)" : "case DK (classical)";
        char what[160];
        double pr[2] = {-3.5, INF}; rc = dfl_engine_region(1, pr, S, -5, DFL_ENGINE_AUTO, v, err, sizeof err);
        snprintf(what, sizeof what, "router: the region below the line at b = Inf is 1, %s", nm); NFV(!rc && v[0] == 1, what);
        pr[1] = -INF; rc = dfl_engine_region(1, pr, S, -5, DFL_ENGINE_AUTO, v, err, sizeof err);
        snprintf(what, sizeof what, "router: the region below the line at b = -Inf is 0, %s", nm); NFV(!rc && v[0] == 0, what);
        pr[1] = NAN; rc = dfl_engine_region(1, pr, S, -5, DFL_ENGINE_AUTO, v, err, sizeof err);
        snprintf(what, sizeof what, "router: the region below the line at b = NaN is NaN, %s", nm); NFV(!rc && isnan(v[0]), what);
        pr[0] = INF; pr[1] = -0.26;
        snprintf(what, sizeof what, "router: the region below the line refuses a = Inf, %s", nm);
        NFE(dfl_engine_region(1, pr, S, -5, DFL_ENGINE_AUTO, v, err, sizeof err), "slope a of the line", what);
    }
#undef NFE
#undef NFV
    dfl_design_free(&DU); dfl_design_free(&DT); dfl_design_free(&PA);
    dfl_det_empty_caches();
}
static det_prep *prep_of(int id, double cb, char *err)
{
    det_prep *P = NULL;
    if (det_prep_make(&DS[id], cb, &P, err, 512)) return NULL;
    return P;
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/det_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    VERB = getenv("DFL_VERBOSE") != NULL;
    static char line[1 << 16]; char err[512], what[256];

    double t0 = dfl_wtime();
    selftest_guard();
    while (fgets(line, sizeof line, f)) {

        char tag[32]; int nc = 0;
        if (sscanf(line, "%31s%n", tag, &nc) != 1) continue;
        char *p = line + nc;
        if (!strcmp(tag, "design")) { int id = inum(&p); read_design(p, id); continue; }
        if (!strcmp(tag, "prep")) {
            int id = inum(&p); double cb = num(&p);
            int N = inum(&p), pp = inum(&p), nj = inum(&p), r = inum(&p), rp = inum(&p), rm = inum(&p);
            double Dv[64]; for (int k = 0; k < r; k++) Dv[k] = num(&p);
            double tol = num(&p);
            det_prep *P = prep_of(id, cb, err);
            if (!P) { grp("prep")->fail++; printf("FAIL prep %d: %s\n", id, err); continue; }
            snprintf(what, sizeof what, "design %d cbar %g: N p nj r rp rm", id, cb);
            check("prep", P->N * 1e10 + P->p * 1e8 + P->nj * 1e6 + P->r * 1e4 + P->rp * 1e2 + P->rm, N * 1e10 + pp * 1e8 + nj * 1e6 + r * 1e4 + rp * 1e2 + rm, 0, what);
            if (P->r == r) for (int k = 0; k < r; k++) { snprintf(what, sizeof what, "design %d cbar %g: Dv[%d]", id, cb, k); check("prep", P->Dv[k], Dv[k], tol * fmax(1, fabs(Dv[k])), what); }
            det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "model")) {
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p);
            int N = inum(&p), r = inum(&p), rp = inum(&p), rm = inum(&p); double a = num(&p), ka = num(&p), tol = num(&p);
            det_prep *P = prep_of(id, cb, err); det_model Mo;
            if (!P || det_model_make(P, st, g0, &Mo, err, sizeof err)) { grp("model")->fail++; printf("FAIL model %d: %s\n", id, err); continue; }
            snprintf(what, sizeof what, "design %d start %d: N r rp rm a kappa", id, st);
            check("model", Mo.N * 1e6 + Mo.r * 1e4 + Mo.rp * 1e2 + Mo.rm, N * 1e6 + r * 1e4 + rp * 1e2 + rm, 0, what);
            check("model", Mo.a, a, tol, what); check("model", Mo.kappa, ka, tol, what);
            det_model_free(&Mo); det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "dval")) {
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p), c = num(&p);
            double sr = num(&p), si = num(&p), tr = num(&p), ti = num(&p), dr = num(&p), di = num(&p), tol = num(&p);
            det_prep *P = prep_of(id, cb, err); det_model Mo;
            if (!P || det_model_make(P, st, g0, &Mo, err, sizeof err)) { grp("dval")->fail++; printf("FAIL dval %d: %s\n", id, err); continue; }
            det_work W; char buf[DET_WORK_STACK]; det_work_init(&W, &Mo, buf, sizeof buf);
            cdouble d = det_dval_st(&W, c, CPX(sr, si), CPX(tr, ti));
            snprintf(what, sizeof what, "design %d cbar %g start %d c %g s %g%+gi t %g%+gi", id, cb, st, c, sr, si, tr, ti);
            check("dval", cabs(csub(d, CPX(dr, di))), 0, tol, what);
            det_work_free(&W); det_model_free(&Mo); det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "roots")) {
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p), c = num(&p);
            double tr = num(&p), ti = num(&p), ldr = num(&p), cfr = num(&p), cfi = num(&p); int r = inum(&p);
            cdouble mu[64]; for (int k = 0; k < r; k++) { double a = num(&p), b = num(&p); mu[k] = CPX(a, b); }
            double tol = num(&p);
            det_prep *P = prep_of(id, cb, err); det_model Mo;
            if (!P || det_model_make(P, st, g0, &Mo, err, sizeof err)) { grp("roots")->fail++; printf("FAIL roots %d: %s\n", id, err); continue; }
            det_work W; char buf[DET_WORK_STACK]; det_work_init(&W, &Mo, buf, sizeof buf);
            det_fixt *fx = st == 1 ? det_fixt_alloc(Mo.r) : NULL;
            cdouble th, ld, m[64];
            int rc = det_mgfparts(&W, c, CPX(tr, ti), &th, &ld, m, fx);
            snprintf(what, sizeof what, "design %d cbar %g start %d c %g t %g%+gi: log|D_V|", id, cb, st, c, tr, ti);
            check("roots", rc ? NAN : creal(ld), ldr, tol * fmax(1, fabs(ldr)), what);
            double mx = 0; for (int k = 0; k < r; k++) mx = fmax(mx, cabs(mu[k]));
            int used[64] = {0};
            for (int k = 0; k < r && Mo.r == r; k++) {
                int best = -1; double bd = INFINITY;
                for (int l = 0; l < r; l++) if (!used[l] && cabs(csub(m[l], mu[k])) < bd) { bd = cabs(csub(m[l], mu[k])); best = l; }
                used[best] = 1;
                snprintf(what, sizeof what, "design %d cbar %g start %d c %g t %g%+gi: mu[%d]", id, cb, st, c, tr, ti, k);
                check("roots", bd, 0, tol * mx, what);
            }
            if (fx) { snprintf(what, sizeof what, "design %d start 1 c %g t %g%+gi: cf", id, c, tr, ti); check("roots", cabs(csub(fx->cf, CPX(cfr, cfi))), 0, tol * fmax(1, cabs(CPX(cfr, cfi))), what); free(fx); }
            det_work_free(&W); det_model_free(&Mo); det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "tstar")) {
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p), c = num(&p), ts = num(&p); int mult = inum(&p); double tol = num(&p);
            det_prep *P = prep_of(id, cb, err); det_model Mo;
            if (!P || det_model_make(P, st, g0, &Mo, err, sizeof err)) { grp("tstar")->fail++; printf("FAIL tstar %d: %s\n", id, err); continue; }
            double tsc; int mc;
            double t1 = dfl_wtime();
            int rc = det_tstar(&Mo, c, &tsc, &mc, err, sizeof err);
            if (VERB) printf("     tstar %.2f ms\n", 1e3 * (dfl_wtime() - t1));
            snprintf(what, sizeof what, "design %d cbar %g start %d c %g: t_* (multiplicity %d vs %d)", id, cb, st, c, mc, mult);
            check("tstar", rc ? NAN : tsc, ts, tol, what);
            if (mc != mult) { grp("tstar")->fail++; printf("FAIL tstar multiplicity %d vs %d (%s)\n", mc, mult, what); }
            det_model_free(&Mo); det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "path")) {
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p), c = num(&p);
            double tau = num(&p), ts = num(&p), a = num(&p), b = num(&p), h = num(&p), phi = num(&p); int K = inum(&p); double mnp = num(&p), mxn = num(&p); int nk = inum(&p);
            det_path *CP = det_path_cached(&DS[id], cb, st, g0, c, NAN, 1.0, err, sizeof err);
            snprintf(what, sizeof what, "design %d cbar %g start %d g0 %g c %g", id, cb, st, g0, c);
            if (!CP) { grp("path")->fail++; printf("FAIL path %s: %s\n", what, err); continue; }
            int r = CP->r;
            double vals[64 * 1024]; int nv = 0;
            for (int j = 0; j < nk; j++) { int k = inum(&p); vals[nv++] = k; vals[nv++] = num(&p); vals[nv++] = num(&p); for (int i = 0; i < 2 * r; i++) vals[nv++] = num(&p); }
            double tol = num(&p);
            char w2[300];
            snprintf(w2, sizeof w2, "%s: K (%d vs %d), phi (%g vs %g)", what, CP->K, K, CP->phi, phi);
            check("path", CP->K, K, 0, w2); check("path", CP->phi, phi, 1e-14, w2);
            snprintf(w2, sizeof w2, "%s: tau, t_*, a, b, h", what);
            check("path", CP->tau, tau, tol * fmax(1, tau), w2); check("path", CP->tstar, ts, tol * fmax(1, ts), w2);
            check("path", CP->a, a, tol * fmax(1, fabs(a)), w2); check("path", CP->b, b, tol * fmax(1, fabs(b)), w2);
            if (!(isnan(h) && isnan(CP->h))) check("path", CP->h, h, tol * fmax(1, fabs(h)), w2);
            snprintf(w2, sizeof w2, "%s: minpos, maxneg", what);
            check("path", CP->minpos, mnp, tol * fmax(1, fabs(mnp)), w2); if (isfinite(mxn)) check("path", CP->maxneg, mxn, tol * fmax(1, fabs(mxn)), w2);
            if (CP->K == K) for (int j = 0, q = 0; j < nk; j++) {
                int k = (int) vals[q++]; double lr = vals[q++], li = vals[q++];
                snprintf(w2, sizeof w2, "%s: log P and roots at node %d", what, k);
                double sc = fmax(1, hypot(lr, li));
                check("path", cabs(csub(CP->lP[k], CPX(lr, li))), 0, tol * sc, w2);
                for (int i = 0; i < r; i++) { double rr = vals[q++], ri = vals[q++]; check("path", cabs(csub(CP->R[(size_t) k * r + i], CPX(rr, ri))), 0, tol * fmax(1, hypot(rr, ri)), w2); }
            }
            det_path_release(CP);
            continue;
        }
        if (!strcmp(tag, "val")) {
            char fn[32]; int k2; sscanf(p, "%31s%n", fn, &k2); p += k2;
            int id = inum(&p); double cb = num(&p); int st = inum(&p); double g0 = num(&p), c = num(&p), c0 = num(&p); int ver = inum(&p); double par = num(&p); int n = inum(&p);
            double *x = malloc(n * sizeof(double)), *want = malloc((n + 8) * sizeof(double)), *got = malloc((n + 8) * sizeof(double));
            for (int i = 0; i < n; i++) x[i] = num(&p);
            int nw = inum(&p);
            for (int i = 0; i < nw; i++) want[i] = num(&p);
            double tol = num(&p);
            int brittle = strstr(p, "brittle") != NULL;
            double t1 = dfl_wtime(); int rc = val_call(fn, id, cb, st, g0, c, c0, ver, par, n, x, got, err);
            double dt = dfl_wtime() - t1;
            char gname[48]; snprintf(gname, sizeof gname, "%s%s", fn, brittle ? "*" : "");
            snprintf(what, sizeof what, "%s design %d cbar %g start %d g0 %g c %g c0 %g%s", fn, id, cb, st, g0, c, c0, ver ? " verified" : "");
            if (rc) { grp(gname)->fail++; grp(gname)->n++; printf("FAIL %s: %s\n", what, err); }
            else for (int i = 0; i < nw; i++) {
                char w2[400]; snprintf(w2, sizeof w2, "%s at %g", what, x[i < n ? i : 0]);
                double tl = tol;
                if (!strcmp(fn, "mom")) tl = tol * fmax(1, fabs(want[i]));
                check(gname, got[i], want[i], tl, w2);
            }
            if (VERB) printf("     %.1f ms %s%s%s\n", 1e3 * dt, what, dfl_det_warning() ? " warning: " : "", dfl_det_warning() ? dfl_det_warning() : "");
            free(x); free(want); free(got);
            dfl_det_empty_caches();
            continue;
        }
        if (!strcmp(tag, "err") || !strcmp(tag, "rerr")) {
            /* a call that must fail with an error whose message contains the text after '|' */
            char fn[32], sp[16] = ""; int k2; sscanf(p, "%31s%n", fn, &k2); p += k2;
            int r = !strcmp(tag, "rerr"), e = 0, id, st; double cb = NAN, gl = NAN, g0, c, c0, par;
            if (r) { e = inum(&p); sscanf(p, "%15s%n", sp, &k2); p += k2; id = inum(&p); st = inum(&p); g0 = num(&p); gl = num(&p); c = num(&p); c0 = num(&p); par = num(&p); }
            else { id = inum(&p); cb = num(&p); st = inum(&p); g0 = num(&p); c = num(&p); c0 = num(&p); (void) inum(&p); par = num(&p); }
            int n = inum(&p); double x[64], got[72];
            for (int i = 0; i < n; i++) x[i] = num(&p);
            char *bar = strchr(p, '|'), txt[256] = "";
            if (bar) { snprintf(txt, sizeof txt, "%s", bar + 2); size_t l = strlen(txt); while (l && (txt[l - 1] == '\n' || txt[l - 1] == '\r')) txt[--l] = 0; }
            err[0] = 0;
            double t1 = dfl_wtime();
            int rc = r ? route_call(fn, e, sp, id, st, g0, gl, c, c0, par, n, x, got, err) : val_call(fn, id, cb, st, g0, c, c0, 0, par, n, x, got, err);
            group_t *g = grp(r ? "route-err" : "errors"); g->n++;
            snprintf(what, sizeof what, "%s %s %s %d start %d g0 %g gls/cbar %g c %g c0 %g x[0] %g (%.2f s)", tag, fn, sp, id, st, g0, r ? gl : cb, c, c0, n ? x[0] : NAN, dfl_wtime() - t1);
            if (!rc || !strstr(err, txt)) { g->fail++; printf("FAIL %s: %s (wanted an error with \"%s\")\n", what, rc ? err : "no error", txt); }
            else if (VERB) printf("ok   %s: %s\n", what, err);
            dfl_det_empty_caches();
            continue;
        }
        if (!strcmp(tag, "lmtstar")) {
            /* the edge of the LM statistic's transform: t_* and the sign of D_V on (0.01, 0.99) t_* and at 1.01 t_* */
            int id = inum(&p); double c = num(&p), ts = num(&p), tol = num(&p), a1, tsc;
            det_prep *P; det_model Mo;
            if (det_prep_lm(&DS[id], &P, &a1, err, sizeof err) || det_model_make(P, 0, 0.0, &Mo, err, sizeof err)) { grp("lmtstar")->fail++; printf("FAIL lmtstar %d: %s\n", id, err); continue; }
            int rc = det_tstar(&Mo, c, &tsc, NULL, err, sizeof err);
            snprintf(what, sizeof what, "LM statistic, design %d c %g: t_*", id, c);
            check("lmtstar", rc ? NAN : tsc, ts, tol, what);
            det_work W; char buf[DET_WORK_STACK]; det_work_init(&W, &Mo, buf, sizeof buf);
            int bad = 0; cdouble th, ld, mu[64];
            for (int i = 0; i <= 200 && !rc; i++) {
                double t = i < 200 ? (0.01 + 0.98 * i / 199.0) * tsc : 1.01 * tsc;
                det_mgfparts(&W, c, CPX(t, 0), &th, &ld, mu, NULL);
                double dv = creal(cexp(ld));
                if (i < 200 ? !(dv > 0) : !(dv < 0)) bad = 1;
            }
            snprintf(what, sizeof what, "LM statistic, design %d c %g: D_V > 0 below t_* and < 0 at 1.01 t_*", id, c);
            check("lmtstar", bad, 0, 0, what);
            det_work_free(&W); det_model_free(&Mo); det_prep_free(P);
            continue;
        }
        if (!strcmp(tag, "auto")) {
            int fi = inum(&p), m = inum(&p), st = inum(&p); double g0 = num(&p), gl = num(&p), c = num(&p), c0 = num(&p), x = num(&p); int ej = inum(&p);
            static const int FM[7] = {0, DFL_F_DIST, DFL_F_ICDF, DFL_F_MOM, DFL_F_LINCDF, DFL_F_ADFMAX, DFL_F_ENVELOPE};
            int f2 = FM[fi], ec = dfl_autoengine(f2, m, st == 2, g0, m == 4 ? (isnan(gl) ? DFL_CBAR_GLS : gl) : gl, c, c0, x) == DFL_ENGINE_GENERAL;
            group_t *g = grp(ec == ej ? "rule" : "rule-moved"); g->n++;
            snprintf(what, sizeof what, "rule f %d model %d stat %d g0 %g gls %g c %g c0 %g x %g: C %s, Julia %s", fi, m, st == 2, g0, gl, c, c0, x, ec ? "general" : "classical", ej ? "general" : "classical");
            if (ec != ej && !moved(f2, m, st == 2, g0, gl, c, c0, x)) { g->fail++; printf("FAIL %s\n", what); }
            else if (VERB && ec != ej) printf("moved %s\n", what);
            continue;
        }
        if (!strcmp(tag, "route")) {
            char fn[32], sp[16]; int k2; sscanf(p, "%31s%n", fn, &k2); p += k2;
            int e = inum(&p); sscanf(p, "%15s%n", sp, &k2); p += k2;
            int a = inum(&p), st = inum(&p); double g0 = num(&p), gl = num(&p), c = num(&p), c0 = num(&p), par = num(&p); int n = inum(&p);
            double *x = malloc(n * sizeof(double)), *want = malloc(64 * sizeof(double)), *got = malloc(64 * sizeof(double));
            for (int i = 0; i < n; i++) x[i] = num(&p);
            int nw = inum(&p); for (int i = 0; i < nw; i++) want[i] = num(&p);
            double tol = num(&p);
            int rc = route_call(fn, e, sp, a, st, g0, gl, c, c0, par, n, x, got, err);
            char gname[48]; snprintf(gname, sizeof gname, "route-%s", e == 0 ? "auto" : e == 1 ? "classical" : "general");
            snprintf(what, sizeof what, "%s engine %d %s %d start %d g0 %g gls %g c %g c0 %g", fn, e, sp, a, st, g0, gl, c, c0);
            if (rc) { grp(gname)->fail++; grp(gname)->n++; printf("FAIL %s: %s\n", what, err); }
            else for (int i = 0; i < nw; i++) { char w2[400]; snprintf(w2, sizeof w2, "%s [%d]", what, i); check(gname, got[i], want[i], !strcmp(fn, "mom") ? tol * fmax(1, fabs(want[i])) : tol, w2); }
            free(x); free(want); free(got);
            dfl_det_empty_caches();
            continue;
        }
    }
    fclose(f);
    selftest_warnings();
    selftest_nonfinite();
    int fails = 0, n = 0;
    printf("%-10s %8s %8s %12s %10s\n", "group", "values", "failed", "largest dev", "dev/tol");
    for (int i = 0; i < NG; i++) { printf("%-10s %8d %8d %12.2e %10.2e\n", G[i].name, G[i].n, G[i].fail, G[i].worst, G[i].worst_rel_tol); fails += G[i].fail; n += G[i].n; }
    printf("%d comparisons, %d failed, %.1f s (%d threads)\n", n, fails, dfl_wtime() - t0, dfl_nworkers());
    return fails ? 1 : 0;
}
