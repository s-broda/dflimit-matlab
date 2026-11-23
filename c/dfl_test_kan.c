/* dfl_test_kan.c: the devices of Bao and Kan (dfl_kan.c, dfl_det_kan.c, the router's defaults in dfl_det_route.c) against
 * the values of the Julia package DFLimit.jl (test/kan_ref.txt, written by test/kan_ref.jl): the closed-form transform of the
 * six classical cases at points and in Taylor-mode arithmetic, the scale, the single inversion of the coefficient statistic
 * (distribution function, density, quantiles), of the linear form U - aV and of the power envelope, the moments of the
 * t-ratio and of the coefficient statistic from jets, the closed form at zero; each kind of line with its own tolerance
 * (at least ten times the largest difference measured on macOS, November 2026), and the largest difference per kind is
 * printed.  The program ends with a line that counts the checks and the failures. */
#include "dfl_det_int.h"
#include "dfl_kan_int.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { const char *tag; double tol, worst; int n, fail; } kind_t;
static kind_t KINDS[] = {
    {"lcf", 1e-12, 0, 0, 0}, {"jet", 1e-11, 0, 0, 0}, {"scale", 1e-12, 0, 0, 0}, {"coef", 1e-11, 0, 0, 0}, {"coefpdf", 1e-11, 0, 0, 0},
    {"coefq", 1e-8, 0, 0, 0}, {"lin", 1e-11, 0, 0, 0}, {"env", 1e-11, 0, 0, 0}, {"tmom", 1e-10, 0, 0, 0}, {"cmom", 1e-10, 0, 0, 0},
    {"ray", 1e-13, 0, 0, 0}, {"rcoef", 1e-11, 0, 0, 0}, {"rcoefpdf", 1e-11, 0, 0, 0}, {"rcoefq", 1e-8, 0, 0, 0}, {"rcoefmom", 1e-10, 0, 0, 0},
    {"rmom", 1e-10, 0, 0, 0}, {"rlin", 1e-11, 0, 0, 0}, {"renv", 1e-11, 0, 0, 0}, {"self", 0, 0, 0, 0},
    {"dcf", 1e-10, 0, 0, 0}, {"ptcdf", 1e-11, 0, 0, 0}, {"ptpdf", 1e-11, 0, 0, 0}, {"ptq", 1e-10, 0, 0, 0}, {"ptpow", 1e-10, 0, 0, 0},
    {"ptenv", 1e-11, 0, 0, 0}, {"gp", 1e-11, 0, 0, 0},
};
static kind_t *kind(const char *tag) { for (size_t i = 0; i < sizeof KINDS / sizeof KINDS[0]; i++) if (!strcmp(KINDS[i].tag, tag)) return &KINDS[i]; return NULL; }
static int nfail = 0, ntot = 0;
static void check(const char *tag, double d, const char *what)
{
    kind_t *k = kind(tag); ntot++;
    if (!k) { printf("unknown kind %s\n", tag); nfail++; return; }
    k->n++; if (d > k->worst || isnan(d)) k->worst = isnan(d) ? INFINITY : fmax(k->worst, d);
    if (!(d <= k->tol)) { k->fail++; nfail++; printf("FAIL %s %s: difference %.3g (tolerance %.1g)\n", tag, what, d, k->tol); }
}
static int rd(FILE *f, double *x) { return fscanf(f, "%lf", x) == 1; }
static int rdc(FILE *f, cdouble *z) { double a, b; if (fscanf(f, "%lf %lf", &a, &b) != 2) return 0; *z = CPX(a, b); return 1; }
static double cdist(cdouble a, cdouble b) { return cabs(csub(a, b)); }

/* the lines of the router's functions with their defaults: 0, or 1 at the end of the file */
static int router_line(FILE *f, const char *tag)
{
    char err[512], what[256]; double x[64], v[64], got[64], est[64];
    if (!strcmp(tag, "renv")) {
        double alpha; int n; if (!(rd(f, &alpha) && fscanf(f, "%d", &n) == 1 && n <= 64)) return 1;
        for (int k = 0; k < n; k++) if (!rd(f, &x[k])) return 1;
        for (int k = 0; k < n; k++) if (!rd(f, &v[k])) return 1;
        int rc = dfl_engine_envelope(n, x, 1, alpha, DFL_ENGINE_AUTO, got, err, sizeof err);
        for (int k = 0; k < n; k++) { snprintf(what, sizeof what, "alpha = %g c = %g: %.15g vs %.15g %s", alpha, x[k], got[k], v[k], rc ? err : ""); check(tag, rc ? INFINITY : fabs(got[k] - v[k]), what); }
        return 0;
    }
    int a1, st, n; double g0, c, par;
    if (fscanf(f, "%d %d %lf %lf %lf %d", &a1, &st, &g0, &c, &par, &n) != 6 || n > 64) return 1;
    int mom = !strcmp(tag, "rcoefmom") || !strcmp(tag, "rmom"), nv = n;
    for (int k = 0; k < n && !mom; k++) if (!rd(f, &x[k])) return 1;
    for (int k = 0; k < nv; k++) if (!rd(f, &v[k])) return 1;
    if (mom) for (int k = 0; k < n; k++) x[k] = k + 1;
    int rc;
    snprintf(what, sizeof what, "%d start %d g0 = %g c = %g (%g)", a1, st, g0, c, par);
    size_t wl = strlen(what);
    dfl_spec S = {NULL, a1, st, g0, NAN};
    if (!strcmp(tag, "rcoef")) rc = dfl_engine_coef(n, x, &S, c, par, DFL_ENGINE_AUTO, got, err, sizeof err);
    else if (!strcmp(tag, "rcoefpdf")) rc = dfl_engine_coef_pdf(n, x, &S, c, par, DFL_ENGINE_AUTO, NULL, got, est, err, sizeof err);
    else if (!strcmp(tag, "rcoefq")) rc = dfl_engine_coef_icdf(n, x, &S, c, par, DFL_ENGINE_AUTO, 1e-9, got, err, sizeof err);
    else if (!strcmp(tag, "rcoefmom")) rc = dfl_engine_coef_mom(n, &S, c, par, DFL_ENGINE_AUTO, NULL, got, est, err, sizeof err);
    else if (!strcmp(tag, "rmom")) rc = dfl_engine_mom(n, &S, c, par, DFL_ENGINE_AUTO, got, err, sizeof err);
    else if (!strcmp(tag, "rlin")) {
        dfl_design D; dfl_design_trend(a1, &D, err, sizeof err);
        dfl_spec SD = {&D, 0, st, g0, NAN};
        rc = dfl_engine_lincdf_x(n, x, par, &SD, c, NAN, DFL_ENGINE_AUTO, NULL, got, est, err, sizeof err);
        if (!rc && a1 == -1 && st == 0) {
            double g2[64]; int r2 = dfl_engine_lincdf(n, x, par, 1, c, NAN, DFL_ENGINE_AUTO, g2, err, sizeof err);
            for (int k = 0; k < n; k++) check("self", r2 || g2[k] != got[k] ? 1 : 0, "the linear form of case DK by the case code and by trend(-1)");
        }
        dfl_design_free(&D);
    } else { printf("unknown tag %s\n", tag); return 1; }
    for (int k = 0; k < nv; k++) {
        double d = fabs(got[k] - v[k]);
        if (!strcmp(tag, "rcoefpdf")) d /= fmax(1.0, fabs(v[k]));
        if (!strcmp(tag, "rcoefq")) d /= fmax(1.0, fabs(v[k]));
        if (!strcmp(tag, "rcoefmom") || !strcmp(tag, "rmom")) d /= fmax(fabs(v[k]), nv > 1 ? pow(fabs(v[1]), (k + 1) / 2.0) : 1e-300);
        if (isnan(v[k]) && isnan(got[k])) d = 0;
        snprintf(what + wl, sizeof what - wl, " %s #%d (%g): %.15g vs %.15g %s", tag, k + 1, x[k], got[k], v[k], rc ? err : "");
        check(tag, rc ? INFINITY : d, what);
    }
    return 0;
}

/* a design of the lines of the general engine: D kind np par.. (kind 0 trend(par), 1, 2, 3 Perron's models A, B, C) */
static int read_design(FILE *f, dfl_design *D, char *desc, size_t dl)
{
    char t[8]; int kind, np; double par[8]; char err[256];
    if (fscanf(f, "%7s %d %d", t, &kind, &np) != 3 || strcmp(t, "D") || np < 1 || np > 8) return 1;
    for (int k = 0; k < np; k++) if (!rd(f, &par[k])) return 1;
    int rc = kind == 0 ? dfl_design_trend((int) par[0], D, err, sizeof err) : dfl_design_perron("ABC"[kind - 1], np, par, D, err, sizeof err);
    if (rc) { printf("design: %s\n", err); return 1; }
    if (kind == 0) snprintf(desc, dl, "trend(%g)", par[0]); else snprintf(desc, dl, "perron(%c, %g%s)", "ABC"[kind - 1], par[0], np > 1 ? ", ..." : "");
    return 0;
}
/* the lines of the general engine: 0, or 1 at the end of the file */
static int general_line(FILE *f, const char *tag)
{
    char err[512], what[300], ds[64]; dfl_design D; double x[64], v[64], got[64], est[64];
    if (read_design(f, &D, ds, sizeof ds)) return 1;
    int rc = 0, n = 0;
    if (!strcmp(tag, "dcf")) {
        int st; double gls, c, t1, t2; cdouble L;
        if (!(fscanf(f, "%d", &st) == 1 && rd(f, &gls) && rd(f, &c) && rd(f, &t1) && rd(f, &t2) && rdc(f, &L))) { dfl_design_free(&D); return 1; }
        det_prep *P = det_prep_cached(&D, gls, err, sizeof err); det_model Mo;
        if (!P || det_model_make(P, st, 0.0, &Mo, err, sizeof err)) { printf("FAIL dcf: %s\n", err); nfail++; dfl_design_free(&D); return 0; }
        det_cf *F = det_cf_new(&Mo, c, 0);
        cdouble Lc = det_cf_logcf(F, t1, t2);
        snprintf(what, sizeof what, "%s start %d gls %g c = %g at (%g, %g): %.15g%+.15gi vs %.15g%+.15gi", ds, st, gls, c, t1, t2, creal(Lc), cimag(Lc), creal(L), cimag(L));
        check(tag, cdist(Lc, L) / fmax(1.0, cabs(L)), what);
        det_cf_free(F); det_model_free(&Mo); det_prep_release(P); dfl_design_free(&D);
        return 0;
    }
    if (!strcmp(tag, "ptcdf") || !strcmp(tag, "ptpdf")) {
        double cbar, c; if (!(rd(f, &cbar) && rd(f, &c) && fscanf(f, "%d", &n) == 1 && n <= 64)) { dfl_design_free(&D); return 1; }
        for (int k = 0; k < n; k++) rd(f, &x[k]);
        for (int k = 0; k < n; k++) rd(f, &v[k]);
        rc = dfl_det_pt_x(tag[2] == 'c' ? 2 : 1, &D, n, x, cbar, c, 0, 0, 0, NAN, got, est, err, sizeof err);
        for (int k = 0; k < n; k++) {
            snprintf(what, sizeof what, "%s cbar %g c = %g q = %g: %.15g vs %.15g (estimate %.2g) %s", ds, cbar, c, x[k], got[k], v[k], est[k], rc ? err : "");
            check(tag, rc ? INFINITY : fabs(got[k] - v[k]) / (tag[2] == 'p' ? fmax(1.0, fabs(v[k])) : 1.0), what);
        }
    } else if (!strcmp(tag, "ptq")) {
        double cbar, p, q, qc, e; if (!(rd(f, &cbar) && rd(f, &p) && rd(f, &q))) { dfl_design_free(&D); return 1; }
        rc = dfl_det_pt_quantile_x(&D, p, cbar, 0, 0, 0, NAN, &qc, &e, err, sizeof err);
        snprintf(what, sizeof what, "%s cbar %g p = %g: %.15g vs %.15g %s", ds, cbar, p, qc, q, rc ? err : "");
        check(tag, rc ? INFINITY : fabs(qc - q) / fmax(1.0, fabs(q)), what);
    } else if (!strcmp(tag, "ptpow")) {
        double cbar, c, lev, b, P, bc, Pc; if (!(rd(f, &cbar) && rd(f, &c) && rd(f, &lev) && rd(f, &b) && rd(f, &P))) { dfl_design_free(&D); return 1; }
        rc = dfl_det_pt_power(&D, cbar, c, lev, 0, 0, &bc, &Pc, err, sizeof err);
        snprintf(what, sizeof what, "%s cbar %g c = %g: b %.15g vs %.15g, power %.15g vs %.15g %s", ds, cbar, c, bc, b, Pc, P, rc ? err : "");
        check(tag, rc ? INFINITY : fmax(fabs(bc - b) / fmax(1.0, fabs(b)), fabs(Pc - P)), what);
    } else if (!strcmp(tag, "ptenv")) {
        double lev; if (!(rd(f, &lev) && fscanf(f, "%d", &n) == 1 && n <= 64)) { dfl_design_free(&D); return 1; }
        for (int k = 0; k < n; k++) rd(f, &x[k]);
        for (int k = 0; k < n; k++) rd(f, &v[k]);
        rc = dfl_det_pt_envelope(&D, n, x, lev, 0, 0, got, err, sizeof err);
        for (int k = 0; k < n; k++) {
            snprintf(what, sizeof what, "%s level %g c = %g: %.15g vs %.15g %s", ds, lev, x[k], got[k], v[k], rc ? err : "");
            check(tag, rc ? INFINITY : (isnan(v[k]) && isnan(got[k]) ? 0.0 : fabs(got[k] - v[k])), what);
        }
    } else if (!strcmp(tag, "gp")) {
        int st; double gls, c, c0; if (!(fscanf(f, "%d", &st) == 1 && rd(f, &gls) && rd(f, &c) && rd(f, &c0) && fscanf(f, "%d", &n) == 1 && n <= 64)) { dfl_design_free(&D); return 1; }
        for (int k = 0; k < n; k++) rd(f, &x[k]);
        for (int k = 0; k < n; k++) rd(f, &v[k]);
        dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0; o.start = st; o.gls = gls;
        rc = dfl_det_coef_cdf_gp_x(&D, n, x, c, &o, 0, 0, 0, NAN, got, est, err, sizeof err);
        for (int k = 0; k < n; k++) {
            snprintf(what, sizeof what, "%s start %d gls %g c = %g c0 = %g y = %g: %.15g vs %.15g (estimate %.2g) %s", ds, st, gls, c, c0, x[k], got[k], v[k], est[k], rc ? err : "");
            check(tag, rc ? INFINITY : (isnan(v[k]) && isnan(got[k]) ? 0.0 : fabs(got[k] - v[k])), what);
        }
    }
    dfl_design_free(&D);
    return 0;
}
/* the router's rules without reference values: the regions of the defaults (the error estimates NaN where the engines
 * compute), method DFL_METHOD_ENGINE as the computation of before (bit for bit), the refusals */
#define SELF(cond, what) check("self", (cond) ? 0.0 : 1.0, what)
static void selftests(void)
{
    char err[512]; double v[8], e[8], w[8];
    dfl_spec DK = {NULL, 1, 0, 0.0, NAN}, DU = {NULL, 2, 0, 0.0, NAN}, SU = {NULL, 2, 2, 0.0, NAN}, DG = {NULL, 4, 0, 0.0, NAN};
    dfl_inv_opts io; dfl_inv_defaults(&io);
    double y[2] = {-9.0, -3.0}, x0 = 0;
    /* the coefficient statistic: the region c <= 20, the engines beyond and for case DG, method engine = the t-ratio at 0 */
    int rc = dfl_engine_coef_x(2, y, &DK, 20.0, 0.0, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isfinite(e[0]) && e[0] <= DFL_KAN_TOL, "coef_cdf at c = 20: the single inversion with its estimate");
    rc = dfl_engine_coef_x(1, y, &DK, 20.5, 0.0, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isnan(e[0]), "coef_cdf at c = 20.5: the engines (estimate NaN)");
    io.method = DFL_METHOD_ENGINE;
    rc = dfl_engine_coef_x(1, y, &DU, -5.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err) || dfl_engine_dist(2, 1, &x0, &DU, -5.0, -9.0, DFL_ENGINE_AUTO, 0, 0, w, err, sizeof err);
    SELF(!rc && v[0] == w[0] && isnan(e[0]), "coef_cdf with method engine is the t-ratio at 0 with the null value c0 + y");
    rc = dfl_engine_coef_x(1, y, &DU, -5.0, 0.0, DFL_ENGINE_CLASSICAL, NULL, w + 1, e, err, sizeof err);
    SELF(!rc && w[1] == w[0], "coef_cdf with engine classical: the engines");
    io.method = DFL_METHOD_INVERSION;
    rc = dfl_engine_coef_x(1, y, &DU, 25.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(rc || isfinite(e[0]), "coef_cdf with method inversion beyond c = 20");
    rc = dfl_engine_coef_x(1, y, &DU, -5.0, 0.0, DFL_ENGINE_GENERAL, &io, v, e, err, sizeof err);
    SELF(rc && strstr(err, "does not apply with engine general"), "coef_cdf: method inversion with engine general refused");
    rc = dfl_engine_coef_x(1, y, &DG, -5.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(rc && strstr(err, "applies to the cases DK, DU, DT, SK, SU, ST"), "coef_cdf: method inversion in case DG refused");
    io.method = 9;
    rc = dfl_engine_coef_x(1, y, &DU, -5.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(rc && strstr(err, "method must be"), "coef_cdf: an unknown method refused");
    /* verified keeps the values whose estimate meets the tolerance */
    dfl_inv_defaults(&io); io.verified = 1;
    rc = dfl_engine_coef_x(2, y, &SU, -10.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err) || dfl_engine_coef(2, y, &SU, -10.0, 0.0, DFL_ENGINE_AUTO, w, err, sizeof err);
    SELF(!rc && v[0] == w[0] && v[1] == w[1], "coef_cdf with verified: the values of the single inversion where their estimates meet the tolerance");
    /* the density with engine classical (which has none): the single inversion; the moments beyond c = 20 refused */
    dfl_inv_defaults(&io);
    rc = dfl_engine_coef_pdf(1, y, &DK, -5.0, 0.0, DFL_ENGINE_CLASSICAL, NULL, v, e, err, sizeof err) || dfl_engine_coef_pdf(1, y, &DK, -5.0, 0.0, DFL_ENGINE_AUTO, NULL, w, e + 1, err, sizeof err);
    SELF(!rc && v[0] == w[0] && isfinite(e[0]), "coef_pdf with engine classical: the single inversion");
    rc = dfl_engine_coef_mom(2, &DK, 21.0, 0.0, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(rc && strstr(err, "c must be at most 20"), "coef_mom at c = 21 refused with method 0");
    rc = dfl_engine_coef_mom(2, &DG, 0.0, 0.0, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(rc && strstr(err, "cases DK, DU, DT, SK, SU, ST"), "coef_mom in case DG refused");
    rc = dfl_engine_coef_mom(2, &DK, 0.0, 0.0, DFL_ENGINE_GENERAL, NULL, v, e, err, sizeof err);
    SELF(rc && strstr(err, "no moments of the coefficient statistic"), "coef_mom with engine general refused");
    /* the moments of the t-ratio: jets up to c = 20, the circles beyond and with method engine */
    rc = dfl_engine_mom_x(4, &DU, -5.0, 0.0, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isfinite(e[0]) && e[3] <= 1e-10 * fabs(v[3]), "ltratio_mom: the jets with their estimates");
    io.method = DFL_METHOD_ENGINE;
    rc = dfl_engine_mom_x(4, &DU, -5.0, 0.0, DFL_ENGINE_AUTO, &io, w, e, err, sizeof err);
    { double m2[4]; int r2 = dfl_engine_mom(4, &DU, -5.0, 0.0, DFL_ENGINE_GENERAL, m2, err, sizeof err);
      SELF(!rc && !r2 && isnan(e[0]) && w[0] == m2[0] && fabs(w[3] / v[3] - 1) < 1e-9, "ltratio_mom with method engine: the engines' circles (the rule takes the general engine at c0 = 0)"); }
    dfl_inv_defaults(&io);
    /* the linear form: case DK up to c = 4, the other cases through their designs up to c = 2, a = 0 and a = NaN by the engines */
    double b[2] = {-0.2639, 0.5};
    rc = dfl_engine_lincdf_x(1, b, -1.0, &DK, 4.0, NAN, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isfinite(e[0]), "ltratio_lincdf case DK at c = 4: the single inversion");
    rc = dfl_engine_lincdf_x(1, b, -1.0, &DK, 4.5, NAN, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isnan(e[0]), "ltratio_lincdf case DK at c = 4.5: the engines");
    rc = dfl_engine_lincdf_x(1, b, 0.0, &DK, -5.0, NAN, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err) || dfl_lincdf(1, b, 0.0, 1, -5.0, DFL_CBAR_GLS, NULL, w, err, sizeof err);
    SELF(!rc && isnan(e[0]) && v[0] == w[0], "ltratio_lincdf at a = 0: the engines' closed form");
    { dfl_design T0; dfl_design_trend(0, &T0, err, sizeof err); dfl_spec ST0 = {&T0, 0, 0, 0.0, NAN};
      rc = dfl_engine_lincdf_x(1, b, -3.5, &ST0, 2.0, NAN, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
      SELF(!rc && isfinite(e[0]), "ltratio_lincdf trend(0) at c = 2: the single inversion");
      rc = dfl_engine_lincdf_x(1, b, -3.5, &ST0, 2.5, NAN, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
      dfl_det_opts o; dfl_det_defaults(&o); int r2 = dfl_det_lin(2, &T0, 1, b, -3.5, 2.5, &o, w, err, sizeof err);
      SELF(!rc && !r2 && isnan(e[0]) && v[0] == w[0], "ltratio_lincdf trend(0) at c = 2.5: the general engine");
      dfl_design_free(&T0); }
    io.method = DFL_METHOD_INVERSION;
    rc = dfl_engine_lincdf_x(1, b, NAN, &DK, -5.0, NAN, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(rc && strstr(err, "needs a finite a"), "ltratio_lincdf: method inversion with a = NaN refused");
    rc = dfl_engine_lincdf_x(1, b, -3.5, &DG, -5.0, NAN, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(rc && strstr(err, "not to case DG"), "ltratio_lincdf: method inversion in case DG refused");
    dfl_inv_defaults(&io);
    /* the envelope: c <= -1e-4 by the single inversion, -1e-4 < c < 0 by the engines */
    double cs[3] = {-1e-4, -5e-5, -7.0};
    rc = dfl_engine_envelope_x(3, cs, 1, 0.05, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isfinite(e[0]) && isnan(e[1]) && isfinite(e[2]), "power_envelope: the single inversion for c <= -1e-4, the engines at c = -5e-5");
    rc = dfl_engine_envelope_x(1, cs + 2, 4, 0.05, DFL_ENGINE_AUTO, NULL, v, e, err, sizeof err);
    SELF(!rc && isnan(e[0]), "power_envelope case DG: the engines");
    /* the warning of the single inversion where an estimate exceeds the tolerance */
    io.tol = 1e-20;
    rc = dfl_engine_coef_x(1, y, &DK, -5.0, 0.0, DFL_ENGINE_AUTO, &io, v, e, err, sizeof err);
    SELF(!rc && dfl_det_warning() && strstr(dfl_det_warning(), "exceeds the tolerance"), "coef_cdf: the warning of an estimate above the tolerance");
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "test/kan_ref.txt", "r");
    if (!f) { printf("no reference file\n"); return 1; }
    char tag[32], err[512], what[256];
    double t0 = dfl_wtime();
    while (fscanf(f, "%31s", tag) == 1) {
        if (!strcmp(tag, "env")) {
            double alpha, c, P, e, Pc, ec;
            if (!(rd(f, &alpha) && rd(f, &c) && rd(f, &P) && rd(f, &e))) break;
            int rc = dfl_kan_envelope1(c, alpha, DFL_KAN_TOL, &Pc, &ec, err, sizeof err);
            snprintf(what, sizeof what, "alpha = %g, c = %g: %.15g vs %.15g", alpha, c, Pc, P);
            check(tag, rc ? INFINITY : fabs(Pc - P), what);
            continue;
        }
        if (tag[0] == 'r' && strcmp(tag, "ray")) { if (router_line(f, tag)) break; continue; }
        if (!strcmp(tag, "dcf") || !strncmp(tag, "pt", 2) || !strcmp(tag, "gp")) { if (general_line(f, tag)) break; continue; }
        int m, st; double c, g0;
        if (fscanf(f, "%d %d %lf %lf", &m, &st, &c, &g0) != 4) break;
        dfl_limcf F;
        if (dfl_limcf_init(&F, m, st, c, g0, err, sizeof err)) { printf("FAIL %s: %s\n", tag, err); nfail++; break; }
        snprintf(what, sizeof what, "model %d stat %d c = %g g0 = %g", m, st, c, g0);
        size_t wl = strlen(what);
        if (!strcmp(tag, "lcf")) {
            cdouble t1, t2, L;
            if (!(rdc(f, &t1) && rdc(f, &t2) && rdc(f, &L))) break;
            cdouble Lc = dfl_kan_logcf(&F, t1, t2);
            snprintf(what + wl, sizeof what - wl, " at (%g%+gi, %g%+gi): %.15g%+.15gi vs %.15g%+.15gi", creal(t1), cimag(t1), creal(t2), cimag(t2), creal(Lc), cimag(Lc), creal(L), cimag(L));
            check(tag, cdist(Lc, L) / fmax(1.0, cabs(L)), what);
        } else if (!strcmp(tag, "jet")) {
            int n; cdouble a[DFL_KAN_JMAX], b[DFL_KAN_JMAX], L[DFL_KAN_JMAX], Lc[DFL_KAN_JMAX]; int ok = fscanf(f, "%d", &n) == 1 && n >= 1 && n <= DFL_KAN_JMAX;
            for (int k = 0; ok && k < n; k++) ok = rdc(f, &a[k]);
            for (int k = 0; ok && k < n; k++) ok = rdc(f, &b[k]);
            for (int k = 0; ok && k < n; k++) ok = rdc(f, &L[k]);
            if (!ok) break;
            dfl_kan_logcf_jet(&F, n, a, b, Lc);
            double d = 0; for (int k = 0; k < n; k++) d = fmax(d, cdist(Lc[k], L[k]) / fmax(1.0, cabs(L[k])));
            snprintf(what + wl, sizeof what - wl, " jet of order %d at (%g, %g)", n, creal(a[0]), creal(b[0]));
            check(tag, d, what);
        } else if (!strcmp(tag, "scale")) {
            double s[4], sc[4]; if (!(rd(f, &s[0]) && rd(f, &s[1]) && rd(f, &s[2]) && rd(f, &s[3]))) break;
            dfl_kan_scale(&F, sc);
            double d = 0; for (int k = 0; k < 4; k++) d = fmax(d, fabs(sc[k] - s[k]) / fmax(1.0, fabs(s[k])));
            check(tag, d, what);
        } else if (!strcmp(tag, "coef") || !strcmp(tag, "coefpdf")) {
            double c0, y, v, e, vc, ec, ft; if (!(rd(f, &c0) && rd(f, &y) && rd(f, &v) && rd(f, &e))) break;
            int rc = tag[4] ? dfl_kan_coef_pdf1(&F, y, c0, DFL_KAN_TOL, &vc, &ec, &ft) : dfl_kan_coef_cdf1(&F, y, c0, DFL_KAN_TOL, &vc, &ec);
            snprintf(what + wl, sizeof what - wl, " c0 = %g y = %.10g: %.15g vs %.15g (estimates %.2g, %.2g)", c0, y, vc, v, ec, e);
            check(tag, rc ? INFINITY : fabs(vc - v) / (tag[4] ? fmax(1.0, fabs(v)) : 1.0), what);
        } else if (!strcmp(tag, "coefq")) {
            double c0, p, y, yc, ec; if (!(rd(f, &c0) && rd(f, &p) && rd(f, &y))) break;
            int rc = dfl_kan_coef_icdf1(&F, p, c0, 1e-9, DFL_KAN_TOL, &yc, &ec, err, sizeof err);
            snprintf(what + wl, sizeof what - wl, " p = %g: %.15g vs %.15g %s", p, yc, y, rc ? err : "");
            check(tag, rc ? INFINITY : fabs(yc - y) / fmax(1.0, fabs(y)), what);
        } else if (!strcmp(tag, "lin")) {
            double a, b, P, e, Pc, ec; if (!(rd(f, &a) && rd(f, &b) && rd(f, &P) && rd(f, &e))) break;
            int rc = dfl_kan_lin_cdf1(&F, b, a, DFL_KAN_TOL, &Pc, &ec);
            snprintf(what + wl, sizeof what - wl, " a = %g b = %.10g: %.15g vs %.15g", a, b, Pc, P);
            check(tag, rc ? INFINITY : fabs(Pc - P), what);
        } else if (!strcmp(tag, "tmom") || !strcmp(tag, "cmom")) {
            double c0, mr[DFL_KAN_JMAX], mc[DFL_KAN_JMAX], ec[DFL_KAN_JMAX]; int k;
            if (!(rd(f, &c0) && fscanf(f, "%d", &k) == 1 && k >= 1 && k < DFL_KAN_JMAX)) break;
            int ok = 1; for (int j = 0; ok && j < k; j++) ok = rd(f, &mr[j]);
            if (!ok) break;
            int rc = dfl_kan_ratio_mom(&F, k, c0 - F.c, tag[0] == 't' ? 0.5 : 1.0, 1e-10, mc, ec, err, sizeof err);
            double d = 0;
            for (int j = 0; j < k; j++) { double sc = fmax(fabs(mr[j]), k > 1 ? pow(fabs(mr[1]), (j + 1) / 2.0) : 0.0); d = fmax(d, fabs(mc[j] - mr[j]) / fmax(sc, 1e-300)); }
            snprintf(what + wl, sizeof what - wl, " c0 = %g: m_1 %.15g vs %.15g %s", c0, mc[0], mr[0], rc ? err : "");
            check(tag, rc ? INFINITY : d, what);
        } else if (!strcmp(tag, "ray")) {
            double P, Pc; if (!rd(f, &P)) break;
            int rc = dfl_kan_ray_cdf(m, st, c, g0, &Pc, err, sizeof err);
            snprintf(what + wl, sizeof what - wl, ": %.15g vs %.15g", Pc, P);
            check(tag, rc ? INFINITY : fabs(Pc - P), what);
        } else { printf("unknown tag %s\n", tag); return 1; }
    }
    fclose(f);
    selftests();
    for (size_t i = 0; i < sizeof KINDS / sizeof KINDS[0]; i++)
        if (KINDS[i].n) printf("  %-10s %5d checks, %d failed, worst difference %.2e (tolerance %.0e)\n", KINDS[i].tag, KINDS[i].n, KINDS[i].fail, KINDS[i].worst, KINDS[i].tol);
    printf("%d checks, %d failed (%.1f s, %d threads)\n", ntot, nfail, dfl_wtime() - t0, dfl_nworkers());
    return nfail ? 1 : 0;
}
