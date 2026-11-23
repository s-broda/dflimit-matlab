/* dfl_test_kan.c: the devices of Bao and Kan (dfl_kan.c, dfl_det_kan.c, the router's defaults in dfl_det_route.c) against
 * the values of the Julia package DFLimit.jl (test/kan_ref.txt, written by test/kan_ref.jl): the closed-form transform of the
 * six classical cases at points and in Taylor-mode arithmetic, the scale, the single inversion of the coefficient statistic
 * (distribution function, density, quantiles), of the linear form U - aV and of the power envelope, the moments of the
 * t-ratio and of the coefficient statistic from jets, the closed form at zero; each kind of line with its own tolerance
 * (at least ten times the largest difference measured on macOS, November 2026), and the largest difference per kind is
 * printed.  The program ends with a line that counts the checks and the failures. */
#include "dfl_kan_int.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { const char *tag; double tol, worst; int n, fail; } kind_t;
static kind_t KINDS[] = {
    {"lcf", 1e-12, 0, 0, 0}, {"jet", 1e-11, 0, 0, 0}, {"scale", 1e-12, 0, 0, 0}, {"coef", 1e-11, 0, 0, 0}, {"coefpdf", 1e-11, 0, 0, 0},
    {"coefq", 1e-8, 0, 0, 0}, {"lin", 1e-11, 0, 0, 0}, {"env", 1e-11, 0, 0, 0}, {"tmom", 1e-10, 0, 0, 0}, {"cmom", 1e-10, 0, 0, 0},
    {"ray", 1e-13, 0, 0, 0},
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
    for (size_t i = 0; i < sizeof KINDS / sizeof KINDS[0]; i++)
        if (KINDS[i].n) printf("  %-10s %5d checks, %d failed, worst difference %.2e (tolerance %.0e)\n", KINDS[i].tag, KINDS[i].n, KINDS[i].fail, KINDS[i].worst, KINDS[i].tol);
    printf("%d checks, %d failed (%.1f s, %d threads)\n", ntot, nfail, dfl_wtime() - t0, dfl_nworkers());
    return nfail ? 1 : 0;
}
