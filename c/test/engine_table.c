/* engine_table.c: the measurements behind the choice of engine of the C router (dfl_autoengine): the cells of
 * dev/engine_table.jl of DFLimit.jl, each computed by the classical C engine and by the general C engine (dfl_det_*), the
 * warm time of one call (the minimum of 7 timed calls after a warm-up, of 3 where a call takes 1 s or more, the first call
 * alone, marked *, where it takes more than 60 s; the general engine's outer path cached), the largest difference of
 * the two, and the engine of dfl_autoengine in the cell.  Writes the table to stdout (c/test/engine_table_c.txt).
 *   cc -std=c11 -O2 -I c -c c/test/engine_table.c ... (link with the objects of c/build.sh); DFL_NTHREADS=4 ./engine_table */
#include "dflimit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

enum { T_DIST, T_MOM, T_ICDF, T_COEF, T_COEFQ, T_LIN, T_ENV, T_ADF, T_REGION };
typedef struct { int type, model, stat; double c, c0, g0, a, cbar; int n; double x[8]; int rtype; double par[2]; } cell_t;
static char ERR[512];
static int run(const cell_t *e, int engine, double *out, int *nout)
{
    dfl_spec S = {NULL, e->model, e->stat ? DFL_START_STAT : (e->g0 != 0 ? DFL_START_FIXED : DFL_START_ZERO), e->g0, NAN};
    int rc = 0, n = e->n;
    switch (e->type) {
    case T_DIST:
        rc = dfl_engine_dist(2, n, e->x, &S, e->c, e->c0, engine, 0, 0, out, ERR, sizeof ERR) || dfl_engine_dist(1, n, e->x, &S, e->c, e->c0, engine, 0, 0, out + n, ERR, sizeof ERR);
        *nout = 2 * n; break;
    case T_MOM: rc = dfl_engine_mom(4, &S, e->c, e->c0, engine, out, ERR, sizeof ERR); *nout = 4; break;
    case T_ICDF: rc = dfl_engine_icdf(n, e->x, &S, e->c, e->c0, engine, 0, 1e-10, out, ERR, sizeof ERR); *nout = n; break;
    case T_COEF: rc = dfl_engine_coef(n, e->x, &S, e->c, e->c0, engine, out, ERR, sizeof ERR); *nout = n; break;
    case T_COEFQ: rc = dfl_engine_coef_icdf(n, e->x, &S, e->c, e->c0, engine, 1e-9, out, ERR, sizeof ERR); *nout = n; break;
    case T_LIN: rc = dfl_engine_lincdf(n, e->x, e->a, e->model, e->c, e->cbar, engine, out, ERR, sizeof ERR); *nout = n; break;
    case T_ENV: rc = dfl_engine_envelope(n, e->x, e->model, 0.05, engine, out, ERR, sizeof ERR); *nout = n; break;
    case T_ADF:
        rc = dfl_engine_adfmax(2, n, e->x, &S, e->c, engine, out, ERR, sizeof ERR) || dfl_engine_adfmax(1, n, e->x, &S, e->c, engine, out + n, ERR, sizeof ERR);
        *nout = 2 * n; break;
    case T_REGION: rc = dfl_engine_region(e->rtype, e->par, &S, e->c, engine, out, ERR, sizeof ERR); *nout = 1; break;
    }
    return rc;
}
static const char *casename(int model, int stat) { return model == 4 ? "DG" : stat ? (model == 1 ? "SK" : model == 2 ? "SU" : "ST") : (model == 1 ? "DK" : model == 2 ? "DU" : "DT"); }
static void row(const char *label, const cell_t *e)
{
    double vc[16], vg[16], t; int nc, ng;
    char lab[200]; int k = snprintf(lab, sizeof lab, "%s (%s, c = %g", label, casename(e->model, e->stat), e->c);
    if (e->c0 != 0) k += snprintf(lab + k, sizeof lab - k, ", c0 = %g", e->c0);
    if (e->g0 != 0) k += snprintf(lab + k, sizeof lab - k, ", g0 = %g", e->g0);
    if (e->type == T_LIN) k += snprintf(lab + k, sizeof lab - k, ", a = %g%s", e->a, e->model == 4 ? "" : "");
    if (e->type == T_LIN && e->model == 4) k += snprintf(lab + k, sizeof lab - k, ", cbar = %g", e->cbar);
    snprintf(lab + k, sizeof lab - k, ")");
    double t0c = dfl_wtime(); int rcc = run(e, DFL_ENGINE_CLASSICAL, vc, &nc); t0c = dfl_wtime() - t0c;
    /* (DFL_TABLE_SKIP_GENERAL: the general engine is not run in the cells whose label contains this text; the clean pass of
     * November 26, 2026 skipped the bulk of case DU at c = 10, where one call of the general engine had not finished after
     * 16 minutes, single-threaded in its line integral's provider, against 23 s in the table of November 22) */
    const char *skg = getenv("DFL_TABLE_SKIP_GENERAL");
    if (skg && *skg && strstr(lab, skg)) {
        double tc = INFINITY;
        for (int r = 0, nr = 7; r < nr; r++) { t = dfl_wtime(); run(e, DFL_ENGINE_CLASSICAL, vc, &nc); t = dfl_wtime() - t; if (t < tc) tc = t; if (t >= 1) nr = 3; }
        printf("  %-66s classical %9.2f  general   skipped (DFL_TABLE_SKIP_GENERAL)\n", lab, 1e3 * tc); fflush(stdout); (void) rcc; return;
    }
    double t0g = dfl_wtime(); int rcg = run(e, DFL_ENGINE_GENERAL, vg, &ng); t0g = dfl_wtime() - t0g;
    if (rcc || rcg) { printf("  %-66s error: %s\n", lab, ERR); fflush(stdout); return; }
    double d = 0;
    for (int i = 0; i < nc; i++) {
        double di = fabs(vc[i] - vg[i]);
        if (e->type == T_MOM) di = fabs(vg[i] / vc[i] - 1);
        if (!(di <= d)) d = di;
    }
    /* (a call of more than 60 s is timed once, by the first call, and marked *, as in engine_table_matlab.m) */
    double tc = INFINITY, tg = INFINITY;
    if (t0c > 60) tc = t0c; else for (int r = 0, nr = 7; r < nr; r++) { t = dfl_wtime(); run(e, DFL_ENGINE_CLASSICAL, vc, &nc); t = dfl_wtime() - t; if (t < tc) tc = t; if (t >= 1) nr = 3; }
    if (t0g > 60) tg = t0g; else for (int r = 0, nr = 7; r < nr; r++) { t = dfl_wtime(); run(e, DFL_ENGINE_GENERAL, vg, &ng); t = dfl_wtime() - t; if (t < tg) tg = t; if (t >= 1) nr = 3; }
    int f = e->type == T_DIST || e->type == T_COEF || e->type == T_COEFQ ? DFL_F_DIST : e->type == T_MOM ? DFL_F_MOM : e->type == T_ICDF ? DFL_F_ICDF :
            e->type == T_LIN ? DFL_F_LINCDF : e->type == T_ENV ? DFL_F_ENVELOPE : e->type == T_ADF ? DFL_F_ADFMAX : DFL_F_REGION;
    double xa = e->type == T_DIST ? e->x[e->n / 2] : e->type == T_LIN ? e->a : 0.0;
    double c0a = e->type == T_COEF ? e->c0 - 10.0 : e->c0;
    /* (the regions as dfl_engine_region asks the rule: the parabola of case DK with x = 1 and its c0 = par[1], the line and the
     * other boundaries by the classical code; until November 27, 2026 the table asked with x = 0 and labelled the parabola
     * classical, which the router takes by the general engine at c0 = 0 for -20 <= c <= 0) */
    int au = e->type == T_REGION
        ? ((e->model == 1 && !e->stat && e->rtype == 2 && dfl_autoengine(DFL_F_REGION, 1, 0, 0.0, NAN, e->c, e->par[1], 1.0) == DFL_ENGINE_GENERAL) ? DFL_ENGINE_GENERAL : DFL_ENGINE_CLASSICAL)
        : dfl_autoengine(f, e->model, e->stat, e->g0, e->type == T_LIN ? e->cbar : (e->model == 4 ? DFL_CBAR_GLS : NAN), e->type == T_ENV ? e->x[0] : e->c, c0a, xa);
    printf("  %-66s classical %9.2f%s general %9.2f%s ratio g/c %7.2f  diff %8.1e  auto %s\n", lab, 1e3 * tc, t0c > 60 ? "*" : " ", 1e3 * tg, t0g > 60 ? "*" : " ", tg / tc, d, au == DFL_ENGINE_GENERAL ? "general" : "classical");
    fflush(stdout);
}
/* F and f at mean + (-2, -1, 0, 1, 2) sd (the classical engine's moments) */
static void tcell(const char *label, int model, int stat, double c, double c0, double g0, int nks, const double *ks, const double *xs)
{
    cell_t e; memset(&e, 0, sizeof e); e.type = T_DIST; e.model = model; e.stat = stat; e.c = c; e.c0 = c0; e.g0 = g0;
    if (xs) { e.n = 5; memcpy(e.x, xs, 5 * sizeof(double)); }
    else {
        double m[2], me[2];
        if (dfl_moments(model, c, c0, g0, stat, 2, m, me, ERR, sizeof ERR)) { printf("  %s: moments: %s\n", label, ERR); return; }
        double sd = sqrt(fmax(m[1] - m[0] * m[0], 1e-12));
        e.n = nks; for (int i = 0; i < nks; i++) e.x[i] = m[0] + ks[i] * sd;
    }
    row(label, &e);
}
static void borders(void);
int main(int argc, char **argv)
{
    double k5[5] = {-2, -1, 0, 1, 2}, k3[3] = {-1, 0, 1};
    if (argc > 1 && !strcmp(argv[1], "borders")) { borders(); return 0; }
    printf("The choice of engine of the C router: warm times in ms (the minimum of 7 timed calls, 3 from 1 s on; * one call of more than 60 s), %d threads; diff the largest difference of the two\n", dfl_nworkers());
    printf("C engines (relative for the moments); auto: the engine of dfl_autoengine in the cell (the cells of dev/engine_table.txt of DFLimit.jl).\n\n");
    printf("(1) the t-ratio, F and f at x = mean + (-2, -1, 0, 1, 2) sd\n");
    const char *L1 = "zero start, c0 = 0, -100 <= c <= 5";
    for (int m = 1; m <= 3; m++) { double cs[3] = {0, -10, -50}; for (int i = 0; i < 3; i++) tcell(L1, m, 0, cs[i], 0, 0, 5, k5, NULL); }
    { double x1[5] = {-3.5, -2.8, -2.0, -1.2, -0.5}, x2[5] = {-4.5, -3.5, -2.5, -1.5, -0.5}; tcell(L1, 4, 0, 0, 0, 0, 5, k5, x1); tcell(L1, 4, 0, -10, 0, 0, 5, k5, x2); }
    tcell("zero start, c0 = 0, 0 < c <= 5", 3, 0, 3, 0, 0, 5, k5, NULL); tcell("zero start, c0 = 0, 0 < c <= 5", 3, 0, 5, 0, 0, 5, k5, NULL);
    for (int m = 1; m <= 3; m++) tcell("zero start, c0 = 0, c < -100", m, 0, -300, 0, 0, 5, k5, NULL);
    { int ms[6] = {2, 1, 3, 2, 1, 3}; double cs[6] = {-10, -5, -5, -50, -20, -100}; for (int i = 0; i < 6; i++) tcell("zero start, c0 < 0, -100 <= c <= 0", ms[i], 0, cs[i], cs[i], 0, 5, k5, NULL); }
    { int ms[8] = {3, 2, 3, 1, 1, 1, 2, 3}; double cs[8] = {-20, -5, -50, -10, -2, 0, 0, 1}, c0s[8] = {3, 2, 5, 2, 1, 2, 2, 3}; for (int i = 0; i < 8; i++) tcell("zero start, c0 > 0, -100 < c <= 1", ms[i], 0, cs[i], c0s[i], 0, 5, k5, NULL); }
    { int ms[3] = {2, 2, 1}; double cs[3] = {-100, 2, 2}, c0s[3] = {3, 1, 2}; for (int i = 0; i < 3; i++) tcell("zero start, c0 > 0, otherwise", ms[i], 0, cs[i], c0s[i], 0, 5, k5, NULL); }
    for (int m = 1; m <= 3; m++) tcell("zero start, c0 != 0, c < -100", m, 0, -300, -300, 0, 5, k5, NULL);
    tcell("zero start, c0 != 0, c < -100", 2, 0, -480, -480, 0, 5, k5, NULL);
    tcell("stationary start, c0 = 0", 2, 1, -10, 0, 0, 5, k5, NULL); tcell("stationary start, c0 = 0", 3, 1, -50, 0, 0, 5, k5, NULL);
    { int ms[5] = {1, 2, 3, 2, 1}; double cs[5] = {-5, -10, -50, -100, -300}; for (int i = 0; i < 5; i++) tcell("stationary start, c0 = c, c <= -5", ms[i], 1, cs[i], cs[i], 0, 5, k5, NULL); }
    { int ms[6] = {1, 2, 1, 2, 2, 3}; double cs[6] = {-1, -0.5, -0.1, -0.1, -3e-3, -3e-3}; for (int i = 0; i < 6; i++) tcell("stationary start, c0 = c, -5 < c <= -3e-3", ms[i], 1, cs[i], cs[i], 0, 5, k5, NULL); }
    { int ms[3] = {2, 3, 1}; double cs[3] = {-1e-3, -1e-3, -1e-4}; for (int i = 0; i < 3; i++) tcell("stationary start, c0 = c, c > -3e-3", ms[i], 1, cs[i], cs[i], 0, 5, k5, NULL); }
    tcell("stationary start, c0 > 0", 3, 1, -1, 2, 0, 5, k5, NULL);
    { int ms[4] = {1, 2, 3, 1}; double cs[4] = {5, 5, 5, 4}, c0s[4] = {-3, -3, -3, -0.5}; for (int i = 0; i < 4; i++) tcell("above the unit root, c0 < 0, c >= 4", ms[i], 0, cs[i], c0s[i], 0, 3, k3, NULL); }
    { int ms[3] = {3, 1, 2}; double cs[3] = {7, 7, 10}; for (int i = 0; i < 3; i++) tcell("above the unit root, c >= 7, the bulk", ms[i], 0, cs[i], 0, 0, 3, k3, NULL); }
    { double cs[5] = {0, -0.5, -5, -50, 1}, gs[5] = {1, 0.5, 2, 0.5, 1}; for (int i = 0; i < 5; i++) tcell("fixed start, case DK, c <= 1, (c g0)^2/2 <= 400", 1, 0, cs[i], 0, gs[i], 5, k5, NULL); }
    tcell("fixed start, case DK, c <= 1, (c g0)^2/2 <= 400", 1, 0, -10, -2, 1, 5, k5, NULL);
    tcell("fixed start, case DK, c >= -100, g0 <= 0.5", 1, 0, -100, 0, 0.5, 5, k5, NULL); tcell("fixed start, case DK, c >= -100, g0 <= 0.5", 1, 0, -80, 0, 0.3, 5, k5, NULL);
    { double cs[3] = {2, -20, -100}, gs[3] = {1, 2, 1}; for (int i = 0; i < 3; i++) tcell("fixed start, case DK, otherwise", 1, 0, cs[i], 0, gs[i], 5, k5, NULL); }
    { int ms[3] = {2, 3, 2}; double cs[3] = {-10, 0, -50}, gs[3] = {1, 2, 1.5}; for (int i = 0; i < 3; i++) tcell("fixed start, cases DU, DT", ms[i], 0, cs[i], 0, gs[i], 5, k5, NULL); }

    printf("\n(2) the moments E t^k, k = 1..4 (diff: relative)\n");
    {
        int ms[18] = {1, 2, 3, 3, 1, 2, 3, 2, 2, 2, 1, 1, 3, 2, 1, 2, 3, 2}, st[18] = {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0};
        double cs[18] = {-5, -5, -5, 0, -5, -5, -5, -5, -50, -100, -480, -300, -200, -5, -100, -480, -300, -5}, gs[18] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double c0s[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -300, -200, -5, -100, -480, -300, -5};
        for (int i = 0; i < 18; i++) {
            cell_t e; memset(&e, 0, sizeof e); e.type = T_MOM; e.model = ms[i]; e.stat = st[i]; e.c = cs[i]; e.c0 = c0s[i]; e.g0 = gs[i];
            row(i == 7 ? "moments, fixed start" : i >= 11 && i <= 13 ? "moments, stationary, c0 != 0" : i >= 14 ? "moments, zero start, c0 != 0" : "moments", &e);
        }
    }
    printf("\n(3) quantiles and the coefficient statistic\n");
    {
        int ms[11] = {2, 2, 3, 3, 1, 2, 1, 2, 2, 3, 4}, st[11] = {0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0};
        double cs[11] = {0, -10, -5, -10, -5, -5, 0, -50, -5, 3, 0}, gs[11] = {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
        for (int i = 0; i < 11; i++) {
            cell_t e; memset(&e, 0, sizeof e); e.type = T_ICDF; e.model = ms[i]; e.stat = st[i]; e.c = cs[i]; e.g0 = gs[i]; e.n = 3; e.x[0] = 0.01; e.x[1] = 0.05; e.x[2] = 0.5;
            row("ltratio_icdf at p = 0.01, 0.05, 0.5", &e);
        }
        int cm[3] = {2, 2, 3}, cst[3] = {0, 1, 0}; double cc[3] = {0, -5, -10};
        for (int i = 0; i < 3; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_COEF; e.model = cm[i]; e.stat = cst[i]; e.c = cc[i]; e.n = 3; e.x[0] = -20; e.x[1] = -10; e.x[2] = -5; row("coef_cdf at y = -20, -10, -5", &e); }
        cell_t e; memset(&e, 0, sizeof e); e.type = T_COEFQ; e.model = 2; e.c = 0; e.n = 1; e.x[0] = 0.05; row("coef_icdf at p = 0.05", &e);
    }
    printf("\n(4) the linear form P(U - aV <= b) and the law of the denominator P(V <= b)\n");
    {
        struct { int m; double a, c, cbar; int n; double b[4]; } L[] = {
            {1, -3.5, 0, -13.5, 3, {-0.6, -0.2639, 0.3}}, {1, -3.5, -7, -13.5, 3, {-0.6, -0.2639, 0.3}}, {1, -1, 0, -13.5, 3, {-0.4, 0, 0.5}}, {1, 2, -5, -13.5, 3, {-2, 0, 0.5}},
            {1, 0, 0, -13.5, 4, {-0.499, -0.4, 0, 1}}, {1, -0.1, 0, -13.5, 4, {-0.49, -0.4, 0, 0.5}}, {1, 0, 3, -13.5, 3, {-0.4, 0.5, 3}}, {4, -7.3, 0, -13.5, 3, {-0.9, -0.5, 0}},
            {4, 0, 0, -13.5, 3, {-0.49, -0.3, 0.5}}, {1, NAN, -5, -13.5, 3, {0.05, 0.2, 0.6}}, {1, NAN, 2, -13.5, 3, {0.05, 0.2, 0.6}}, {1, NAN, 5, -13.5, 3, {0.02, 0.06, 0.1}},
            {1, NAN, 7, -13.5, 3, {0.02, 0.06, 0.1}}, {1, -1, 5, -13.5, 3, {-0.3, 0, 50}}, {1, 2.5, 5, -13.5, 3, {-160, -0.49, 0}},
            {4, -2, 0, -5, 4, {-0.9, -0.5, 0, 0.5}}, {4, -7.3, -5, -2, 4, {-0.9, -0.5, 0, 0.5}}, {4, -2, 0, -1, 4, {-0.9, -0.5, 0, 0.5}}, {4, -2, 2, -1, 4, {-0.9, -0.5, 0, 0.5}} };
        for (size_t i = 0; i < sizeof L / sizeof L[0]; i++) {
            cell_t e; memset(&e, 0, sizeof e); e.type = T_LIN; e.model = L[i].m; e.a = L[i].a; e.c = L[i].c; e.cbar = L[i].cbar; e.n = L[i].n; memcpy(e.x, L[i].b, sizeof L[i].b);
            row(isnan(L[i].a) ? "ltratio_lincdf, a = NaN (P(V <= b))" : "ltratio_lincdf", &e);
        }
    }
    printf("\n(5) the power envelope, ADFmax, regions\n");
    {
        struct { int m, n; double c[3]; } E[] = { {1, 3, {-5, -10, -20}}, {1, 2, {-0.01, -0.5}}, {4, 3, {-5, -10, -20}}, {4, 2, {-1, -2}} };
        for (int i = 0; i < 4; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_ENV; e.model = E[i].m; e.n = E[i].n; memcpy(e.x, E[i].c, sizeof E[i].c); row("power_envelope", &e); }
        int am[6] = {2, 3, 2, 3, 2, 3}, ast[6] = {0, 0, 1, 1, 1, 1}; double ac[6] = {0, -10, -10, -5, -0.1, -50};
        for (int i = 0; i < 6; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_ADF; e.model = am[i]; e.stat = ast[i]; e.c = ac[i]; e.n = 3; e.x[0] = -3; e.x[1] = -2.4; e.x[2] = -1.5; row("adfmax_cdf, adfmax_pdf at x = -3, -2.4, -1.5", &e); }
        cell_t e; memset(&e, 0, sizeof e); e.type = T_REGION; e.model = 1; e.c = -5; e.rtype = 1; e.par[0] = -3.5; e.par[1] = -0.2639; row("ltratio_region, line a = -3.5, b = -0.2639", &e);
        e.c = -7; e.rtype = 2; e.par[0] = -1.94087; e.par[1] = 0; row("ltratio_region, parabola x = -1.94087", &e);
    }
    borders();
    return 0;
}
/* (6) the cells that place the borders where the C timings move the rule */
static void borders(void)
{
    double k5[5] = {-2, -1, 0, 1, 2};
    printf("\n(6) the borders of the C rule\n");
    { double cs[4] = {-20, -30, -40, -100}; for (int m = 1; m <= 3; m++) for (int i = 0; i < 4; i++) tcell("zero start, c0 = 0, far from the unit root", m, 0, cs[i], 0, 0, 5, k5, NULL); }
    { double cs[4] = {-20, -30, -50, -100}; for (int m = 1; m <= 3; m++) for (int i = 0; i < 4; i++) tcell("stationary start, c0 = 0, far from the unit root", m, 1, cs[i], 0, 0, 5, k5, NULL); }
    { double cs[3] = {-2, -20, -30}; for (int m = 2; m <= 3; m++) for (int i = 0; i < 3; i++) tcell("zero start, c0 = c, cases DU, DT", m, 0, cs[i], cs[i], 0, 5, k5, NULL); }
    { int ms[4] = {1, 1, 2, 3}, st[4] = {1, 0, 1, 0}; double cs[4] = {-5, -10, -1, 2};
      for (int i = 0; i < 4; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_ICDF; e.model = ms[i]; e.stat = st[i]; e.c = cs[i]; e.n = 3; e.x[0] = 0.01; e.x[1] = 0.05; e.x[2] = 0.5; row("ltratio_icdf at p = 0.01, 0.05, 0.5", &e); } }
    { int ms[3] = {2, 1, 2}, st[3] = {1, 0, 0}; double cs[3] = {-50, -5, -50}, c0s[3] = {-50, -5, -50};
      for (int i = 0; i < 3; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_MOM; e.model = ms[i]; e.stat = st[i]; e.c = cs[i]; e.c0 = c0s[i]; row("moments, c0 != 0", &e); } }
    { double cs[5] = {1, 2, 3, 4, 2}, cb[5] = {-13.5, -13.5, -13.5, -13.5, -5};
      for (int i = 0; i < 5; i++) { cell_t e; memset(&e, 0, sizeof e); e.type = T_LIN; e.model = 4; e.a = -7.3; e.c = cs[i]; e.cbar = cb[i]; e.n = 3; e.x[0] = -0.9; e.x[1] = -0.5; e.x[2] = 0; row("ltratio_lincdf", &e); } }
}
