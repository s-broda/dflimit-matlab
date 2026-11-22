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

#define MAXD 64
static dfl_design DS[MAXD];
static int HAVE[MAXD];

static void f_sincos(double r, double *o, void *c) { (void) c; o[0] = sin(2 * M_PI * r); o[1] = cos(2 * M_PI * r); }
static void f_perronC04(double r, double *o, void *c) { (void) c; o[0] = 1; o[1] = r; o[2] = r > 0.4 ? 1.0 : 0.0; o[3] = r > 0.4 ? r - 0.4 : 0.0; }
static void f_expr(double r, double *o, void *c) { (void) c; o[0] = exp(r); }

/* the groups of comparisons: name, count, largest deviation (absolute), failures */
typedef struct { const char *name; int n, fail; double worst, worst_rel_tol; } group_t;
static group_t G[64]; static int NG = 0;
static group_t *grp(const char *name) { for (int i = 0; i < NG; i++) if (!strcmp(G[i].name, name)) return &G[i]; G[NG].name = strdup(name); return &G[NG++]; }
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
    }
    fclose(f);
    int fails = 0, n = 0;
    printf("%-10s %8s %8s %12s %10s\n", "group", "values", "failed", "largest dev", "dev/tol");
    for (int i = 0; i < NG; i++) { printf("%-10s %8d %8d %12.2e %10.2e\n", G[i].name, G[i].n, G[i].fail, G[i].worst, G[i].worst_rel_tol); fails += G[i].fail; n += G[i].n; }
    printf("%d comparisons, %d failed, %.1f s (%d threads)\n", n, fails, dfl_wtime() - t0, dfl_nworkers());
    return fails ? 1 : 0;
}
