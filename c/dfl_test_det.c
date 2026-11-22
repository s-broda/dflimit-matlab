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
            dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0; o.start = st; o.gamma0 = g0; o.gls = cb; o.verified = ver;
            double t1 = dfl_wtime(); int rc = 1;
            if (!strcmp(fn, "cdf")) rc = dfl_det_dist(2, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "pdf")) rc = dfl_det_dist(1, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "quant")) rc = dfl_det_quantile(&DS[id], n, x, c, par, 1e-9, &o, got, err, sizeof err);
            else if (!strcmp(fn, "coefcdf")) rc = dfl_det_coef(2, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "coefpdf")) rc = dfl_det_coef(1, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "coefq")) rc = dfl_det_coef_quantile(&DS[id], n, x, c, par, 1e-9, &o, got, err, sizeof err);
            else if (!strcmp(fn, "coefgp")) rc = dfl_det_coef_cdf_gp(&DS[id], n, x, c, &o, 60000, 4000, got, err, sizeof err);
            else if (!strcmp(fn, "mom")) rc = dfl_det_moments(&DS[id], c, n, &o, got, NULL, err, sizeof err);
            else if (!strcmp(fn, "lpow")) { o.c0 = 0; rc = dfl_det_local_power(&DS[id], n, x, 0.05, (int) par, &o, got, err, sizeof err); }
            else if (!strcmp(fn, "lincdf")) rc = dfl_det_lin(2, &DS[id], n, x, par, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "linpdf")) rc = dfl_det_lin(1, &DS[id], n, x, par, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "linq")) rc = dfl_det_lin_quantile(&DS[id], x[0], par, c, 0.0, 1e-10, &o, got, err, sizeof err);
            else if (!strcmp(fn, "adfcdf")) rc = dfl_det_adfmax(2, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "adfpdf")) rc = dfl_det_adfmax(1, &DS[id], n, x, c, &o, got, err, sizeof err);
            else if (!strcmp(fn, "adfq")) rc = dfl_det_adfmax_quantile(&DS[id], x[0], c, par, 1e-10, &o, got, err, sizeof err);
            else if (!strcmp(fn, "env")) rc = dfl_det_power_envelope(&DS[id], n, x, par, got, err, sizeof err);
            else if (!strcmp(fn, "ptcdf")) rc = dfl_det_pt(2, &DS[id], n, x, par, c, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "ptpdf")) rc = dfl_det_pt(1, &DS[id], n, x, par, c, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "ptq")) rc = dfl_det_pt_quantile(&DS[id], x[0], par, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "ptpow")) rc = dfl_det_pt_power(&DS[id], par, c, x[0], 0, 0, &got[0], &got[1], err, sizeof err);
            else if (!strcmp(fn, "ptenv")) rc = dfl_det_pt_envelope(&DS[id], n, x, par, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "lmcdf")) rc = dfl_det_lm_cdf(&DS[id], n, x, c, got, err, sizeof err);
            else if (!strcmp(fn, "lmq")) rc = dfl_det_lm_quantile(&DS[id], x[0], c, got, err, sizeof err);
            else if (!strcmp(fn, "lmpow")) rc = dfl_det_lm_power(&DS[id], n, x, par, &got[0], got + 1, err, sizeof err);
            else if (!strcmp(fn, "mzacdf")) rc = dfl_det_mza_cdf(&DS[id], n, x, c, par, NULL, got, err, sizeof err);
            else if (!strcmp(fn, "mzaq")) rc = dfl_det_mza_quantile(&DS[id], x[0], c, par, NAN, 1e-9, NULL, got, err, sizeof err);
            else if (!strcmp(fn, "msbcdf")) rc = dfl_det_msb_cdf(&DS[id], n, x, c, par, got, err, sizeof err);
            else if (!strcmp(fn, "msbq")) rc = dfl_det_msb_quantile(&DS[id], x[0], c, par, got, err, sizeof err);
            else { printf("unknown value function %s\n", fn); }
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
        if (!strcmp(tag, "route")) {
            char fn[32], sp[16]; int k2; sscanf(p, "%31s%n", fn, &k2); p += k2;
            int e = inum(&p); sscanf(p, "%15s%n", sp, &k2); p += k2;
            int a = inum(&p), st = inum(&p); double g0 = num(&p), gl = num(&p), c = num(&p), c0 = num(&p), par = num(&p); int n = inum(&p);
            double *x = malloc(n * sizeof(double)), *want = malloc(64 * sizeof(double)), *got = malloc(64 * sizeof(double));
            for (int i = 0; i < n; i++) x[i] = num(&p);
            int nw = inum(&p); for (int i = 0; i < nw; i++) want[i] = num(&p);
            double tol = num(&p);
            dfl_spec S; S.D = strcmp(sp, "design") ? NULL : &DS[a]; S.model = S.D ? 0 : a; S.start = st; S.g0 = g0; S.gls = S.D ? gl : NAN;
            int rc = 1;
            if (!strcmp(fn, "cdf")) rc = dfl_engine_dist(2, n, x, &S, c, c0, e, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "pdf")) rc = dfl_engine_dist(1, n, x, &S, c, c0, e, 0, 0, got, err, sizeof err);
            else if (!strcmp(fn, "icdf")) rc = dfl_engine_icdf(n, x, &S, c, c0, e, 0, 1e-10, got, err, sizeof err);
            else if (!strcmp(fn, "mom")) rc = dfl_engine_mom((int) x[0], &S, c, c0, e, got, err, sizeof err);
            else if (!strcmp(fn, "coef")) rc = dfl_engine_coef(n, x, &S, c, c0, e, got, err, sizeof err);
            else if (!strcmp(fn, "coefq")) rc = dfl_engine_coef_icdf(n, x, &S, c, c0, e, 1e-9, got, err, sizeof err);
            else if (!strcmp(fn, "lin")) rc = dfl_engine_lincdf(n, x, par, a, c, gl, e, got, err, sizeof err);
            else if (!strcmp(fn, "env")) rc = dfl_engine_envelope(n, x, a, par, e, got, err, sizeof err);
            else if (!strcmp(fn, "adf")) rc = dfl_engine_adfmax(2, n, x, &S, c, e, got, err, sizeof err);
            else if (!strcmp(fn, "adfpdf")) rc = dfl_engine_adfmax(1, n, x, &S, c, e, got, err, sizeof err);
            else if (!strcmp(fn, "region1")) { double pr[2] = {par, x[0]}; rc = dfl_engine_region(1, pr, &S, c, e, got, err, sizeof err); }
            else if (!strcmp(fn, "region2")) { double pr[2] = {x[0], 0.0}; rc = dfl_engine_region(2, pr, &S, c, e, got, err, sizeof err); }
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
    int fails = 0, n = 0;
    printf("%-10s %8s %8s %12s %10s\n", "group", "values", "failed", "largest dev", "dev/tol");
    for (int i = 0; i < NG; i++) { printf("%-10s %8d %8d %12.2e %10.2e\n", G[i].name, G[i].n, G[i].fail, G[i].worst, G[i].worst_rel_tol); fails += G[i].fail; n += G[i].n; }
    printf("%d comparisons, %d failed, %.1f s (%d threads)\n", n, fails, dfl_wtime() - t0, dfl_nworkers());
    return fails ? 1 : 0;
}
