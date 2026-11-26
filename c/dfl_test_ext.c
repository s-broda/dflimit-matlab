/* dfl_test_ext.c: the extensions of the limit of the general engine (dfl_det_ext.c, the volatility transform of
 * dfl_det_mgf.c and its paths in dfl_det_path.c) against the Julia engine (reference values of dev/ext_ref.jl of
 * DFLimit.jl in test/ext_ref.txt, first argument; the tags in the header of that file): the volatility transform, its
 * roots and t_*, the distribution functions, densities and quantiles with volatility regimes, with a normal component and
 * of the first seasonal t-ratio, the argument errors; and the limiting cases (one regime is the engine itself, a dummy
 * regime boundary changes nothing, the transform at a constant variance is the engine's, the transfer matrices of the
 * regimes without deterministic terms, the density against the distribution function, the checked computation, the
 * reflection and the endpoints of rho, the residues of the seasonal t-ratio against the doubled transform).  Each value line
 * carries its tolerance; the program reports the number of comparisons and the largest deviation of each group, and fails
 * if any value is outside its tolerance. */
#include "dflimit.h"
#include "dfl_det_int.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct { const char *name; int n, fail; double worst, worst_rel_tol; } group_t;
static group_t G[64]; static int NG = 0; static char GNAMES[64][48];
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
static void checkerr(const char *gname, int rc, const char *err, const char *frag, const char *what)
{
    int ok = rc != 0 && strstr(err, frag) != NULL;
    group_t *g = grp(gname); g->n++;
    if (!ok) { g->fail++; printf("FAIL %-10s %s: rc %d, \"%s\" without \"%s\"\n", gname, what, rc, rc ? err : "", frag); }
    else if (VERB) printf("ok   %-10s %s: %s\n", gname, what, err);
}

#define MAXD 16
static dfl_design DS[MAXD];
typedef struct { int nb; double b[8], om[9]; } volspec;
static volspec VS[MAXD];
static double num(char **p)
{
    while (**p == ' ') (*p)++;
    if (!strncmp(*p, "NaN", 3)) { *p += 3; return NAN; }
    if (!strncmp(*p, "Inf", 3)) { *p += 3; return INFINITY; }
    if (!strncmp(*p, "-Inf", 4)) { *p += 4; return -INFINITY; }
    return strtod(*p, p);
}
static int inum(char **p) { return (int) num(p); }
static int word(char **p, char *w, int len) { int k = 0; while (**p == ' ') (*p)++; while (**p && **p != ' ' && **p != '\n' && k < len - 1) w[k++] = *(*p)++; w[k] = 0; return k; }
static int read_design(char *p, int id)
{
    char err[512], kind[32]; word(&p, kind, sizeof kind);
    int rc = 1;
    if (!strcmp(kind, "trend")) rc = dfl_design_trend(inum(&p), &DS[id], err, sizeof err);
    else if (!strcmp(kind, "perron")) {
        char m[4]; word(&p, m, sizeof m); int q = inum(&p); double lam[8];
        for (int i = 0; i < q; i++) lam[i] = num(&p);
        rc = dfl_design_perron(m[0], q, lam, &DS[id], err, sizeof err);
    }
    if (rc) printf("FAIL design %d: %s\n", id, err);
    return rc;
}
/* the volatility path vid with the c of its regimes (nc values from p) into v */
static int read_dvc(char **p, int *did, dfl_vol *v, double *cs)
{
    *did = inum(p); int vid = inum(p), nc = inum(p);
    for (int i = 0; i < nc; i++) cs[i] = num(p);
    v->nreg = VS[vid].nb + 1; v->breaks = VS[vid].b; v->omega = VS[vid].om; v->c = cs;
    return nc;
}
static int vecs(char **p, double *x) { int n = inum(p); for (int i = 0; i < n; i++) x[i] = num(p); return n; }
static int cmpc(const void *a, const void *b)
{
    const cdouble *x = a, *y = b;
    if (creal(*x) != creal(*y)) return creal(*x) < creal(*y) ? -1 : 1;
    return (cimag(*x) > cimag(*y)) - (cimag(*x) < cimag(*y));
}

/* the determinant without deterministic terms by the transfer matrices of the regimes (the backward equation):
 * D = e_1' T_1 ... T_m (1, s)' */
static cdouble transfer_d0(const dfl_vol *v, cdouble s, cdouble t)
{
    int m = v->nreg; double lam[10], om[9], w = 0;
    lam[0] = 0; for (int i = 0; i < m - 1; i++) lam[i + 1] = v->breaks[i]; lam[m] = 1;
    for (int i = 0; i < m; i++) { om[i] = v->omega[i]; w += om[i] * om[i] * (lam[i + 1] - lam[i]); }
    for (int i = 0; i < m; i++) om[i] /= sqrt(w);
    cdouble x1 = DONE, x2 = s;
    for (int i = m - 1; i >= 0; i--) {
        double l = lam[i + 1] - lam[i], c = v->c[i], w2 = om[i] * om[i];
        cdouble th = csq(CPX(c * c - 2 * creal(t) * w2, -2 * cimag(t) * w2)), tl = cscale(l, th);
        cdouble ch = cmulr(cadd(cexp(tl), cexp(cneg(tl))), 0.5), sh = cfabs(tl) < 1e-6 ? cscale(l, cadd(DONE, cdivr(cmul(tl, tl), 6))) : cdiv(cmulr(csub(cexp(tl), cexp(cneg(tl))), 0.5), th);
        double e = exp(c * l);
        cdouble y1 = cscale(e, csub(cmul(csub(ch, cscale(c, sh)), x1), cmul(cscale(w2, sh), x2)));
        cdouble y2 = cscale(e, cadd(cmul(cmul(cscale(2, t), sh), x1), cmul(cadd(ch, cscale(c, sh)), x2)));
        x1 = y1; x2 = y2;
    }
    return x1;
}

static double CL1[] = {0.05};
static double CL2[] = {1.0};
static double CL3[] = {1.5};
static double CL4[] = {0.3};
static double CL5[] = {0.6};
static double CL6[] = {0, 0};
static double CL7[] = {0.5};
static double CL8[] = {0.5};
static dfl_vol mkvol(int nreg, const double *b, const double *om, const double *c) { dfl_vol v; v.nreg = nreg; v.breaks = b; v.omega = om; v.c = c; return v; }
static int SKIP_RHO = 0, SKIP_SEAS = 0;
int main(int argc, char **argv)
{
    clock_t t0 = clock();
    FILE *f = fopen(argc > 1 ? argv[1] : "test/ext_ref.txt", "r"); if (!f) { printf("no reference file\n"); return 1; }
    for (int i = 2; i < argc; i++) { if (!strcmp(argv[i], "-v")) VERB = 1; if (!strcmp(argv[i], "norho")) SKIP_RHO = 1; if (!strcmp(argv[i], "noseas")) SKIP_SEAS = 1; }
    char line[8192], err[1024], what[600];
    int nerr = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char tag[32], *p = line; word(&p, tag, sizeof tag);
        if (!strcmp(tag, "design")) { int id = inum(&p); if (read_design(p, id)) return 1; }
        else if (!strcmp(tag, "vol")) { int id = inum(&p), nb = inum(&p); VS[id].nb = nb; for (int i = 0; i < nb; i++) VS[id].b[i] = num(&p); for (int i = 0; i <= nb; i++) VS[id].om[i] = num(&p); }
        else if (!strcmp(tag, "vdval")) {
            int did; dfl_vol v; double cs[8]; int nc = read_dvc(&p, &did, &v, cs);
            double sr = num(&p), si = num(&p), tr = num(&p), ti = num(&p), dr = num(&p), di = num(&p), tol = num(&p);
            det_vmodel *VM = det_vol_model(&DS[did], &v, cs, err, sizeof err);
            snprintf(what, sizeof what, "design %d, %d regimes, c[0] %g: D(%g%+gi, %g%+gi)", did, nc, cs[0], sr, si, tr, ti);
            check("vdval", VM ? cfabs(csub(det_vol_dval(VM, CPX(sr, si), CPX(tr, ti)), CPX(dr, di))) : NAN, 0, tol, what);
            det_vol_release(VM);
        } else if (!strcmp(tag, "vroots")) {
            int did; dfl_vol v; double cs[8]; int nc = read_dvc(&p, &did, &v, cs);
            double tr = num(&p), ti = num(&p), ldr = num(&p); int r = inum(&p); cdouble ref[16], mu[16];
            for (int j = 0; j < r; j++) { double a = num(&p), b = num(&p); ref[j] = CPX(a, b); }
            double tol = num(&p);
            det_vmodel *VM = det_vol_model(&DS[did], &v, cs, err, sizeof err);
            snprintf(what, sizeof what, "design %d, %d regimes: t = %g%+gi", did, nc, tr, ti);
            if (!VM) { check("vroots", NAN, 0, tol, what); continue; }
            det_vwork W; det_vwork_init(&W, VM, NULL, 0);
            cdouble ld; int rc = det_vol_parts(&W, CPX(tr, ti), &ld, mu);
            det_vwork_free(&W);
            check("vroots", rc ? NAN : creal(ld), ldr, tol, what);
            check("vroots", VM->r, r, 0, what);
            qsort(mu, r, sizeof(cdouble), cmpc); qsort(ref, r, sizeof(cdouble), cmpc);
            double mx = 0, d = 0; for (int j = 0; j < r; j++) { mx = fmax(mx, cfabs(ref[j])); d = fmax(d, cfabs(csub(mu[j], ref[j]))); }
            check("vroots", d, 0, tol * mx, what);
            det_vol_release(VM);
        } else if (!strcmp(tag, "vtstar")) {
            int did; dfl_vol v; double cs[8]; int nc = read_dvc(&p, &did, &v, cs);
            double ts = num(&p), tol = num(&p), tsc = NAN;
            det_vmodel *VM = det_vol_model(&DS[did], &v, cs, err, sizeof err);
            if (VM && det_vol_tstar(VM, &tsc, NULL, err, sizeof err)) tsc = NAN;
            snprintf(what, sizeof what, "design %d, %d regimes, c[0] %g", did, nc, cs[0]);
            check("vtstar", tsc, ts, tol, what);
            det_vol_release(VM);
        } else if (!strcmp(tag, "vval")) {
            char fn[8]; word(&p, fn, sizeof fn);
            int did; dfl_vol v; double cs[8]; int nc = read_dvc(&p, &did, &v, cs);
            double c0 = num(&p), x[32], ref[32], got[32]; int n = vecs(&p, x); vecs(&p, ref); double tol = num(&p);
            dfl_det_opts o; dfl_det_defaults(&o); o.vol = &v; o.c0 = c0;
            int rc = dfl_det_dist(!strcmp(fn, "cdf") ? 2 : 1, &DS[did], n, x, cs[0], &o, got, err, sizeof err);
            for (int i = 0; i < n; i++) { snprintf(what, sizeof what, "%s design %d, %d regimes, c[0] %g, c0 %g at x = %g%s%s", fn, did, nc, cs[0], c0, x[i], rc ? ": " : "", rc ? err : ""); check(!strcmp(fn, "cdf") ? "vcdf" : "vpdf", rc ? NAN : got[i], ref[i], tol, what); }
        } else if (!strcmp(tag, "vquant")) {
            int did; dfl_vol v; double cs[8]; int nc = read_dvc(&p, &did, &v, cs);
            double c0 = num(&p), pp = num(&p), x = num(&p), tol = num(&p), q = NAN;
            dfl_det_opts o; dfl_det_defaults(&o); o.vol = &v; o.c0 = c0;
            int rc = dfl_det_quantile(&DS[did], 1, &pp, cs[0], -2, 1e-9, &o, &q, err, sizeof err);
            snprintf(what, sizeof what, "design %d, %d regimes, c[0] %g: p = %g%s%s", did, nc, cs[0], pp, rc ? ": " : "", rc ? err : "");
            check("vquant", rc ? NAN : q, x, tol, what);
        } else if (!strcmp(tag, "cval") && !SKIP_RHO) {
            char fn[8], meth[16]; word(&p, fn, sizeof fn);
            int did = inum(&p), st = inum(&p); double g0 = num(&p), c = num(&p), c0 = num(&p), rho = num(&p); word(&p, meth, sizeof meth);
            double x[32], ref[32], got[32]; int n = vecs(&p, x); vecs(&p, ref); double tol = num(&p);
            dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0; o.start = st; o.gamma0 = g0; o.rho = rho;
            o.method = rho == 0 ? DFL_DET_AUTO : !strcmp(meth, "lines") ? DFL_DET_LINES : DFL_DET_MIXTURE;
            int rc = dfl_det_dist(!strcmp(fn, "cdf") ? 2 : 1, &DS[did], n, x, c, &o, got, err, sizeof err);
            for (int i = 0; i < n; i++) { snprintf(what, sizeof what, "%s design %d start %d c %g c0 %g rho %g %s at x = %g%s%s", fn, did, st, c, c0, rho, meth, x[i], rc ? ": " : "", rc ? err : ""); check(!strcmp(fn, "cdf") ? "ccdf" : "cpdf", rc ? NAN : got[i], ref[i], tol, what); }
        } else if (!strcmp(tag, "cquant") && !SKIP_RHO) {
            int did = inum(&p), st = inum(&p); double g0 = num(&p), c = num(&p), c0 = num(&p), rho = num(&p), pp = num(&p), x = num(&p), tol = num(&p), q = NAN;
            dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0; o.start = st; o.gamma0 = g0; o.rho = rho;
            int rc = dfl_det_quantile(&DS[did], 1, &pp, c, -2, 1e-9, &o, &q, err, sizeof err);
            snprintf(what, sizeof what, "design %d c %g c0 %g rho %g: p = %g%s%s", did, c, c0, rho, pp, rc ? ": " : "", rc ? err : "");
            check("cquant", rc ? NAN : q, x, tol, what);
        } else if (!strcmp(tag, "sval") && !SKIP_SEAS) {
            char fn[8]; word(&p, fn, sizeof fn);
            if (strncmp(fn, "t3", 2)) continue;          /* (the second t-ratio and F are not ported) */
            int did = inum(&p); double c = num(&p), c0 = num(&p); int vid = inum(&p);
            double x[32], ref[32], got[32]; int n = vecs(&p, x); vecs(&p, ref); double tol = num(&p);
            dfl_det_opts o; dfl_det_defaults(&o); o.c0 = c0;
            dfl_vol v; double cs[2] = {c, c};
            if (vid) { v.nreg = VS[vid].nb + 1; v.breaks = VS[vid].b; v.omega = VS[vid].om; v.c = cs; o.vol = &v; }
            int rc = dfl_hegy_t3(!strcmp(fn, "t3cdf") ? 2 : 1, &DS[did], n, x, c, &o, got, err, sizeof err);
            for (int i = 0; i < n; i++) { snprintf(what, sizeof what, "%s design %d c %g c0 %g vol %d at x = %g%s%s", fn, did, c, c0, vid, x[i], rc ? ": " : "", rc ? err : ""); check(fn, rc ? NAN : got[i], ref[i], tol, what); }
        } else if (!strcmp(tag, "squant") && !SKIP_SEAS) {
            char fn[8]; word(&p, fn, sizeof fn);
            if (strcmp(fn, "t3")) continue;
            int did = inum(&p); double c = num(&p), pp = num(&p), x = num(&p), tol = num(&p), q = NAN;
            int rc = dfl_hegy_t3_quantile(&DS[did], 1, &pp, c, NAN, NAN, NULL, &q, err, sizeof err);
            snprintf(what, sizeof what, "t3 design %d c %g: p = %g%s%s", did, c, pp, rc ? ": " : "", rc ? err : "");
            check("t3quant", rc ? NAN : q, x, tol, what);
        } else if (!strcmp(tag, "err")) nerr++;
    }
    fclose(f);
    /* ---------------------------------------------------------------------------------------------------------------- */
    /* the argument errors of the reference file (its err lines; the fragments of the messages), as the calls in C */
    {
        double x = -2, cs3[3] = {0, -1, 0}, b5[1] = {0.5}, om11[1] = {1}, om1m[2] = {1, -1}, b12[1] = {1.2}, b53[2] = {0.5, 0.3}, b3[1] = {0.3}, om13[2] = {1, 3}, out[2];
        dfl_vol v; dfl_det_opts o; int rc;
        v = mkvol(2, b12, NULL, NULL); rc = dfl_vol_check(&v, err, sizeof err); checkerr("errors", rc, err, "must lie in (0, 1)", "volpath: a break date at 1.2");
        v = mkvol(3, b53, NULL, NULL); rc = dfl_vol_check(&v, err, sizeof err); checkerr("errors", rc, err, "strictly increasing", "volpath: break dates 0.5, 0.3");
        v = mkvol(2, b5, om1m, NULL); rc = dfl_vol_check(&v, err, sizeof err); checkerr("errors", rc, err, "finite and positive", "volpath: omega -1");
        (void) om11;
        dfl_vol V2 = {2, b3, om13, NULL};
        dfl_det_defaults(&o); o.vol = &V2; rc = dfl_det_dist(2, &DS[2], 1, &x, 2.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "must be at most 0", "dist_cdf vol, c = 2");
        dfl_det_defaults(&o); o.vol = &V2; o.start = DFL_START_STAT; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "zero start only", "dist_cdf vol, stationary start");
        dfl_det_defaults(&o); o.vol = &V2; o.start = DFL_START_FIXED; o.gamma0 = 1; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "zero start only", "dist_cdf vol, fixed start");
        dfl_det_defaults(&o); o.vol = &V2; o.gls = -13.5; rc = dfl_det_dist(2, &DS[3], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "GLS detrending", "dist_cdf vol, gls = -13.5");
        dfl_design F1; dfl_design_detrend_first(&DS[3], &F1, err, sizeof err);
        dfl_det_defaults(&o); o.vol = &V2; rc = dfl_det_dist(2, &F1, 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "detrend_first", "dist_cdf vol, detrend_first");
        dfl_design_free(&F1);
        dfl_det_defaults(&o); o.vol = &V2; o.method = DFL_DET_LINE; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "method must be 'auto' or 'paths'", "dist_cdf vol, method line");
        dfl_det_defaults(&o); o.vol = &V2; rc = dfl_det_dist(2, &DS[2], 1, &x, -500.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "at least -480", "dist_cdf vol, c = -500");
        dfl_vol V2c = {2, b3, om13, cs3}; cs3[1] = 1;
        dfl_det_defaults(&o); o.vol = &V2c; rc = dfl_det_dist(2, &DS[2], 1, &x, 0.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "at most 0", "dist_cdf vol, c = [0, 1]");
        dfl_det_defaults(&o); o.vol = &V2; o.rho = 0.5; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "rho and vol together", "dist_cdf vol and rho = 0.5");
        dfl_det_defaults(&o); o.vol = &V2; o.rho = 0.5; rc = dfl_det_quantile(&DS[2], 1, CL1, -5.0, NAN, NAN, &o, out, err, sizeof err); checkerr("errors", rc, err, "rho and vol together", "dist_quantile vol and rho = 0.5");
        if (!SKIP_RHO) {
            dfl_det_defaults(&o); o.rho = 1.5; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "rho must lie in [-1, 1]", "dist_cdf rho = 1.5");
            dfl_det_defaults(&o); o.rho = NAN; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "rho must lie in [-1, 1]", "dist_cdf rho = NaN");
            dfl_det_defaults(&o); o.rho = 0.5; o.gls = -13.5; rc = dfl_det_dist(2, &DS[3], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "GLS detrending", "dist_cdf rho = 0.5, gls = -13.5");
            dfl_det_defaults(&o); o.rho = 0.5; o.start = DFL_START_FIXED; o.gamma0 = 1; o.method = DFL_DET_LINES; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "zero or a stationary start", "dist_cdf rho = 0.5, method lines, fixed start");
            dfl_det_defaults(&o); o.rho = 0.5; o.method = DFL_DET_LINES; rc = dfl_det_dist(2, &DS[2], 1, &x, 1.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "c <= 0", "dist_cdf rho = 0.5, method lines, c = 1");
            dfl_det_defaults(&o); o.rho = 0.5; o.method = DFL_DET_PATHS; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "method must be", "dist_cdf rho = 0.5, method paths");
            dfl_det_defaults(&o); o.rho = 0.5; o.hfac = 2; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "do not apply", "dist_cdf rho = 0.5, hfac = 2");
            dfl_det_defaults(&o); o.rho = 0.5; o.n = 16; rc = dfl_det_dist(2, &DS[2], 1, &x, -5.0, &o, out, err, sizeof err); checkerr("errors", rc, err, "do not apply", "dist_cdf rho = 0.5, n = 16");
            double c0;
            rc = dfl_cadf_c0(-5.0, 0.0, NAN, &c0, err, sizeof err); checkerr("errors", rc, err, "not be 0", "cadf_c0 rho = 0");
            rc = dfl_predictive_c0(-5.0, 1.5, 0.0, &c0, err, sizeof err); checkerr("errors", rc, err, "not be 0", "predictive_c0 rho = 1.5");
        }
        if (!SKIP_SEAS) {
            rc = dfl_hegy_t3(2, &DS[2], 1, CL2, 1.0, NULL, out, err, sizeof err); checkerr("errors", rc, err, "must be at most 0", "hegy_t3_cdf c = 1");
            rc = dfl_hegy_t3_quantile(&DS[2], 1, CL3, 0.0, NAN, NAN, NULL, out, err, sizeof err); checkerr("errors", rc, err, "p must lie in (0, 1)", "hegy_t3_quantile p = 1.5");
        }
    }
    /* ---------------------------------------------------------------------------------------------------------------- */
    /* the limiting cases of the volatility transform */
    {
        double b4[1] = {0.4}, om11[2] = {1, 1}, b3[1] = {0.3}, om13[2] = {1, 3}, b27[2] = {0.2, 0.7}, om3[3] = {2, 1, 0.5}, c3[3] = {-3, 0, -10};
        double b5[4] = {0.1, 0.3, 0.6, 0.8}, om5[5] = {1, 2, 0.5, 3, 1.5}, c5[5] = {-1, 0, -4, -2, 0};
        cdouble pts[5][2] = {{{0.3, 0}, {-2, 0}}, {{-0.4, 0}, {0.5, 0}}, {{0.2, 0.7}, {0.3, -4}}, {{0.1, 0}, {-50, 20}}, {{0, 0}, {3, 0}}};
        dfl_design Pa, Pc; dfl_design_perron('A', 1, CL4, &Pa, err, sizeof err); dfl_design_perron('C', 1, CL5, &Pc, err, sizeof err);
        const dfl_design *Ds[5] = {&DS[1], &DS[2], &DS[3], &Pa, &Pc};
        /* the transform at a constant variance is the engine's */
        for (int d = 0; d < 5; d++) for (int ci = 0; ci < 2; ci++) {
            double c = ci ? -5.0 : 0.0, cc[2] = {c, c};
            dfl_vol v = {2, b4, om11, cc};
            det_vmodel *VM = det_vol_model(Ds[d], &v, cc, err, sizeof err);
            det_prep *P = det_prep_cached(Ds[d], NAN, err, sizeof err);
            det_model Mo; det_model_make(P, 0, 0, &Mo, err, sizeof err);
            det_work W; det_work_init(&W, &Mo, NULL, 0);
            for (int k = 0; k < 4; k++) {
                cdouble a = det_vol_dval(VM, pts[k][0], pts[k][1]), b = det_dval_st(&W, c, pts[k][0], pts[k][1]);
                snprintf(what, sizeof what, "design %d at a constant variance, c %g: D(%g%+gi, %g%+gi) against the engine's", d, c, creal(pts[k][0]), cimag(pts[k][0]), creal(pts[k][1]), cimag(pts[k][1]));
                check("vconst", cfabs(csub(a, b)) / cfabs(b), 0, 1e-10, what);
            }
            det_work_free(&W); det_model_free(&Mo); det_prep_release(P); det_vol_release(VM);
        }
        /* no deterministic terms: the product of the transfer matrices of the regimes */
        dfl_vol vv[3] = {{2, b3, om13, CL6}, {3, b27, om3, c3}, {5, b5, om5, c5}};
        for (int i = 0; i < 3; i++) {
            det_vmodel *VM = det_vol_model(&DS[1], &vv[i], vv[i].c, err, sizeof err);
            for (int k = 0; k < 5; k++) {
                cdouble a = det_vol_dval(VM, pts[k][0], pts[k][1]), b = transfer_d0(&vv[i], pts[k][0], pts[k][1]);
                snprintf(what, sizeof what, "trend(-1), %d regimes: D(%g%+gi, %g%+gi) against the transfer matrices", vv[i].nreg, creal(pts[k][0]), cimag(pts[k][0]), creal(pts[k][1]), cimag(pts[k][1]));
                check("vtransfer", cfabs(csub(a, b)) / cfabs(b), 0, 1e-11, what);
            }
            det_vol_release(VM);
        }
        /* the reversal at the unit root (the values of the Julia tests) */
        {
            double c2[2] = {0, 0}, b7[1] = {0.7}, om31[2] = {3, 1};
            dfl_vol v = {2, b3, om13, c2}, vr = {2, b7, om31, c2};
            det_vmodel *VM = det_vol_model(&DS[2], &v, c2, err, sizeof err), *VR = det_vol_model(&DS[2], &vr, c2, err, sizeof err), *VT = det_vol_model(&DS[3], &v, c2, err, sizeof err);
            check("vreverse", creal(det_vol_dval(VM, CPX(0.5, 0), CPX(-2, 0))), 1.6815863197, 1e-9, "trend(0), volpath([0.3], [1, 3]): D(0.5, -2)");
            check("vreverse", creal(det_vol_dval(VM, CPX(-0.5, 0), CPX(-2, 0))), 2.0691752124, 1e-9, "D(-0.5, -2)");
            check("vreverse", creal(det_vol_dval(VR, CPX(-0.5, 0), CPX(-2, 0))), 1.6815863197, 1e-9, "the reversed path: D(-0.5, -2)");
            check("vreverse", creal(det_vol_dval(VM, CPX(1, 0), CPX(1, 0))), 0.4215903315, 1e-9, "D(1, 1)");
            check("vreverse", creal(det_vol_dval(VT, CPX(0.5, 0), CPX(-2, 0))), 1.2374269682, 1e-9, "trend(1): D(0.5, -2)");
            check("vreverse", creal(det_vol_dval(VT, CPX(-0.5, 0), CPX(-2, 0))), 1.3528083151, 1e-9, "trend(1): D(-0.5, -2)");
            det_vol_release(VM); det_vol_release(VR); det_vol_release(VT);
        }
        /* t_* at a constant variance is the engine's (model C at 1/2 a double zero) */
        {
            dfl_design C5; dfl_design_perron('C', 1, CL7, &C5, err, sizeof err);
            const dfl_design *Dt[4] = {&DS[2], &DS[3], &C5, &DS[1]}; double ct[4] = {0, -5, 0, -100};
            for (int i = 0; i < 4; i++) {
                double cc[2] = {ct[i], ct[i]}; dfl_vol v = {2, b4, om11, cc};
                det_vmodel *VM = det_vol_model(Dt[i], &v, cc, err, sizeof err);
                det_prep *P = det_prep_cached(Dt[i], NAN, err, sizeof err); det_model Mo; det_model_make(P, 0, 0, &Mo, err, sizeof err);
                double ts = NAN, te = NAN; int m = 0, me = 0;
                det_vol_tstar(VM, &ts, &m, err, sizeof err); det_tstar(&Mo, ct[i], &te, &me, err, sizeof err);
                snprintf(what, sizeof what, "t_* at a constant variance, design %d, c %g, against the engine's", i, ct[i]);
                check("vtstar1", fabs(ts - te) / te, 0, me == 2 ? 1e-6 : 1e-10, what);
                check("vtstar1", m, me, 0, what);
                det_model_free(&Mo); det_prep_release(P); det_vol_release(VM);
            }
            dfl_design_free(&C5);
        }
        /* one regime is the engine itself; a dummy boundary at a constant variance is the engine's paths to rounding */
        {
            double xs[4] = {-4, -2.86, -1, 0.5}, a[4], b[4];
            dfl_design Pb; dfl_design_perron('B', 1, CL8, &Pb, err, sizeof err);
            const dfl_design *D1[2] = {&DS[2], &Pb};
            for (int d = 0; d < 2; d++) for (int ci = 0; ci < 2; ci++) for (int kind = 1; kind <= 2; kind++) {
                double c = ci ? -10.0 : 0.0; dfl_vol v1 = {1, NULL, NULL, NULL};
                dfl_det_opts o; dfl_det_defaults(&o);
                int rc = dfl_det_dist(kind, D1[d], 4, xs, c, &o, b, err, sizeof err);
                o.vol = &v1; rc |= dfl_det_dist(kind, D1[d], 4, xs, c, &o, a, err, sizeof err);
                for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "one regime, design %d, c %g, kind %d, x %g", d, c, kind, xs[i]); check("vone", rc ? NAN : a[i], b[i], 1e-15, what); }
            }
            dfl_design_free(&Pb);
            double pq = 0.05, q1 = NAN, q2 = NAN; dfl_vol v1 = {1, NULL, NULL, NULL}; dfl_det_opts o; dfl_det_defaults(&o);
            dfl_det_quantile(&DS[3], 1, &pq, 0.0, NAN, NAN, &o, &q2, err, sizeof err); o.vol = &v1; dfl_det_quantile(&DS[3], 1, &pq, 0.0, NAN, NAN, &o, &q1, err, sizeof err);
            check("vone", q1, q2, 1e-15, "one regime: the quantile of trend(1)");
            const dfl_design *Dd[4] = {&DS[1], &DS[2], &DS[3], &Pc}; double cd[3] = {0, -10, -100};
            for (int d = 0; d < 4; d++) for (int ci = 0; ci < 3; ci++) for (int kind = 1; kind <= 2; kind++) {
                double cc[2] = {cd[ci], cd[ci]}; dfl_vol v = {2, b4, om11, cc};
                dfl_det_opts o2; dfl_det_defaults(&o2); o2.method = DFL_DET_PATHS;
                int rc = dfl_det_dist(kind, Dd[d], 4, xs, cd[ci], &o2, b, err, sizeof err);
                o2.vol = &v; rc |= dfl_det_dist(kind, Dd[d], 4, xs, cd[ci], &o2, a, err, sizeof err);
                for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "a dummy boundary at 0.4, design %d, c %g, kind %d, x %g%s%s", d, cd[ci], kind, xs[i], rc ? ": " : "", rc ? err : ""); check("vdummy", rc ? NAN : a[i], b[i], 1e-12, what); }
            }
            double b25[2] = {0.25, 0.75}, om111[3] = {1, 1, 1}, c555[3] = {-5, -5, -5}; dfl_vol v3 = {3, b25, om111, c555};
            dfl_det_opts o3; dfl_det_defaults(&o3); o3.method = DFL_DET_PATHS; o3.c0 = -2;
            int rc = dfl_det_dist(2, &DS[2], 4, xs, -5.0, &o3, b, err, sizeof err); o3.vol = &v3; rc |= dfl_det_dist(2, &DS[2], 4, xs, -5.0, &o3, a, err, sizeof err);
            for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "two dummy boundaries, c -5, c0 -2, x %g", xs[i]); check("vdummy", rc ? NAN : a[i], b[i], 1e-12, what); }
        }
        /* the checked computation, the density against the distribution function, non-finite points */
        {
            double b1[1] = {0.1}, om13b[2] = {1, 3}, x3[3] = {-3, -2, -0.5}, a[3], b[3];
            dfl_vol v = {2, b1, om13b, NULL};
            for (int kind = 1; kind <= 2; kind++) {
                dfl_det_opts o; dfl_det_defaults(&o); o.vol = &v;
                int rc = dfl_det_dist(kind, &DS[2], 3, x3, -5.0, &o, a, err, sizeof err); o.verified = 1; rc |= dfl_det_dist(kind, &DS[2], 3, x3, -5.0, &o, b, err, sizeof err);
                for (int i = 0; i < 3; i++) { snprintf(what, sizeof what, "verified, kind %d, x %g", kind, x3[i]); check("vverified", rc ? NAN : a[i], b[i], 1e-9, what); }
            }
            double h = 1e-3, b7[1] = {0.7}, c00[2] = {0, 0};
            dfl_vol vd[2] = {{2, b7, om13b, c00}, {3, b27, om3, c3}}; const dfl_design *Dh[2] = {&DS[1], &DS[3]};
            for (int i = 0; i < 2; i++) {
                double xx[2] = {-2 - h, -2 + h}, F[2], fd, xm = -2;
                dfl_det_opts o; dfl_det_defaults(&o); o.vol = &vd[i];
                int rc = dfl_det_dist(2, Dh[i], 2, xx, 0.0, &o, F, err, sizeof err) | dfl_det_dist(1, Dh[i], 1, &xm, 0.0, &o, &fd, err, sizeof err);
                snprintf(what, sizeof what, "the density against the difference quotient, case %d", i);
                check("vdensity", rc ? NAN : (F[1] - F[0]) / (2 * h), fd, 1e-6, what);
            }
            double xn[3] = {-INFINITY, INFINITY, NAN}, vn[3];
            dfl_det_opts o; dfl_det_defaults(&o); o.vol = &v;
            int rc = dfl_det_dist(2, &DS[2], 3, xn, -5.0, &o, vn, err, sizeof err);
            check("vnonfinite", rc ? NAN : vn[0], 0, 0, "-Inf"); check("vnonfinite", rc ? NAN : vn[1], 1, 0, "Inf"); check("vnonfinite", rc ? 0 : isnan(vn[2]) ? 1 : 0, 1, 0, "NaN");
        }
        dfl_design_free(&Pa); dfl_design_free(&Pc);
    }
    if (!SKIP_RHO) {
        /* the endpoints and the reflection of rho */
        double xs[4] = {-3, -1, 0.5, 2}, a[4], b[4];
        dfl_det_opts o; dfl_det_defaults(&o); o.rho = 0;
        int rc = dfl_det_dist(2, &DS[2], 4, xs, -5.0, &o, a, err, sizeof err);
        for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "rho = 0: Phi(%g)", xs[i]); check("rhoends", rc ? NAN : a[i], 0.5 * erfc(-xs[i] / sqrt(2)), 1e-15, what); }
        dfl_det_defaults(&o); o.rho = -1; rc = dfl_det_dist(2, &DS[2], 4, xs, -5.0, &o, a, err, sizeof err);
        double mx[4]; for (int i = 0; i < 4; i++) mx[i] = -xs[i];
        dfl_det_defaults(&o); rc |= dfl_det_dist(2, &DS[2], 4, mx, -5.0, &o, b, err, sizeof err);
        for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "rho = -1: 1 - F(-x) at x = %g", xs[i]); check("rhoends", rc ? NAN : a[i], 1 - b[i], 1e-15, what); }
        for (int m = 0; m < 2; m++) {
            dfl_det_defaults(&o); o.rho = -0.7; o.method = m ? DFL_DET_MIXTURE : DFL_DET_LINES; rc = dfl_det_dist(2, &DS[2], 4, xs, -5.0, &o, a, err, sizeof err);
            o.rho = 0.7; rc |= dfl_det_dist(2, &DS[2], 4, mx, -5.0, &o, b, err, sizeof err);
            for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "rho = -0.7 by the reflection (%s) at x = %g", m ? "mixture" : "lines", xs[i]); check("rhoends", rc ? NAN : a[i], 1 - b[i], 1e-15, what); }
        }
        /* the two routes against each other (verified), and the density against the distribution function */
        for (int m = 0; m < 2; m++) {
            dfl_det_defaults(&o); o.rho = sqrt(0.5); o.method = m ? DFL_DET_MIXTURE : DFL_DET_LINES;
            rc = dfl_det_dist(2, &DS[3], 4, xs, -3.0, &o, a, err, sizeof err); o.verified = 1; rc |= dfl_det_dist(2, &DS[3], 4, xs, -3.0, &o, b, err, sizeof err);
            for (int i = 0; i < 4; i++) { snprintf(what, sizeof what, "rho^2 = 0.5, %s against verified at x = %g", m ? "mixture" : "lines", xs[i]); check("rhoverified", rc ? NAN : a[i], b[i], 1e-9, what); }
            double h = 1e-3, xx[2] = {-1 - h, -1 + h}, F[2], fd, xm = -1;
            o.verified = 0;
            rc = dfl_det_dist(2, &DS[3], 2, xx, -3.0, &o, F, err, sizeof err) | dfl_det_dist(1, &DS[3], 1, &xm, -3.0, &o, &fd, err, sizeof err);
            snprintf(what, sizeof what, "the density against the difference quotient (%s)", m ? "mixture" : "lines");
            check("rhodensity", rc ? NAN : (F[1] - F[0]) / (2 * h), fd, 1e-6, what);
        }
        double c0;
        rc = dfl_cadf_c0(-5.0, sqrt(0.5), NAN, &c0, err, sizeof err); check("rhohelpers", rc ? NAN : c0, 5.0, 1e-13, "cadf_c0(-5, sqrt(0.5)) = 5");
        rc = dfl_cadf_c0(0.0, 0.3, 2.0, &c0, err, sizeof err); check("rhohelpers", rc ? NAN : c0, 0.0, 0, "cadf_c0(0, ...) = 0");
        rc = dfl_predictive_c0(-5.0, -0.95, 1.0, &c0, err, sizeof err); check("rhohelpers", rc ? NAN : c0, -5.0 + 1 / 0.95, 1e-14, "predictive_c0(-5, -0.95, 1)");
    }
    if (!SKIP_SEAS) {
        /* the residues against the doubled transform (verified), the density against the distribution function */
        double xs[3] = {-3, -1.5, 0.2}, a[3], b[3];
        dfl_det_opts o; dfl_det_defaults(&o);
        int rc = dfl_hegy_t3(2, &DS[3], 3, xs, -2.0, &o, a, err, sizeof err); o.verified = 1; rc |= dfl_hegy_t3(2, &DS[3], 3, xs, -2.0, &o, b, err, sizeof err);
        for (int i = 0; i < 3; i++) { snprintf(what, sizeof what, "t3 trend(1), c -2: the residues against verified at x = %g", xs[i]); check("t3verified", rc ? NAN : a[i], b[i], 1e-9, what); }
        double h = 1e-3, xx[2] = {-2 - h, -2 + h}, F[2], fd, xm = -2;
        dfl_det_defaults(&o);
        rc = dfl_hegy_t3(2, &DS[2], 2, xx, 0.0, &o, F, err, sizeof err) | dfl_hegy_t3(1, &DS[2], 1, &xm, 0.0, &o, &fd, err, sizeof err);
        check("t3density", rc ? NAN : (F[1] - F[0]) / (2 * h), fd, 1e-6, "t3 trend(0): the density against the difference quotient");
        double qq = NAN, pp = 0.05;
        rc = dfl_hegy_t3_quantile(&DS[2], 1, &pp, 0.0, NAN, NAN, NULL, &qq, err, sizeof err);
        double Fq = NAN; if (!rc) rc = dfl_hegy_t3(2, &DS[2], 1, &qq, 0.0, NULL, &Fq, err, sizeof err);
        check("t3density", rc ? NAN : Fq, 0.05, 1e-10, "t3 trend(0): F at the 5% quantile");
    }
    printf("(%d err lines in the reference file; their calls in C are in the group errors)\n", nerr);
    int fails = 0, n = 0;
    printf("%-11s %8s %8s %12s %10s\n", "group", "values", "failed", "largest dev", "dev/tol");
    for (int i = 0; i < NG; i++) { printf("%-11s %8d %8d %12.2e %10.2e\n", G[i].name, G[i].n, G[i].fail, G[i].worst, G[i].worst_rel_tol); fails += G[i].fail; n += G[i].n; }
    printf("%d comparisons, %d failed, %.1f s cpu (%d threads)\n", n, fails, (double) (clock() - t0) / CLOCKS_PER_SEC, dfl_nworkers());
    return fails ? 1 : 0;
}
