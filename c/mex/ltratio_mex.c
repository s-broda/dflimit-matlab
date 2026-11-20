/* ltratio_mex.c: MATLAB gateway to the C engine (c/), one MEX file for every operation of the MATLAB programs; the
 * m-files dfx_core.m, tr_limdist.m, tr_contour.m and tr_limmom.m next to it present the engine under the names and signatures of
 * the MATLAB engines, so that ltratio_cdf, ltratio_pdf, ltratio_icdf, ltratio_mom and ltratio_test run unchanged on
 * top of it, and ltratio.m is a self-contained wrapper.  All arguments are double arrays (logical flags as 0/1):
 *
 *   F = ltratio_mex('single', kind, x, model, c, c0, g0, stat, sopt)    the single integral; kind 1 density, 2 cdf
 *   [F, err, status, nbad, nunres, nrefine] = ltratio_mex('double', kind, x, model, c, c0, g0, stat, dopt)    dopt(14):
 *                     verified (0: one evaluation, the default; 1: the agreement of two paths and the refinements)
 *   [F, err, status, nodes, level] = ltratio_mex('contour', kind, x, model, c, c0, g0, verified)    the contour integral, a
 *                     fixed start g0 ~= 0 in models 2, 3 only (dfl_contour): verified 0 (default) one evaluation (err NaN,
 *                     status 0); verified 1 error estimate, status (0, or 2 where the estimate exceeds 1e-10), outer nodes
 *                     and the level of the inner rules (1..4) per point
 *   [F, nstatus] = ltratio_mex('dist', kind, x, model, c, c0, g0, stat, method, verified)    method 0 automatic, 1 single,
 *                     2 double, 3 contour (a fixed start in models 2, 3; also the automatic choice there); verified 0
 *                     (default) one evaluation by one method, 1 the estimates, fallbacks and tail recomputation (dfl_dist)
 *   [q, F, iter, conv] = ltratio_mex('icdf', p, model, c, c0, g0, stat, method, tol, verified)    as in 'dist'
 *   msg = ltratio_mex('brittle', kind, model, c, c0, g0, stat, method)    the warning of a brittle regime of the
 *                     one-evaluation default ('' if none; dfl_brittle)
 *   [m, err] = ltratio_mex('mom', p, model, c, c0, g0, stat)             the moments 1..p
 *   ts = ltratio_mex('tstar', model, c, stat)                            the end of the domain of the mgf of V
 *   [F, err, status] = ltratio_mex('adfmax', kind, x, model, c, stat, dopt)    the ADFmax statistic (model 2, 3)
 *   [P, err, status] = ltratio_mex('region', type, par, v0, c, dopt)     P{U <= g(V)}, case DK, zero start: type 1 line
 *                     (par = [a b], g = b + a v), 2 parabola ([x c0]), 3 the weighted-average-power test ([k omega])
 *   P = ltratio_mex('lincdf', b, a, model, c, cbar, gopt)               P(U - aV <= b), zero start, model 1 (case DK) or 4
 *                     (case DG, detrended at cbar, default -13.5); a = NaN: P(V <= b) (case DK)
 *   P = ltratio_mex('envelope', c, model, alpha, gopt)                   the power envelope of Elliott, Rothenberg and Stock
 *                     (1996), model 1 (case DK, also DU) or 4 (case DG), level alpha (default 0.05)
 *   n = ltratio_mex('threads'[, n])                                      the number of threads (sets it with n)
 *
 * model 1, 2, 3 (no deterministic terms, intercept, intercept and trend); sopt = [tau hu N nodefac pfac edge_hw
 * uniform_hw gen], dopt = [M Nv Nl maxlev plen otol stepw errtol blo bhi bmin mumin tdec] and gopt = [N umax] (NaN or a
 * shorter vector: the defaults).  Compile with -DDFLIMIT_STUB where MATLAB's headers are not available (mexstub.h). */
#include "../dflimit.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef DFLIMIT_STUB
#include "mexstub.h"
DFLIMIT_DEFINE_APIVERSION
#else
#include "mex.h"
#endif

static double scalar(const mxArray *a, const char *name)
{
    if (!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1) mexErrMsgIdAndTxt("ltratio:args", "%s must be a real scalar", name);
    return mxGetScalar(a);
}
static int intarg(const mxArray *a, const char *name) { double v = scalar(a, name); if (v != floor(v)) mexErrMsgIdAndTxt("ltratio:args", "%s must be an integer", name); return (int) v; }
static const mxArray *vecarg(const mxArray *a, const char *name, int *n)
{
    if (!mxIsDouble(a) || mxIsComplex(a)) mexErrMsgIdAndTxt("ltratio:args", "%s must be a real array", name);
    *n = (int) mxGetNumberOfElements(a); return a;
}
static mxArray *like(const mxArray *a) { return mxCreateDoubleMatrix(mxGetM(a), mxGetN(a), mxREAL); }
/* an optional vector of option values, NaN meaning the default */
static void optvec(const mxArray *a, double *v, int n) { for (int k = 0; k < n; k++) v[k] = NAN; if (a == NULL) return; int m = (int) mxGetNumberOfElements(a); if (!mxIsDouble(a)) mexErrMsgIdAndTxt("ltratio:args", "the option vector must be double"); const double *p = mxGetPr(a); for (int k = 0; k < n && k < m; k++) v[k] = p[k]; }
#define SET(field, v) do { if (!isnan(v)) o.field = (v); } while (0)
static dfl_single_opts sopts(const mxArray *a)
{
    dfl_single_opts o; dfl_single_defaults(&o); double v[8]; optvec(a, v, 8);
    SET(tau, v[0]); SET(hu, v[1]); if (!isnan(v[2])) o.N = (int) v[2]; if (!isnan(v[3])) o.nodefac = (int) v[3]; if (!isnan(v[4])) o.pfac = (int) v[4];
    SET(edge_hw, v[5]); SET(uniform_hw, v[6]); if (!isnan(v[7])) o.gen = v[7] != 0;
    return o;
}
static dfl_double_opts dopts(const mxArray *a)
{
    dfl_double_opts o; dfl_double_defaults(&o); double v[14]; optvec(a, v, 14);
    if (!isnan(v[13])) o.verified = v[13] != 0;
    if (!isnan(v[0])) o.M = (int) v[0]; if (!isnan(v[1])) o.Nv = (int) v[1]; if (!isnan(v[2])) o.Nl = (int) v[2]; if (!isnan(v[3])) o.maxlev = (int) v[3];
    SET(plen, v[4]); SET(otol, v[5]); SET(stepw, v[6]); SET(errtol, v[7]); SET(blo, v[8]); SET(bhi, v[9]); SET(bmin, v[10]); SET(mumin, v[11]); SET(tdec, v[12]);
    return o;
}
static dfl_gp_opts gopts(const mxArray *a)
{
    dfl_gp_opts o; dfl_gp_defaults(&o); double v[2]; optvec(a, v, 2);
    if (!isnan(v[0])) o.N = (int) v[0]; SET(umax, v[1]);
    return o;
}
#undef SET
static void need(int nrhs, int n, const char *usage) { if (nrhs < n) mexErrMsgIdAndTxt("ltratio:args", "usage: %s", usage); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char op[16], err[512] = "";
    if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("ltratio:args", "the first argument is the operation: single, double, contour, dist, icdf, mom, tstar, adfmax, region, lincdf, envelope or threads");
    mxGetString(prhs[0], op, sizeof op);
    if (strcmp(op, "threads") == 0) {
        if (nrhs > 1) dfl_set_nworkers(intarg(prhs[1], "n"));
        plhs[0] = mxCreateDoubleScalar(dfl_nworkers()); return;
    }
    if (strcmp(op, "tstar") == 0) {
        need(nrhs, 3, "ltratio_mex('tstar', model, c, stat)");
        plhs[0] = mxCreateDoubleScalar(dfl_tstar(intarg(prhs[1], "model"), scalar(prhs[2], "c"), nrhs > 3 && scalar(prhs[3], "stat") != 0)); return;
    }
    if (strcmp(op, "mom") == 0) {
        need(nrhs, 3, "ltratio_mex('mom', p, model, c, c0, g0, stat)");
        int p = intarg(prhs[1], "p"); if (p < 1) mexErrMsgIdAndTxt("ltratio:args", "p must be a positive integer");
        plhs[0] = mxCreateDoubleMatrix(1, p, mxREAL); mxArray *e = mxCreateDoubleMatrix(1, p, mxREAL);
        int rc = dfl_moments(intarg(prhs[2], "model"), scalar(prhs[3], "c"), nrhs > 4 ? scalar(prhs[4], "c0") : 0.0, nrhs > 5 ? scalar(prhs[5], "g0") : 0.0,
                             nrhs > 6 && scalar(prhs[6], "stat") != 0, p, mxGetPr(plhs[0]), mxGetPr(e), err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        if (nlhs > 1) plhs[1] = e; return;
    }
    if (strcmp(op, "single") == 0) {
        need(nrhs, 5, "ltratio_mex('single', kind, x, model, c, c0, g0, stat, sopt)");
        int n; const mxArray *xa = vecarg(prhs[2], "x", &n); plhs[0] = like(xa);
        dfl_single_opts o = sopts(nrhs > 8 ? prhs[8] : NULL);
        int rc = dfl_single(intarg(prhs[1], "kind"), n, mxGetPr(xa), intarg(prhs[3], "model"), scalar(prhs[4], "c"), nrhs > 5 ? scalar(prhs[5], "c0") : 0.0,
                            nrhs > 6 ? scalar(prhs[6], "g0") : 0.0, nrhs > 7 && scalar(prhs[7], "stat") != 0, &o, mxGetPr(plhs[0]), err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        return;
    }
    if (strcmp(op, "double") == 0) {
        need(nrhs, 5, "ltratio_mex('double', kind, x, model, c, c0, g0, stat, dopt)");
        int n; const mxArray *xa = vecarg(prhs[2], "x", &n); plhs[0] = like(xa);
        dfl_double_opts o = dopts(nrhs > 8 ? prhs[8] : NULL);
        dfl_double_info *info = malloc((n + 1) * sizeof(dfl_double_info));
        int kind = intarg(prhs[1], "kind"); if (kind != 1 && kind != 2) mexErrMsgIdAndTxt("ltratio:args", "kind must be 1 (density) or 2 (distribution function)");
        int model = intarg(prhs[3], "model"); double c = scalar(prhs[4], "c"), c0 = nrhs > 5 ? scalar(prhs[5], "c0") : 0.0, g0 = nrhs > 6 ? scalar(prhs[6], "g0") : 0.0;
        bool stat = nrhs > 7 && scalar(prhs[7], "stat") != 0;
        if (model < 1 || model > 4) mexErrMsgIdAndTxt("ltratio:args", "model must be 1, 2, 3 or 4");
        if (stat && !(c <= DFL_CSTAT)) mexErrMsgIdAndTxt("ltratio:args", "%s", DFL_STAT_MSG);
        /* finite x only: the m-file handles the rest */
        const double *x = mxGetPr(xa); for (int k = 0; k < n; k++) if (!isfinite(x[k])) mexErrMsgIdAndTxt("ltratio:args", "x must be finite");
        int rc = n > 0 ? dfl_double(kind, n, x, model, c, c0, g0, stat, &o, mxGetPr(plhs[0]), info, err, sizeof err) : 0;
        if (rc) { free(info); mexErrMsgIdAndTxt("ltratio:engine", "%s", err); }
        for (int j = 1; j < 6 && j < nlhs; j++) {
            plhs[j] = like(xa); double *p = mxGetPr(plhs[j]);
            for (int k = 0; k < n; k++) p[k] = j == 1 ? info[k].err : j == 2 ? info[k].status : j == 3 ? info[k].nbad : j == 4 ? info[k].nunres : info[k].nrefine;
        }
        free(info); return;
    }
    if (strcmp(op, "contour") == 0) {
        need(nrhs, 7, "ltratio_mex('contour', kind, x, model, c, c0, g0, verified)");
        int n; const mxArray *xa = vecarg(prhs[2], "x", &n); plhs[0] = like(xa);
        double *e = malloc((n + 1) * sizeof(double)); int *st = malloc((n + 1) * sizeof(int)), *nd = malloc((n + 1) * sizeof(int)), *lv = malloc((n + 1) * sizeof(int));
        int rc = dfl_contour_lev(intarg(prhs[1], "kind"), n, mxGetPr(xa), intarg(prhs[3], "model"), scalar(prhs[4], "c"), scalar(prhs[5], "c0"), scalar(prhs[6], "g0"),
                                 nrhs > 7 && scalar(prhs[7], "verified") != 0, mxGetPr(plhs[0]), e, st, nd, lv, err, sizeof err);
        if (rc) { free(e); free(st); free(nd); free(lv); mexErrMsgIdAndTxt("ltratio:engine", "%s", err); }
        for (int j = 1; j < 5 && j < nlhs; j++) {
            plhs[j] = like(xa); double *p = mxGetPr(plhs[j]);
            for (int k = 0; k < n; k++) p[k] = j == 1 ? e[k] : j == 2 ? st[k] : j == 3 ? nd[k] : lv[k];
        }
        free(e); free(st); free(nd); free(lv); return;
    }
    if (strcmp(op, "dist") == 0) {
        need(nrhs, 5, "ltratio_mex('dist', kind, x, model, c, c0, g0, stat, method, verified)");
        int n; const mxArray *xa = vecarg(prhs[2], "x", &n); plhs[0] = like(xa); int ns = 0;
        dfl_double_opts dv; dfl_double_defaults(&dv); dv.verified = nrhs > 9 && scalar(prhs[9], "verified") != 0;
        int rc = dfl_dist(intarg(prhs[1], "kind"), n, mxGetPr(xa), intarg(prhs[3], "model"), scalar(prhs[4], "c"), nrhs > 5 ? scalar(prhs[5], "c0") : 0.0,
                          nrhs > 6 ? scalar(prhs[6], "g0") : 0.0, nrhs > 7 && scalar(prhs[7], "stat") != 0, nrhs > 8 ? intarg(prhs[8], "method") : 0,
                          NULL, &dv, mxGetPr(plhs[0]), &ns, err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(ns);
        return;
    }
    if (strcmp(op, "brittle") == 0) {
        need(nrhs, 8, "ltratio_mex('brittle', kind, model, c, c0, g0, stat, method)");
        const char *bm = dfl_brittle(intarg(prhs[1], "kind"), intarg(prhs[2], "model"), scalar(prhs[3], "c"), scalar(prhs[4], "c0"), scalar(prhs[5], "g0"),
                                     scalar(prhs[6], "stat") != 0, intarg(prhs[7], "method"));
        plhs[0] = mxCreateString(bm ? bm : "");
        return;
    }
    if (strcmp(op, "icdf") == 0) {
        need(nrhs, 4, "ltratio_mex('icdf', p, model, c, c0, g0, stat, method, tol, verified)");
        int n; const mxArray *pa = vecarg(prhs[1], "p", &n); plhs[0] = like(pa);
        mxArray *Fq = like(pa), *it = like(pa), *cv = like(pa); int *iter = malloc((n + 1) * sizeof(int)), *conv = malloc((n + 1) * sizeof(int));
        dfl_double_opts dvi; dfl_double_defaults(&dvi); dvi.verified = nrhs > 9 && scalar(prhs[9], "verified") != 0;
        int rc = dfl_icdf(n, mxGetPr(pa), intarg(prhs[2], "model"), scalar(prhs[3], "c"), nrhs > 4 ? scalar(prhs[4], "c0") : 0.0, nrhs > 5 ? scalar(prhs[5], "g0") : 0.0,
                          nrhs > 6 && scalar(prhs[6], "stat") != 0, nrhs > 7 ? intarg(prhs[7], "method") : 0, nrhs > 8 ? scalar(prhs[8], "tol") : 1e-10,
                          NULL, &dvi, mxGetPr(plhs[0]), mxGetPr(Fq), iter, conv, err, sizeof err);
        for (int k = 0; k < n; k++) { mxGetPr(it)[k] = iter[k]; mxGetPr(cv)[k] = conv[k]; }
        free(iter); free(conv);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        if (nlhs > 1) plhs[1] = Fq; if (nlhs > 2) plhs[2] = it; if (nlhs > 3) plhs[3] = cv;
        return;
    }
    if (strcmp(op, "adfmax") == 0) {
        need(nrhs, 6, "ltratio_mex('adfmax', kind, x, model, c, stat, dopt)");
        int n; const mxArray *xa = vecarg(prhs[2], "x", &n); plhs[0] = like(xa);
        dfl_double_opts o = dopts(nrhs > 6 ? prhs[6] : NULL);
        dfl_double_info *info = malloc((n + 1) * sizeof(dfl_double_info));
        int rc = dfl_adfmax(intarg(prhs[1], "kind"), n, mxGetPr(xa), intarg(prhs[3], "model"), scalar(prhs[4], "c"), scalar(prhs[5], "stat") != 0, &o, mxGetPr(plhs[0]), info, err, sizeof err);
        if (rc) { free(info); mexErrMsgIdAndTxt("ltratio:engine", "%s", err); }
        for (int j = 1; j < 3 && j < nlhs; j++) { plhs[j] = like(xa); double *p = mxGetPr(plhs[j]); for (int k = 0; k < n; k++) p[k] = j == 1 ? info[k].err : info[k].status; }
        free(info); return;
    }
    if (strcmp(op, "region") == 0) {
        need(nrhs, 5, "ltratio_mex('region', type, par, v0, c, dopt)");
        int np; const mxArray *pa = vecarg(prhs[2], "par", &np); if (np < 2) mexErrMsgIdAndTxt("ltratio:args", "par must have two elements");
        dfl_double_opts o = dopts(nrhs > 5 ? prhs[5] : NULL); dfl_double_info info; double P;
        int rc = dfl_region(intarg(prhs[1], "type"), mxGetPr(pa), scalar(prhs[3], "v0"), 0, NULL, scalar(prhs[4], "c"), &o, &P, &info, err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        plhs[0] = mxCreateDoubleScalar(P); if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(info.err); if (nlhs > 2) plhs[2] = mxCreateDoubleScalar(info.status);
        return;
    }
    if (strcmp(op, "lincdf") == 0) {
        need(nrhs, 5, "ltratio_mex('lincdf', b, a, model, c, cbar, gopt)");
        int n; const mxArray *ba = vecarg(prhs[1], "b", &n); plhs[0] = like(ba);
        dfl_gp_opts o = gopts(nrhs > 6 ? prhs[6] : NULL);
        int rc = dfl_lincdf(n, mxGetPr(ba), scalar(prhs[2], "a"), intarg(prhs[3], "model"), scalar(prhs[4], "c"), nrhs > 5 ? scalar(prhs[5], "cbar") : DFL_CBAR_GLS,
                            &o, mxGetPr(plhs[0]), err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        return;
    }
    if (strcmp(op, "envelope") == 0) {
        need(nrhs, 3, "ltratio_mex('envelope', c, model, alpha, gopt)");
        int n; const mxArray *ca = vecarg(prhs[1], "c", &n); plhs[0] = like(ca);
        dfl_gp_opts o = gopts(nrhs > 4 ? prhs[4] : NULL);
        int rc = dfl_power_envelope(n, mxGetPr(ca), intarg(prhs[2], "model"), nrhs > 3 ? scalar(prhs[3], "alpha") : 0.05, &o, mxGetPr(plhs[0]), err, sizeof err);
        if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
        return;
    }
    mexErrMsgIdAndTxt("ltratio:args", "unknown operation %s", op);
}
