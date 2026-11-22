/* ltratio_mex.c: MATLAB gateway to the C engine (c/), one MEX file for every operation of the MATLAB programs; the
 * m-files dfx_core.m, tr_limdist.m, tr_contour.m and tr_limmom.m next to it present the engine under the names and signatures of
 * the MATLAB engines, so that ltratio_cdf, ltratio_pdf, ltratio_icdf, ltratio_mom and ltratio_test run unchanged on
 * top of it, and ltratio.m is a self-contained wrapper.  All arguments are double arrays (logical flags as 0/1):
 *
 *   F = ltratio_mex('single', kind, x, model, c, c0, g0, stat, sopt)    the single integral; kind 1 density, 2 cdf; at a
 *                     point beyond the bounds of the engine the error dfx:bound:rule or dfx:bound:path (as the MATLAB
 *                     engine), which LTRATIO_DIST routes to the conditional integral
 *   [F, err, status, nbad, nunres, nrefine] = ltratio_mex('double', kind, x, model, c, c0, g0, stat, dopt)    dopt(14):
 *                     verified (0: one evaluation, the default; 1: the agreement of two paths and the refinements)
 *   [F, err, status, nodes, level] = ltratio_mex('contour', kind, x, model, c, c0, g0, verified)    the contour integral, a
 *                     fixed start g0 ~= 0 in models 2, 3 only (dfl_contour): verified 0 (default) one evaluation (err NaN,
 *                     status 0); verified 1 error estimate, status (0, or 2 where the estimate exceeds 1e-10), outer nodes
 *                     and the level of the inner rules (1..4) per point
 *   [F, nstatus] = ltratio_mex('dist', kind, x, model, c, c0, g0, stat, method, verified)    method 0 automatic, 1 single,
 *                     2 double, 3 contour (a fixed start in models 2, 3; also the automatic choice there); verified 0
 *                     (default) one evaluation by one method, 1 the estimates, fallbacks and tail recomputation (dfl_dist);
 *                     the warning ltratio:bound where the line integral refused points whose values are flagged
 *                     (dfl_warning), the error ltratio:bound where it refused a point with method 1
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
 *   [br, C, loc, first] = ltratio_mex('design', kind, ...)               a design of the general engine (dfl_design): kind
 *                     'trend' p, 'perron' model lam (1 A, 2 B, 3 C), 'ltrend' p, 'matrix' X TB, 'first' br C loc first,
 *                     'fun' f K breaks (a function handle); br the break points, C the coefficient matrices side by side
 *   v = ltratio_mex('det', what, br, C, loc, first, ..., opts)          the general engine (dfl_det_*; see det_op below)
 *   v = ltratio_mex('route', what, engine, spec, ...)                    the router over the two engines (dfl_engine_*)
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

/* ---- the general engine (dfl_det_*): a design crosses the gateway as numeric arrays, br (1 x (np + 1)), C ((deg + 1) x
 * (K np): the coefficient matrices of the pieces side by side), loc and first ---- */
static void design_in(const mxArray *const *a, dfl_design *D)
{
    int nb; vecarg(a[0], "br", &nb);
    if (nb < 2) mexErrMsgIdAndTxt("ltratio:args", "a design needs at least two break points (0 and 1)");
    if (!mxIsDouble(a[1]) || mxIsComplex(a[1])) mexErrMsgIdAndTxt("ltratio:args", "the coefficients of a design must be a real matrix");
    D->np = nb - 1; D->deg = (int) mxGetM(a[1]) - 1; int cols = (int) mxGetN(a[1]);
    if (mxGetNumberOfElements(a[1]) == 0) { D->deg = 0; cols = 0; }
    if (cols % D->np) mexErrMsgIdAndTxt("ltratio:args", "the coefficient matrix of a design has K np columns");
    D->K = cols / D->np; D->br = mxGetPr(a[0]); D->C = mxGetPr(a[1]);
    D->loc = scalar(a[2], "loc") != 0; D->first = scalar(a[3], "first") != 0;
}
static void design_out(int nlhs, mxArray *plhs[], const dfl_design *D)
{
    plhs[0] = mxCreateDoubleMatrix(1, D->np + 1, mxREAL); memcpy(mxGetPr(plhs[0]), D->br, (D->np + 1) * sizeof(double));
    if (nlhs > 1) { plhs[1] = mxCreateDoubleMatrix(D->deg + 1, (size_t) D->K * D->np, mxREAL); if (D->K > 0) memcpy(mxGetPr(plhs[1]), D->C, (size_t) (D->deg + 1) * D->K * D->np * sizeof(double)); }
    if (nlhs > 2) plhs[2] = mxCreateDoubleScalar(D->loc);
    if (nlhs > 3) plhs[3] = mxCreateDoubleScalar(D->first);
}
typedef struct { const mxArray *fh; int K; int bad; } mfun_ctx;
static void mfun(double r, double *out, void *v)
{
    mfun_ctx *m = v; mxArray *in[2], *res[1];
    in[0] = (mxArray *) m->fh; in[1] = mxCreateDoubleScalar(r);
    mexCallMATLAB(1, res, 2, in, "feval");
    if (!mxIsDouble(res[0]) || (int) mxGetNumberOfElements(res[0]) != m->K) { m->bad = 1; for (int i = 0; i < m->K; i++) out[i] = 0; return; }
    memcpy(out, mxGetPr(res[0]), m->K * sizeof(double));
}
/* [br, C, loc, first] = ltratio_mex('design', kind, ...): 'trend' p | 'perron' model lam (model 1 A, 2 B, 3 C) | 'ltrend' p |
 * 'matrix' X TB (TB empty: the breaks detected) | 'first' br C loc first | 'fun' f K breaks (f a function handle of r returning
 * the K regressors) */
static void det_design_op(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char kind[16], err[512]; dfl_design D; int rc = 1;
    need(nrhs, 3, "ltratio_mex('design', kind, ...)");
    if (!mxIsChar(prhs[1])) mexErrMsgIdAndTxt("ltratio:args", "the kind of the design must be a string");
    mxGetString(prhs[1], kind, sizeof kind);
    memset(&D, 0, sizeof D);
    if (!strcmp(kind, "trend")) rc = dfl_design_trend(intarg(prhs[2], "p"), &D, err, sizeof err);
    else if (!strcmp(kind, "ltrend")) rc = dfl_design_ltrend(intarg(prhs[2], "p"), &D, err, sizeof err);
    else if (!strcmp(kind, "perron")) {
        need(nrhs, 4, "ltratio_mex('design', 'perron', model, lam)");
        int m = intarg(prhs[2], "model"), q; const mxArray *la = vecarg(prhs[3], "lam", &q);
        rc = dfl_design_perron(m == 1 ? 'A' : m == 2 ? 'B' : m == 3 ? 'C' : '?', q, mxGetPr(la), &D, err, sizeof err);
    } else if (!strcmp(kind, "matrix")) {
        if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) mexErrMsgIdAndTxt("ltratio:args", "X must be a real matrix");
        int n = (int) mxGetM(prhs[2]), K = (int) mxGetN(prhs[2]), nb = -1; int *TB = NULL;
        if (nrhs > 3 && mxGetNumberOfElements(prhs[3]) > 0) {
            vecarg(prhs[3], "TB", &nb); TB = malloc(nb * sizeof(int));
            for (int i = 0; i < nb; i++) TB[i] = (int) mxGetPr(prhs[3])[i];
        }
        rc = dfl_design_from_matrix(n, K, mxGetPr(prhs[2]), nb, TB, NAN, 0, 0, &D, err, sizeof err);
        free(TB);
        if (!rc && dfl_det_warning()) mexWarnMsgIdAndTxt("ltratio:design", "%s", dfl_det_warning());
    } else if (!strcmp(kind, "first")) {
        need(nrhs, 6, "ltratio_mex('design', 'first', br, C, loc, first)");
        dfl_design B; design_in(prhs + 2, &B);
        rc = dfl_design_detrend_first(&B, &D, err, sizeof err);
    } else if (!strcmp(kind, "fun")) {
        need(nrhs, 4, "ltratio_mex('design', 'fun', f, K, breaks)");
        mfun_ctx m = {prhs[2], intarg(prhs[3], "K"), 0}; int nb = 0; const double *b = NULL;
        if (nrhs > 4) { vecarg(prhs[4], "breaks", &nb); b = mxGetPr(prhs[4]); }
        rc = dfl_design_fdesign(mfun, &m, m.K, nb, b, NAN, 0, 0, &D, err, sizeof err);
        if (!rc && m.bad) { dfl_design_free(&D); mexErrMsgIdAndTxt("ltratio:args", "the function of the design must return a real vector of K values"); }
    } else mexErrMsgIdAndTxt("ltratio:args", "the kind of design must be trend, perron, ltrend, matrix, first or fun");
    if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
    design_out(nlhs, plhs, &D);
    dfl_design_free(&D);
}
/* the options of the general engine: [c0 start g0 gls verified] (NaN or a shorter vector: the defaults) */
static dfl_det_opts detopts(const mxArray *a)
{
    dfl_det_opts o; dfl_det_defaults(&o); double v[5]; optvec(a, v, 5);
    if (!isnan(v[0])) o.c0 = v[0]; if (!isnan(v[1])) o.start = (int) v[1]; if (!isnan(v[2])) o.gamma0 = v[2]; o.gls = v[3]; if (!isnan(v[4])) o.verified = v[4] != 0;
    return o;
}
/* v = ltratio_mex('det', what, br, C, loc, first, a1, a2, ..., opts): the general engine (dfl_det_*); what and its arguments
 *   'cdf' 'pdf' x c | 'icdf' p c x0 tol | 'coefcdf' 'coefpdf' y c | 'coeficdf' p c y0 tol | 'mom' pmax c |
 *   'lincdf' 'linpdf' b a c (a = NaN: the law of V) | 'linicdf' p a c b0 tol | 'adfcdf' 'adfpdf' x c | 'adficdf' p c x0 tol |
 *   'power' cs level statistic (0 t-ratio, 1 coefficient) | 'envelope' cs level | 'ptcdf' 'ptpdf' q cbar c | 'pticdf' p cbar |
 *   'ptenvelope' cs level | 'lmcdf' x c | 'lmicdf' p c | 'mzacdf' y c gls | 'mzaicdf' p c gls | 'msbcdf' x c gls | 'msbicdf' p c gls
 * opts = [c0 start g0 gls verified] after the arguments (start 0 zero, 1 fixed, 2 stationary); the warnings of the engine
 * (dfl_det_warning) as the warning ltratio:det */
static void det_op(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char what[16], err[512]; int rc = 1;
    (void) nlhs;
    need(nrhs, 8, "ltratio_mex('det', what, br, C, loc, first, ...)");
    if (!mxIsChar(prhs[1])) mexErrMsgIdAndTxt("ltratio:args", "what must be a string");
    mxGetString(prhs[1], what, sizeof what);
    dfl_design D; design_in(prhs + 2, &D);
    const mxArray *const *a = prhs + 6; int na = nrhs - 6;
#define OPTS(k) detopts(na > (k) ? a[k] : NULL)
#define SC(k, nm, dflt) (na > (k) && mxGetNumberOfElements(a[k]) > 0 ? scalar(a[k], nm) : (dflt))
    int n; const mxArray *xa = vecarg(a[0], "x", &n);
    if (!strcmp(what, "cdf") || !strcmp(what, "pdf")) { dfl_det_opts o = OPTS(2); plhs[0] = like(xa); rc = dfl_det_dist(what[0] == 'c' ? 2 : 1, &D, n, mxGetPr(xa), SC(1, "c", 0.0), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "icdf")) { dfl_det_opts o = OPTS(4); plhs[0] = like(xa); rc = dfl_det_quantile(&D, n, mxGetPr(xa), SC(1, "c", 0.0), SC(2, "x0", -2.0), SC(3, "tol", 1e-9), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "coefcdf") || !strcmp(what, "coefpdf")) { dfl_det_opts o = OPTS(2); plhs[0] = like(xa); rc = dfl_det_coef(what[4] == 'c' ? 2 : 1, &D, n, mxGetPr(xa), SC(1, "c", 0.0), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "coeficdf")) { dfl_det_opts o = OPTS(4); plhs[0] = like(xa); rc = dfl_det_coef_quantile(&D, n, mxGetPr(xa), SC(1, "c", 0.0), SC(2, "y0", -10.0), SC(3, "tol", 1e-9), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "mom")) { dfl_det_opts o = OPTS(2); int p = intarg(a[0], "pmax"); plhs[0] = mxCreateDoubleMatrix(1, p, mxREAL); rc = dfl_det_moments(&D, SC(1, "c", 0.0), p, &o, mxGetPr(plhs[0]), NULL, err, sizeof err); }
    else if (!strcmp(what, "lincdf") || !strcmp(what, "linpdf")) { dfl_det_opts o = OPTS(3); plhs[0] = like(xa); rc = dfl_det_lin(what[3] == 'c' ? 2 : 1, &D, n, mxGetPr(xa), SC(1, "a", 0.0), SC(2, "c", 0.0), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "linicdf")) {
        dfl_det_opts o = OPTS(5); plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_lin_quantile(&D, mxGetPr(xa)[k], SC(1, "a", 0.0), SC(2, "c", 0.0), SC(3, "b0", 0.0), SC(4, "tol", 1e-10), &o, mxGetPr(plhs[0]) + k, err, sizeof err);
    } else if (!strcmp(what, "adfcdf") || !strcmp(what, "adfpdf")) { dfl_det_opts o = OPTS(2); plhs[0] = like(xa); rc = dfl_det_adfmax(what[3] == 'c' ? 2 : 1, &D, n, mxGetPr(xa), SC(1, "c", 0.0), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "adficdf")) {
        dfl_det_opts o = OPTS(4); plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_adfmax_quantile(&D, mxGetPr(xa)[k], SC(1, "c", 0.0), SC(2, "x0", -2.0), SC(3, "tol", 1e-10), &o, mxGetPr(plhs[0]) + k, err, sizeof err);
    } else if (!strcmp(what, "power")) { dfl_det_opts o = OPTS(3); plhs[0] = like(xa); rc = dfl_det_local_power(&D, n, mxGetPr(xa), SC(1, "level", 0.05), (int) SC(2, "statistic", 0.0), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "envelope")) { plhs[0] = like(xa); rc = dfl_det_power_envelope(&D, n, mxGetPr(xa), SC(1, "level", 0.05), mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "ptcdf") || !strcmp(what, "ptpdf")) { plhs[0] = like(xa); rc = dfl_det_pt(what[2] == 'c' ? 2 : 1, &D, n, mxGetPr(xa), SC(1, "cbar", -7.0), SC(2, "c", 0.0), 0, 0, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "pticdf")) {
        plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_pt_quantile(&D, mxGetPr(xa)[k], SC(1, "cbar", -7.0), 0, 0, mxGetPr(plhs[0]) + k, err, sizeof err);
    } else if (!strcmp(what, "ptenvelope")) { plhs[0] = like(xa); rc = dfl_det_pt_envelope(&D, n, mxGetPr(xa), SC(1, "level", 0.05), 0, 0, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "lmcdf")) { plhs[0] = like(xa); rc = dfl_det_lm_cdf(&D, n, mxGetPr(xa), SC(1, "c", 0.0), mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "lmicdf")) {
        plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_lm_quantile(&D, mxGetPr(xa)[k], SC(1, "c", 0.0), mxGetPr(plhs[0]) + k, err, sizeof err);
    } else if (!strcmp(what, "mzacdf")) { dfl_det_opts o = OPTS(3); plhs[0] = like(xa); rc = dfl_det_mza_cdf(&D, n, mxGetPr(xa), SC(1, "c", 0.0), SC(2, "gls", -13.5), &o, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "mzaicdf")) {
        dfl_det_opts o = OPTS(3); plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_mza_quantile(&D, mxGetPr(xa)[k], SC(1, "c", 0.0), SC(2, "gls", -13.5), NAN, 1e-9, &o, mxGetPr(plhs[0]) + k, err, sizeof err);
    } else if (!strcmp(what, "msbcdf")) { plhs[0] = like(xa); rc = dfl_det_msb_cdf(&D, n, mxGetPr(xa), SC(1, "c", 0.0), SC(2, "gls", -13.5), mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "msbicdf")) {
        plhs[0] = like(xa); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_msb_quantile(&D, mxGetPr(xa)[k], SC(1, "c", 0.0), SC(2, "gls", -13.5), mxGetPr(plhs[0]) + k, err, sizeof err);
    } else mexErrMsgIdAndTxt("ltratio:args", "unknown function %s of the general engine", what);
#undef OPTS
#undef SC
    if (rc) mexErrMsgIdAndTxt("ltratio:engine", "%s", err);
    if (dfl_det_warning()) mexWarnMsgIdAndTxt("ltratio:det", "%s", dfl_det_warning());
}
/* v = ltratio_mex('route', what, engine, spec, a1, ...): one interface over the two engines (dfl_engine_*); engine 0 auto,
 * 1 classical, 2 general; spec = [model start g0] of a classical case, or {br, C, loc, first, start, g0, gls} as the
 * numeric arguments br, C, loc, first, start, g0, gls (spec given as the scalar 0 then these seven); what:
 *   'cdf' 'pdf' x c c0 | 'icdf' p c c0 tol | 'mom' pmax c c0 | 'coefcdf' y c c0 | 'coeficdf' p c c0 tol |
 *   'adfcdf' 'adfpdf' x c | 'lincdf' b c a cbar (the model of spec: 1 case DK, 4 case DG detrended at cbar, NaN for
 *   -13.5) | 'envelope' cs alpha (the model of spec, 1 or 4) | 'region' type c p1 p2 (P{U <= g(V)}: type 1 the line
 *   p2 + p1 v, 2 the parabola p1 sqrt(v) + p2 v, 3 the weighted-average-power boundary with k = p1, omega = p2)
 * The warnings: dfl_warning as ltratio:bound, dfl_det_warning as ltratio:det, and for a classical case with engine auto
 * the warning ltratio:brittle of a regime in which the classical engine's one evaluation is brittle, where the router
 * takes the classical engine at a point (dfl_brittle, once per message; as the classical entry points give it) */
static int brittle_add(const char **bm, int nb, char store[][320], const char *m)
{
    if (!m || nb >= 8) return nb;
    for (int j = 0; j < nb; j++) if (!strcmp(bm[j], m)) return nb;
    snprintf(store[nb], 320, "%s", m); bm[nb] = store[nb];
    return nb + 1;
}
static void route_op(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char what[16], err[512]; int rc = 1;
    (void) nlhs;
    need(nrhs, 5, "ltratio_mex('route', what, engine, spec, ...)");
    mxGetString(prhs[1], what, sizeof what);
    int engine = intarg(prhs[2], "engine");
    dfl_spec S; dfl_design D; int k0;
    int ns; vecarg(prhs[3], "spec", &ns);
    if (ns == 3) { const double *sp = mxGetPr(prhs[3]); S.D = NULL; S.model = (int) sp[0]; S.start = (int) sp[1]; S.g0 = sp[2]; S.gls = NAN; k0 = 4; }
    else {
        need(nrhs, 11, "ltratio_mex('route', what, engine, 0, br, C, loc, first, start, g0, gls, ...)");
        design_in(prhs + 4, &D); S.D = &D; S.model = 0; S.start = intarg(prhs[8], "start"); S.g0 = scalar(prhs[9], "g0"); S.gls = scalar(prhs[10], "gls"); k0 = 11;
    }
    const mxArray *const *a = prhs + k0; int na = nrhs - k0;
    if (na < 2) mexErrMsgIdAndTxt("ltratio:args", "the route needs the points and c");
#define SC(k, nm, dflt) (na > (k) && mxGetNumberOfElements(a[k]) > 0 ? scalar(a[k], nm) : (dflt))
    int n; const mxArray *xa = vecarg(a[0], "x", &n);
    double c = SC(1, "c", 0.0), c0 = SC(2, "c0", 0.0);
    /* (a classical case with engine auto: the brittle regimes of the points that the router gives the classical engine) */
    int cls = S.D == NULL && engine == DFL_ENGINE_AUTO && S.model >= 1 && S.model <= 4, stat = S.start == 2;
    double g0 = stat ? 0.0 : S.g0, gls = S.model == 4 ? DFL_CBAR_GLS : NAN;
    const char *bm[8]; char bstore[8][320]; int nb = 0;
    if (!strcmp(what, "cdf") || !strcmp(what, "pdf")) {
        int kind = what[0] == 'c' ? 2 : 1; const double *x = mxGetPr(xa);
        plhs[0] = like(xa); rc = dfl_engine_dist(kind, n, x, &S, c, c0, engine, 0, 0, mxGetPr(plhs[0]), err, sizeof err);
        if (!rc && cls) {
            const double *v = mxGetPr(plhs[0]); int any = 0;
            for (int k = 0; k < n && !any; k++) {
                if (!isfinite(x[k])) continue;
                if (dfl_autoengine(DFL_F_DIST, S.model, stat, g0, gls, c, c0, x[k]) == DFL_ENGINE_CLASSICAL) any = 1;
                else if (v[k] < 1e-3 && dfl_automethod(S.model, c, c0, g0, stat, x[k]) == 2) any = 1;   /* (a tail sent back) */
            }
            if (any) nb = brittle_add(bm, nb, bstore, dfl_brittle(kind, S.model, c, c0, g0, stat, 0));
        }
    } else if (!strcmp(what, "icdf")) {
        plhs[0] = like(xa); rc = dfl_engine_icdf(n, mxGetPr(xa), &S, c, c0, engine, 0, SC(3, "tol", 1e-10), mxGetPr(plhs[0]), err, sizeof err);
        if (!rc && cls && dfl_autoengine(DFL_F_ICDF, S.model, stat, g0, gls, c, c0, 0.0) == DFL_ENGINE_CLASSICAL) nb = brittle_add(bm, nb, bstore, dfl_brittle(2, S.model, c, c0, g0, stat, 0));
    }
    else if (!strcmp(what, "mom")) { int p = intarg(a[0], "pmax"); plhs[0] = mxCreateDoubleMatrix(1, p, mxREAL); rc = dfl_engine_mom(p, &S, c, c0, engine, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "coefcdf") || !strcmp(what, "coeficdf")) {
        int q = what[4] == 'i';
        plhs[0] = like(xa);
        rc = q ? dfl_engine_coef_icdf(n, mxGetPr(xa), &S, c, c0, engine, SC(3, "tol", 1e-9), mxGetPr(plhs[0]), err, sizeof err)
               : dfl_engine_coef(n, mxGetPr(xa), &S, c, c0, engine, mxGetPr(plhs[0]), err, sizeof err);
        if (!rc && cls) for (int k = 0; k < n; k++) {
            /* (the t-ratio at 0 with the null value c0 + y: its engine and its brittle regime at each y, or at each quantile) */
            double y = q ? mxGetPr(plhs[0])[k] : mxGetPr(xa)[k], c0y = c0 + y;
            if (isfinite(y) && dfl_autoengine(DFL_F_DIST, S.model, stat, g0, gls, c, c0y, 0.0) == DFL_ENGINE_CLASSICAL) nb = brittle_add(bm, nb, bstore, dfl_brittle(2, S.model, c, c0y, g0, stat, 0));
        }
    }
    else if (!strcmp(what, "adfcdf") || !strcmp(what, "adfpdf")) { plhs[0] = like(xa); rc = dfl_engine_adfmax(what[3] == 'c' ? 2 : 1, n, mxGetPr(xa), &S, c, engine, mxGetPr(plhs[0]), err, sizeof err); }
    else if (!strcmp(what, "lincdf")) {
        if (S.D) mexErrMsgIdAndTxt("ltratio:args", "the linear form of the router takes a case (spec [model start g0], model 1 or 4)");
        plhs[0] = like(xa); rc = dfl_engine_lincdf(n, mxGetPr(xa), SC(2, "a", 0.0), S.model, c, SC(3, "cbar", NAN), engine, mxGetPr(plhs[0]), err, sizeof err);
    } else if (!strcmp(what, "envelope")) {
        if (S.D) mexErrMsgIdAndTxt("ltratio:args", "the envelope of the router takes a case (spec [model start g0], model 1 or 4)");
        plhs[0] = like(xa); rc = dfl_engine_envelope(n, mxGetPr(xa), S.model, SC(1, "alpha", 0.05), engine, mxGetPr(plhs[0]), err, sizeof err);
    } else if (!strcmp(what, "region")) {
        int type = intarg(a[0], "type"); double par[2] = {SC(2, "p1", NAN), SC(3, "p2", NAN)};
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); rc = dfl_engine_region(type, par, &S, c, engine, mxGetPr(plhs[0]), err, sizeof err);
    }
    else mexErrMsgIdAndTxt("ltratio:args", "unknown function %s of the router", what);
#undef SC
    if (rc) mexErrMsgIdAndTxt(rc > 10 ? "ltratio:bound" : "ltratio:engine", "%s", err);
    for (int j = 0; j < nb; j++) mexWarnMsgIdAndTxt("ltratio:brittle", "%s", bm[j]);
    if (dfl_warning()) mexWarnMsgIdAndTxt("ltratio:bound", "%s", dfl_warning());
    if (dfl_det_warning()) mexWarnMsgIdAndTxt("ltratio:det", "%s", dfl_det_warning());
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char op[16], err[512] = "";
    if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("ltratio:args", "the first argument is the operation: single, double, contour, dist, icdf, mom, tstar, adfmax, region, lincdf, envelope, threads, design, det or route");
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
        /* a point that the line integral refuses (its paths would exceed the bounds of the engine): the identifiers of the
         * MATLAB engine, dfx:bound:rule and dfx:bound:path, so that LTRATIO_DIST routes it as there */
        if (rc == 10 + DFL_BOUND_RULE) mexErrMsgIdAndTxt("dfx:bound:rule", "%s", err);
        if (rc == 10 + DFL_BOUND_PATH) mexErrMsgIdAndTxt("dfx:bound:path", "%s", err);
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
        if (rc) mexErrMsgIdAndTxt(rc > 10 ? "ltratio:bound" : "ltratio:engine", "%s", err);
        if (dfl_warning()) mexWarnMsgIdAndTxt("ltratio:bound", "%s", dfl_warning());   /* (refused points whose values are flagged) */
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
        if (dfl_warning()) mexWarnMsgIdAndTxt("ltratio:bound", "%s", dfl_warning());
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
    if (strcmp(op, "design") == 0) { det_design_op(nlhs, plhs, nrhs, prhs); return; }
    if (strcmp(op, "det") == 0) { det_op(nlhs, plhs, nrhs, prhs); return; }
    if (strcmp(op, "route") == 0) { route_op(nlhs, plhs, nrhs, prhs); return; }
    mexErrMsgIdAndTxt("ltratio:args", "unknown operation %s", op);
}
