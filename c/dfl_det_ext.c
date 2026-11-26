/* dfl_det_ext.c: the extensions of the limit of the general engine (a port of src/DFDet_vol.jl, src/DFDet_cov.jl and
 * src/DFDet_seas.jl of DFLimit.jl): the checks of a volatility path; a normal component, T = rho tau_c0 + sqrt(1 - rho^2) Z,
 * the limit of the covariate-augmented t-ratio and of the t-ratio of the slope in a predictive regression with a regressor
 * local to unity (the inversion on two lines with the normal component in the kernel, the normal average of the t-ratio's
 * distribution function, the local null values of the two applications); the first seasonal t-ratio of the regression of
 * Hylleberg, Engle, Granger and Yoo at a harmonic frequency (the sum of the residues on the outer contour).  The transform
 * and the paths of the volatility regimes are in dfl_det_mgf.c and dfl_det_path.c, their user-level functions in
 * dfl_det_api.c (dfl_det_dist and dfl_det_quantile with o->vol). */
#include "dfl_det_int.h"

int dfl_vol_check(const dfl_vol *v, char *err, size_t errlen)
{
    if (!v || v->nreg < 1) { snprintf(err, errlen, "volpath: a volatility path has at least one regime (nreg >= 1)"); return 1; }
    for (int i = 0; i < v->nreg - 1; i++) if (!(isfinite(v->breaks[i]) && v->breaks[i] > 0 && v->breaks[i] < 1)) { snprintf(err, errlen, "volpath: the break dates must lie in (0, 1)"); return 1; }
    for (int i = 0; i + 1 < v->nreg - 1; i++) if (!(v->breaks[i] < v->breaks[i + 1])) { snprintf(err, errlen, "volpath: the break dates must be strictly increasing"); return 1; }
    if (v->omega) for (int i = 0; i < v->nreg; i++) if (!(isfinite(v->omega[i]) && v->omega[i] > 0)) { snprintf(err, errlen, "volpath: the standard deviations omega must be finite and positive"); return 1; }
    return 0;
}
/* ---------- the local null values of the two applications of a normal component ---------- */
int dfl_cadf_c0(double c, double rho, double R, double *c0, char *err, size_t errlen)
{
    if (isnan(R)) R = fabs(rho);
    if (!(isfinite(rho) && rho >= -1 && rho <= 1 && rho != 0)) { snprintf(err, errlen, "cadf_c0: rho must lie in [-1, 1] and not be 0"); return 1; }
    if (!(isfinite(R) && R > 0)) { snprintf(err, errlen, "cadf_c0: R must be finite and positive"); return 1; }
    if (!isfinite(c)) { snprintf(err, errlen, "cadf_c0: c must be finite"); return 1; }
    *c0 = c == 0 ? 0.0 : c * (1 - 1 / (rho * R));
    return 0;
}
int dfl_predictive_c0(double c, double rho, double b, double *c0, char *err, size_t errlen)
{
    if (!(isfinite(rho) && rho >= -1 && rho <= 1 && rho != 0)) { snprintf(err, errlen, "predictive_c0: rho must lie in [-1, 1] and not be 0"); return 1; }
    if (!(isfinite(c) && isfinite(b))) { snprintf(err, errlen, "predictive_c0: c and b must be finite"); return 1; }
    *c0 = c - b / rho;
    return 0;
}

/* (stubs, filled in below by the later stages) */
int det_cov_dist(int kind, const dfl_design *D, int n, const double *x, double c, const dfl_det_opts *o, double *out, char *err, size_t errlen)
{ (void) kind; (void) D; (void) n; (void) x; (void) c; (void) o; (void) out; snprintf(err, errlen, "rho: not yet"); return 1; }
int det_cov_quantile(const dfl_design *D, int n, const double *p, double c, double x0, double tol, const dfl_det_opts *o, double *q, char *err, size_t errlen)
{ (void) D; (void) n; (void) p; (void) c; (void) x0; (void) tol; (void) o; (void) q; snprintf(err, errlen, "rho: not yet"); return 1; }
int dfl_hegy_t3(int kind, const dfl_design *D, int n, const double *x, double c, const dfl_det_opts *o, double *out, char *err, size_t errlen)
{ (void) kind; (void) D; (void) n; (void) x; (void) c; (void) o; (void) out; snprintf(err, errlen, "t3: not yet"); return 1; }
int dfl_hegy_t3_quantile(const dfl_design *D, int n, const double *p, double c, double x0, double tol, const dfl_det_opts *o, double *q, char *err, size_t errlen)
{ (void) D; (void) n; (void) p; (void) c; (void) x0; (void) tol; (void) o; (void) q; snprintf(err, errlen, "t3: not yet"); return 1; }
