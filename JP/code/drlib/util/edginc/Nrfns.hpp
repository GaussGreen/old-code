//----------------------------------------------------------------------------
//
//      File           : Nrfns.cpp
//
//      Description    : Headers for some Numerical Recipes functions.
//						Use arrays with ptr-1.
//
//
//----------------------------------------------------------------------------

#ifndef NRFNS_H
#define NRFNS_H

#include "edginc/NRException.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

bool causeIsNR(exception* e);

UTIL_DLL float gasdev(long *idum);  // only this is kept float

UTIL_DLL double gammln(double xx);
UTIL_DLL void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
UTIL_DLL void polcof(double xa[], double ya[], int n, double cof[]);
//// Versions using iterators - needed for structer ansi compliance
#if defined(__GNUC__) && ( __GNUC__ >= 3)
void spline(vector<double>::iterator x,
            vector<double>::iterator y, 
            int                      n, 
            double                   yp1,
            double                   ypn, 
            vector<double>::iterator y2);

void splint(vector<double>::iterator xa, 
            vector<double>::iterator ya, 
            vector<double>::iterator y2a,
            int n, 
            double x,
            double *y);
#endif
UTIL_DLL void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
UTIL_DLL void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
UTIL_DLL void hpsort(unsigned long n, double ra[]);
UTIL_DLL double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2, double xacc);
UTIL_DLL double zbrent(double (*func)(double), double x1, double x2, double tol);
UTIL_DLL double zbrentUseful(double (*func)(double , void *), void *params, double x1, double x2, double tol);
UTIL_DLL double zbrentUsefulBoundary(double (*func)(double , void *), void *params, double x1, double x2, double tol, double fx1, double fx2);
UTIL_DLL double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
UTIL_DLL int zbrac(double (*func)(double), double *x1, double *x2);
typedef enum {ZBRAC_SUCCESS, ZBRAC_BRAC_ERROR, ZBRAC_FUNC_ERROR} ZbracReturn;
UTIL_DLL ZbracReturn zbracUseful(double (*func)(double, void*), void *params, double *x1, double *x2, double *fx1, double *fx2);
UTIL_DLL double nrselect(unsigned long k, unsigned long n, double arr[]);
UTIL_DLL void dfpmin(double p[], int n, double gtol, int *iter, double *fret, double(*func)(double []), void (*dfunc)(double [], double []));
UTIL_DLL void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double []));
UTIL_DLL void locate(double xx[], unsigned long n, double x, unsigned long *j);
UTIL_DLL void hunt(double xx[], unsigned long n, double x, unsigned long *jlo);
UTIL_DLL void jacobi(double **a, int n, double d[], double **v, int *nrot);
UTIL_DLL void eigsrt(double d[], double **v, int n);
UTIL_DLL void svdcmp(double **a, int m, int n, double w[], double **v);
UTIL_DLL void choldc(double **a, int n, double p[]);
UTIL_DLL void ludcmp(double **a, int n, int *indx, double *d);
UTIL_DLL void lubksb(double **a, int n, int *indx, double b[]);
UTIL_DLL void fdjac(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []));
UTIL_DLL void mnewt(int ntrial, double x[], int n, double tolx, double tolf,
           int* k, double* errx, double* errf,
           void (*usrfun)(double *x,int n,double *fvec,double **fjac));
/* Here MAXITS is the maximum number of iterations; TOLF sets the convergence criterion on
function values; TOLMIN sets the criterion for deciding whether spurious convergence to a
minimum of fmin has occurred; TOLX is the convergence criterion on äx; STPMX is the scaled
maximum step length allowed in line searches. */
#define MAXITS_default 200
#define TOLF_default 1.0e-4
#define TOLMIN_default 1.0e-6
#define TOLX_default 1.0e-7
#define STPMX_default 100.0
UTIL_DLL void newt(double x[], int n,
          int *its, double *errf, double *errx, int *check,
          void (*vecfunc)(int, double [], double []),
          void (*jac)(int, double [], double [], double **,
                      void (*)(int, double [], double [])) = fdjac,
          int MAXITS = MAXITS_default, double TOLF = TOLF_default,
          double TOLMIN = TOLMIN_default, double TOLX = TOLX_default,
          double STPMX = STPMX_default);
#undef MAXITS_default
#undef TOLF_default
#undef TOLMIN_default
#undef TOLX_default
#undef STPMX_default
UTIL_DLL void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double));

/** Given mm, kk, and cof[0..mm+kk], evaluate and return the rational function (cof[0] +
cof[1]x + · · · + cof[mm]xmm)/(1 + cof[mm+1]x + · · · + cof[mm+kk]xkk). */
UTIL_DLL double ratval(double x, double cof[], int mm, int kk);

/** Pade approximation */
UTIL_DLL void pade(double cof[], int n, double *resid);

/** Improves a solution vector x[1..n] of the linear set of equations A · X = B. */
UTIL_DLL void mprove(float **a, float **alud, int n, int indx[], float b[], float x[]);

//// undocumented
UTIL_DLL float ran1(long *idum);
DRLIB_END_NAMESPACE

#endif
