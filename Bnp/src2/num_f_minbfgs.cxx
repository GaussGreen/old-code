/* -------------------------------------------------------------

   AUTHOR: E. FOURNIE

   DATE : MAY 99

   FILENAME: num_f_minbfgs.cxx

   PURPOSE:	Quasi-Newton and Conjugate gradient as in numerical recipes

   MODIFICATION: JAN 00

 ----------------------------------------------------------------------- */
#include "math.h"
#include "num_h_minbfgs.h"
#include "srt_h_all.h"

#define NRANSI
#define ITMAX 200
#define EPSILOC 3.0e-8
#define TOLX (4 * EPSILOC)
#define STPMX 100.0
#define ALF 1.0e-4

#define FFMAX(a, b) ((a) > (b) ? (a) : (b))
#define FREEALL                                                                \
  free_dvector(xi, 1, n);                                                      \
  free_dvector(pnew, 1, n);                                                    \
  free_matrix(hessin, 1, n, 1, n);                                             \
  free_dvector(hdg, 1, n);                                                     \
  free_dvector(g, 1, n);                                                       \
  free_dvector(dg, 1, n);

Err dfpmin(DFPMinCriteriumType returnCriterium, double p[], int n,
           double func_gtol, double dfunc_tol, int *iter, double *fret,
           Err (*func)(double[], double *), Err (*dfunc)(double[], double[])) {
  Err err = NULL;
  Err lnsrch1(int n, double xold[], double fold, double g[], double p[],
              double x[], double *f, double stpmax, int *check,
              Err (*func)(double[], double *));

  int check, i, its, j;
  double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test;
  double *dg, *g, *hdg, **hessin, *pnew, *xi;

  dg = dvector(1, n);
  g = dvector(1, n);
  hdg = dvector(1, n);
  hessin = matrix(1, n, 1, n);
  pnew = dvector(1, n);
  xi = dvector(1, n);

  err = (*func)(p, &fp);
  if (err)
    return err;

  err = (*dfunc)(p, g);
  if (err)
    return err;

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++)
      hessin[i][j] = 0.0;
    hessin[i][i] = 1.0;
    xi[i] = -g[i];
    sum += p[i] * p[i];
  }
  stpmax = STPMX * FFMAX(sqrt(sum), (double)n);
  for (its = 1; its <= ITMAX; its++) {
    *iter = its;
    err = lnsrch1(n, p, fp, g, xi, pnew, fret, stpmax, &check, func);
    if (err)
      return err;
    /* stop on func close to zero (gtol) */
    if (returnCriterium == FUNCCRITERIUM)
      if (fabs(*fret) < func_gtol)
        return err;
    /* else criterium zero gradient */
    fp = *fret;
    for (i = 1; i <= n; i++) {
      xi[i] = pnew[i] - p[i];
      p[i] = pnew[i];
    }
    test = 0.0;
    for (i = 1; i <= n; i++) {
      temp = fabs(xi[i]) / FFMAX(fabs(p[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX) {
      FREEALL
      return err;
    }
    for (i = 1; i <= n; i++)
      dg[i] = g[i];
    err = (*dfunc)(p, g);
    if (err)
      return err;
    test = 0.0;
    den = FFMAX(*fret, 1.0);
    for (i = 1; i <= n; i++) {
      temp = fabs(g[i]) * FFMAX(fabs(p[i]), 1.0) / den;
      if (temp > test)
        test = temp;
    }
    if (test < dfunc_tol) {
      FREEALL
      return err;
    }
    for (i = 1; i <= n; i++)
      dg[i] = g[i] - dg[i];
    for (i = 1; i <= n; i++) {
      hdg[i] = 0.0;
      for (j = 1; j <= n; j++)
        hdg[i] += hessin[i][j] * dg[j];
    }
    fac = fae = sumdg = sumxi = 0.0;
    for (i = 1; i <= n; i++) {
      fac += dg[i] * xi[i];
      fae += dg[i] * hdg[i];
      sumdg += SQR(dg[i]);
      sumxi += SQR(xi[i]);
    }
    if (fac * fac > EPSILOC * sumdg * sumxi) {
      fac = 1.0 / fac;
      fad = 1.0 / fae;
      for (i = 1; i <= n; i++)
        dg[i] = fac * xi[i] - fad * hdg[i];
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
          hessin[i][j] +=
              fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
        }
      }
    }
    for (i = 1; i <= n; i++) {
      xi[i] = 0.0;
      for (j = 1; j <= n; j++)
        xi[i] -= hessin[i][j] * g[j];
    }
  }
  smessage("too many iterations in dfpmin");
  FREEALL

  return err;
}

static Err lnsrch1(int n, double xold[], double fold, double g[], double p[],
                   double x[], double *f, double stpmax, int *check,
                   Err (*func)(double[], double *)) {
  Err err = NULL;
  int i;
  double a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope, sum,
      temp, test, tmplam;

  *check = 0;
  for (sum = 0.0, i = 1; i <= n; i++)
    sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++)
      p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  test = 0.0;
  for (i = 1; i <= n; i++) {
    temp = fabs(p[i]) / FFMAX(fabs(xold[i]), 1.0);
    if (temp > test)
      test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;) {
    for (i = 1; i <= n; i++)
      x[i] = xold[i] + alam * p[i];
    err = (*func)(x, f);
    if (err)
      return err;
    if (alam < alamin) {
      for (i = 1; i <= n; i++)
        x[i] = xold[i];
      *check = 1;
      return err;
    } else if (*f <= fold + ALF * alam * slope)
      return err;
    else {
      if (alam == 1.0)
        tmplam = -slope / (2.0 * (*f - fold - slope));
      else {
        rhs1 = *f - fold - alam * slope;
        rhs2 = f2 - fold2 - alam2 * slope;
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) /
            (alam - alam2);
        if (a == 0.0)
          tmplam = -slope / (2.0 * b);
        else {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0.0)
            smessage("Roundoff problem in lnsrch.");
          else
            tmplam = (-b + sqrt(disc)) / (3.0 * a);
        }
        if (tmplam > 0.5 * alam)
          tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = *f;
    fold2 = fold;
    alam = FFMAX(tmplam, 0.1 * alam);
  }
}

#undef ITMAX
#undef EPSILOC
#undef TOLX
#undef STPMX
#undef FREEALL
#undef ALF
#undef TOLX
#undef FFMAX
#undef NRANSI

/* ------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------
 */

#define NRANSI
#define ITMAX 200
#define EPSILOC 1.0e-10
#define FREEALLFRPMN                                                           \
  free_dvector(xi, 1, n);                                                      \
  free_dvector(h, 1, n);                                                       \
  free_dvector(g, 1, n);

Err frprmn(double p[], int n, double ftol, int *iter, double *fret,
           Err (*func)(double[], double *), Err (*dfunc)(double[], double[])) {
  void linmin(double p[], double xi[], int n, double *fret,
              Err (*func)(double[], double *));
  int j, its;
  double gg, gam, fp, dgg;
  double *g, *h, *xi;
  Err err = NULL;

  g = dvector(1, n);
  h = dvector(1, n);
  xi = dvector(1, n);
  err = (*func)(p, &fp);
  if (err)
    return err;

  err = (*dfunc)(p, xi);
  if (err)
    return err;

  for (j = 1; j <= n; j++) {
    g[j] = -xi[j];
    xi[j] = h[j] = g[j];
  }
  for (its = 1; its <= ITMAX; its++) {
    *iter = its;
    linmin(p, xi, n, fret, func);
    if (2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + EPSILOC)) {
      FREEALLFRPMN
      return err;
    }
    err = (*func)(p, &fp);
    if (err)
      return err;

    (*dfunc)(p, xi);
    dgg = gg = 0.0;
    for (j = 1; j <= n; j++) {
      gg += g[j] * g[j];
      dgg += (xi[j] + g[j]) * xi[j];
    }
    if (gg == 0.0) {
      FREEALLFRPMN
      return err;
    }
    gam = dgg / gg;
    for (j = 1; j <= n; j++) {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j] + gam * h[j];
    }
  }
  smessage("Too many iterations in frprmn");

  return err;
}
#undef ITMAX
#undef EPSILOC
#undef FREEALL

/* ------------------------------------------------------------------------------------
 */

#define TOL 2.0e-4
int ncom;
double *pcom, *xicom;
Err (*nrfunc)(double[], double *);

static void linmin(double p[], double xi[], int n, double *fret,
                   Err (*func)(double[], double *)) {
  double brent(double ax, double bx, double cx, double (*f)(double), double tol,
               double *xmin);
  double f1dim(double x);
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
              double *fc, double (*func)(double));
  int j;
  double xx, xmin, fx, fb, fa, bx, ax;

  ncom = n;
  pcom = dvector(1, n);
  xicom = dvector(1, n);
  nrfunc = func;

  for (j = 1; j <= n; j++) {
    pcom[j] = p[j];
    xicom[j] = xi[j];
  }
  ax = 0.0;
  xx = 1.0;
  mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
  *fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
  for (j = 1; j <= n; j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom, 1, n);
  free_dvector(pcom, 1, n);
}
#undef TOL

static double f1dim(double x) {
  int j;
  double f, *xt;
  Err err = NULL;

  xt = dvector(1, ncom);
  for (j = 1; j <= ncom; j++)
    xt[j] = pcom[j] + x * xicom[j];
  err = (*nrfunc)(xt, &f);
  free_dvector(xt, 1, ncom);
  return f;
}

/* ------------------------------------------------------------------------------------
 */

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a, b, c, d)                                                       \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);

static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
                   double *fc, double (*func)(double)) {
  double ulim, u, r, q, fu, dum;

  *fa = (*func)(*ax);
  *fb = (*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum, *ax, *bx, dum)
    SHFT(dum, *fb, *fa, dum)
  }
  *cx = (*bx) + GOLD * (*bx - *ax);
  *fc = (*func)(*cx);
  while (*fb > *fc) {
    r = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
                    (2.0 * SIGN(FMAX(fabs(q - r), TINY), q - r));
    ulim = (*bx) + GLIMIT * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0) {
      fu = (*func)(u);
      if (fu < *fc) {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
        return;
      } else if (fu > *fb) {
        *cx = u;
        *fc = fu;
        return;
      }
      u = (*cx) + GOLD * (*cx - *bx);
      fu = (*func)(u);
    } else if ((*cx - u) * (u - ulim) > 0.0) {
      fu = (*func)(u);
      if (fu < *fc) {
        SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
        SHFT(*fb, *fc, fu, (*func)(u))
      }
    } else if ((u - ulim) * (ulim - *cx) >= 0.0) {
      u = ulim;
      fu = (*func)(u);
    } else {
      u = (*cx) + GOLD * (*cx - *bx);
      fu = (*func)(u);
    }
    SHFT(*ax, *bx, *cx, u)
    SHFT(*fa, *fb, *fc, fu)
  }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

/* ------------------------------------------------------------------------------------
 */

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPSILOC 1.0e-10
#define SHFT(a, b, c, d)                                                       \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);

static double brent(double ax, double bx, double cx, double (*f)(double),
                    double tol, double *xmin) {
  int iter;
  double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
  double e = 0.0;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = (*f)(x);
  for (iter = 1; iter <= ITMAX; iter++) {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPSILOC);
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) ||
          p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));
      else {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      }
    } else {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }
    u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
    fu = (*f)(u);
    if (fu <= fx) {
      if (u >= x)
        a = x;
      else
        b = x;
      SHFT(v, w, x, u)
      SHFT(fv, fw, fx, fu)
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  smessage("Too many iterations in brent");
  *xmin = x;
  return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPSILOC
#undef SHFT
#undef NRANSI

/* ========= END OF FILE =================================================== */
