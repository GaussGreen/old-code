/* =============================================================
   FILENAME:              num_f_broydn.c
   ============================================================= */
#include "math.h"
#include "utallhdr.h"
#define NRANSI
#define MAXITS 200
#define EPS_BROYDN 1.0e-8
#define TOLF 1.0E-4
#define TOLX 1.0E-4
#define STPMX 100.0
#define TOLMIN 1.0e-6
#define FACTOR 1.6
#define NTRY 50
#define FREERETURN                                                             \
  {                                                                            \
    free_dvector(fvec, 1, n);                                                  \
    free_dvector(xold, 1, n);                                                  \
    free_dvector(w, 1, n);                                                     \
    free_dvector(t, 1, n);                                                     \
    free_dvector(s, 1, n);                                                     \
    free_dmatrix(r, 1, n, 1, n);                                               \
    free_dmatrix(qt, 1, n, 1, n);                                              \
    free_dvector(p, 1, n);                                                     \
    free_dvector(g, 1, n);                                                     \
    free_dvector(fvcold, 1, n);                                                \
    free_dvector(d, 1, n);                                                     \
    free_dvector(c, 1, n);                                                     \
    return (error);                                                            \
  }

EXTERND int nn;
EXTERND double *fvec;
EXTERND void (*nrfuncv)(int n, double v[], double f[]);

char *broydn(double x[], int n, int *check,
             void (*vecfunc)(int, double[], double[])) {
  void fdjac(int n, double x[], double fvec[], double **df,
             void (*vecfunc)(int, double[], double[]));
  double fmin(double x[]);
  char *lnsrch(int n, double xold[], double fold, double g[], double p[],
               double x[], double *f, double stpmax, int *check,
               double (*func)(double[]));
  void qrdcmp(double **a, int n, double *c, double *d, int *sing);
  void qrupdt(double **r, double **qt, int n, double u[], double v[]);
  void rsolv(double **a, int n, double d[], double b[]);
  int i, its, j, k, restrt, sing, skip;
  double den, f, fold, stpmax, sum, temp, test, *c, *d, *fvcold;
  double *g, *p, **qt, **r, *s, *t, *w, *xold;
  char *error = 0;

  c = dvector(1, n);
  d = dvector(1, n);
  fvcold = dvector(1, n);
  g = dvector(1, n);
  p = dvector(1, n);
  qt = dmatrix(1, n, 1, n);
  r = dmatrix(1, n, 1, n);
  s = dvector(1, n);
  t = dvector(1, n);
  w = dvector(1, n);
  xold = dvector(1, n);
  fvec = dvector(1, n);
  nn = n;
  nrfuncv = vecfunc;
  f = fmin(x);
  test = 0.0;
  for (i = 1; i <= n; i++)
    if (fabs(fvec[i]) > test)
      test = fabs(fvec[i]);
  if (test < 0.01 * TOLF)
    FREERETURN
  for (sum = 0.0, i = 1; i <= n; i++)
    sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt(sum), (double)n);
  restrt = 1;
  for (its = 1; its <= MAXITS; its++) {
    if (restrt) {
      fdjac(n, x, fvec, r, vecfunc);
      qrdcmp(r, n, c, d, &sing);
      if (sing) {
        *check = -1;
      }
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++)
          qt[i][j] = 0.0;
        qt[i][i] = 1.0;
      }
      for (k = 1; k < n; k++) {
        if (c[k]) {
          for (j = 1; j <= n; j++) {
            sum = 0.0;
            for (i = k; i <= n; i++)
              sum += r[i][k] * qt[i][j];
            sum /= c[k];
            for (i = k; i <= n; i++)
              qt[i][j] -= sum * r[i][k];
          }
        }
      }
      for (i = 1; i <= n; i++) {
        r[i][i] = d[i];
        for (j = 1; j < i; j++)
          r[i][j] = 0.0;
      }
    } else {
      for (i = 1; i <= n; i++)
        s[i] = x[i] - xold[i];
      for (i = 1; i <= n; i++) {
        for (sum = 0.0, j = i; j <= n; j++)
          sum += r[i][j] * s[j];
        t[i] = sum;
      }
      skip = 1;
      for (i = 1; i <= n; i++) {
        for (sum = 0.0, j = 1; j <= n; j++)
          sum += qt[j][i] * t[j];
        w[i] = fvec[i] - fvcold[i] - sum;
        if (fabs(w[i]) >= EPS_BROYDN * (fabs(fvec[i]) + fabs(fvcold[i])))
          skip = 0;
        else
          w[i] = 0.0;
      }
      if (!skip) {
        for (i = 1; i <= n; i++) {
          for (sum = 0.0, j = 1; j <= n; j++)
            sum += qt[i][j] * w[j];
          t[i] = sum;
        }
        for (den = 0.0, i = 1; i <= n; i++)
          den += SQR(s[i]);
        for (i = 1; i <= n; i++)
          s[i] /= den;
        qrupdt(r, qt, n, t, s);
        for (i = 1; i <= n; i++) {
          if (r[i][i] == 0.0) {
            error = "Error: R is singular in broydn.";
            FREERETURN
          };
          d[i] = r[i][i];
        }
      }
    }
    for (i = 1; i <= n; i++) {
      for (sum = 0.0, j = 1; j <= n; j++)
        sum += qt[i][j] * fvec[j];
      g[i] = sum;
    }
    for (i = n; i >= 1; i--) {
      for (sum = 0.0, j = 1; j <= i; j++)
        sum += r[j][i] * g[j];
      g[i] = sum;
    }
    for (i = 1; i <= n; i++) {
      xold[i] = x[i];
      fvcold[i] = fvec[i];
    }
    fold = f;
    for (i = 1; i <= n; i++) {
      for (sum = 0.0, j = 1; j <= n; j++)
        sum += qt[i][j] * fvec[j];
      p[i] = -sum;
    }
    rsolv(r, n, d, p);
    error = lnsrch(n, xold, fold, g, p, x, &f, stpmax, check, fmin);
    if (error)
      FREERETURN
    test = 0.0;
    for (i = 1; i <= n; i++)
      if (fabs(fvec[i]) > test)
        test = fabs(fvec[i]);
    if (test < TOLF) {
      *check = 0;
      FREERETURN
    }
    if (*check) {
      if (restrt)
        FREERETURN
      else {
        test = 0.0;
        den = FMAX(f, 0.5 * n);
        for (i = 1; i <= n; i++) {
          temp = fabs(g[i]) * FMAX(fabs(x[i]), 1.0) / den;
          if (temp > test)
            test = temp;
        }
        if (test < TOLMIN)
          FREERETURN
        else
          restrt = 1;
      }
    } else {
      restrt = 0;
      test = 0.0;
      for (i = 1; i <= n; i++) {
        temp = (fabs(x[i] - xold[i])) / FMAX(fabs(x[i]), 1.0);
        if (temp > test)
          test = temp;
      }
      if (test < TOLX)
        FREERETURN
    }
  }
  error = "Error: MAXITS exceeded in broydn.";
  FREERETURN
}

int zbrac(double (*func)(double), double *x1, double *x2) {
  void nrerror(char error_text[]);
  int j;
  double f1, f2;

  if (*x1 == *x2)
    serror("Bad initial range in zbrac");
  f1 = (*func)(*x1);
  f2 = (*func)(*x2);
  for (j = 1; j <= NTRY; j++) {
    if (f1 * f2 < 0.0)
      return 1;
    if (fabs(f1) < fabs(f2))
      f1 = (*func)(*x1 += FACTOR * (*x1 - *x2));
    else
      f2 = (*func)(*x2 += FACTOR * (*x2 - *x1));
  }
  return 0;
}
#undef FACTOR
#undef NTRY
#undef MAXITS
#undef EPS_BROYDN
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
