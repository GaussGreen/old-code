
#include "SRT_H_ALL.H>
#define NRANSI
#define MAXITS 200
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0
#define TINY 1.0e-20;

static void lubksb(double **a, int n, int *indx, double b[]) {
  int i, ii = 0, ip, j;
  double sum;

  for (i = 1; i <= n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n; i >= 1; i--) {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software #>2$. */

static void ludcmp(double **a, int n, int *indx, double *d) {
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = dvector(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i][j])) > big)
        big = temp;
    if (big == 0.0)
      serror("Singular matrix in routine ludcmp");
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = a[i][j];
      for (k = 1; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = a[i][j];
      for (k = 1; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 1; k <= n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i <= n; i++)
        a[i][j] *= dum;
    }
  }
  free_vector(vv, 1, n);
}

int nn;
double *fvec;
void (*nrfuncv)(int n, double v[], double f[]);
#define FREERETURN                                                             \
  {                                                                            \
    free_dvector(fvec, 1, n);                                                  \
    free_dvector(xold, 1, n);                                                  \
    free_dvector(p, 1, n);                                                     \
    free_dvector(g, 1, n);                                                     \
    free_matrix(fjac, 1, n, 1, n);                                             \
    free_ivector(indx, 1, n);                                                  \
    return;                                                                    \
  }

void newt(double x[], int n, int *check,
          void (*vecfunc)(int, double[], double[])) {
  void fdjac(int n, double x[], double fvec[], double **df,
             void (*vecfunc)(int, double[], double[]));
  double fmin(double x[]);
  void lnsrch(int n, double xold[], double fold, double g[], double p[],
              double x[], double *f, double stpmax, int *check,
              double (*func)(double[]));
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  int i, its, j, *indx;
  double d, den, f, fold, stpmax, sum, temp, test, **fjac, *g, *p, *xold;

  indx = ivector(1, n);
  fjac = dmatrix(1, n, 1, n);
  g = dvector(1, n);
  p = dvector(1, n);
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
  for (its = 1; its <= MAXITS; its++) {
    fdjac(n, x, fvec, fjac, vecfunc);
    for (i = 1; i <= n; i++) {
      for (sum = 0.0, j = 1; j <= n; j++)
        sum += fjac[j][i] * fvec[j];
      g[i] = sum;
    }
    for (i = 1; i <= n; i++)
      xold[i] = x[i];
    fold = f;
    for (i = 1; i <= n; i++)
      p[i] = -fvec[i];
    ludcmp(fjac, n, indx, &d);
    lubksb(fjac, n, indx, p);
    lnsrch(n, xold, fold, g, p, x, &f, stpmax, check, fmin);
    test = 0.0;
    for (i = 1; i <= n; i++)
      if (fabs(fvec[i]) > test)
        test = fabs(fvec[i]);
    if (test < TOLF) {
      *check = 0;
      FREERETURN
    }
    if (*check) {
      test = 0.0;
      den = FMAX(f, 0.5 * n);
      for (i = 1; i <= n; i++) {
        temp = fabs(g[i]) * FMAX(fabs(x[i]), 1.0) / den;
        if (temp > test)
          test = temp;
      }
      *check = (test < TOLMIN ? 1 : 0);
      FREERETURN
    }
    test = 0.0;
    for (i = 1; i <= n; i++) {
      temp = (fabs(x[i] - xold[i])) / FMAX(fabs(x[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX)
      FREERETURN
  }
  serror("MAXITS exceeded in newt");
}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
#undef TINY
/* (C) Copr. 1986-92 Numerical Recipes Software #>2$. */
