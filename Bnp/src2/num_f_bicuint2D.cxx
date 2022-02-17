/************************************************************
 * NR implementation of a bicubic interpolation
 * Copied from NR book
 * C. Godart
 * Date: 24/01/2000
 ************************************************************/
/* Here we compute the c_ij such that:
   y = Sum(i=1        ,i=4)Sum(j=1        ,j=4) c_ij t^(i-1)u^(j-1)
   dy/dx1 = Sum(i=1        ,i=4)Sum(j=1        ,j=4) (i-1) c_ij t^(i-2)u^(j-1)
   dy/dx2 = Sum(i=1        ,i=4)Sum(j=1        ,j=4) (i-1) c_ij t^(i-1)u^(j-2)
   d2y/dx1dx2 = Sum(i=1        ,i=4)Sum(j=1        ,j=4) (i-1)(j-1) c_ij
   t^(i-2)u^(j-2)
   */
/**********************************************************************
  This is the numbering convention for all the function in this file:
    Cell: Counterclockwise from lower left.
                    4       3
                                 x     x
                   <---
                                  ----|
                                 x     x
                                1       2
 **************************************************************************/

#include "utallhdr.h"

void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
            double d2, double **c) {
  static int wt[16][16] = {
      1,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  1,  0, 0,  0,  0,  0,  0,  0, -3, 0,  0,  3,  0,  0,
      0,  0,  -2, 0,  0,  -1, 0, 0,  0,  0,  2,  0,  0, -2, 0,  0,  0,  0,  1,
      0,  0,  1,  0,  0,  0,  0, 0,  0,  0,  0,  1,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  1,  0,  0,
      0,  0,  0,  0,  0,  -3, 0, 0,  3,  0,  0,  0,  0, -2, 0,  0,  -1, 0,  0,
      0,  0,  2,  0,  0,  -2, 0, 0,  0,  0,  1,  0,  0, 1,  -3, 3,  0,  0,  -2,
      -1, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      -3, 3,  0,  0,  -2, -1, 0, 0,  9,  -9, 9,  -9, 6, 3,  -3, -6, 6,  -6, -3,
      3,  4,  2,  1,  2,  -6, 6, -6, 6,  -4, -2, 2,  4, -3, 3,  3,  -3, -2, -1,
      -1, -2, 2,  -2, 0,  0,  1, 1,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0, 2,  -2, 0,  0,  1,  1, 0,  0,  -6, 6,  -6, 6,
      -3, -3, 3,  3,  -4, 4,  2, -2, -2, -2, -1, -1, 4, -4, 4,  -4, 2,  2,  -2,
      -2, 2,  -2, -2, 2,  1,  1, 1,  1};
  int l, k, j, i;
  double xx, d1d2, cl[16], x[16];

  d1d2 = d1 * d2;
  for (i = 1; i <= 4; i++) { /* pack a temporary vector x */
    x[i - 1] = y[i];
    x[i + 3] = y1[i] * d1;
    x[i + 7] = y2[i] * d2;
    x[i + 11] = y12[i] * d1d2;
  }
  for (i = 0; i <= 15; i++) { /* matrix multiply by the stored table */
    xx = 0.0;
    for (k = 0; k <= 15; k++)
      xx += wt[i][k] * x[k];
    cl[i] = xx;
  }
  l = 0;
  for (i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      c[i][j] = cl[l++];
}

Err bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
           double x1u, double x2l, double x2u, double x1, double x2,
           double *ansy, double *ansy1, double *ansy2) {
  int i;
  double t, u, d1, d2, **c;

  c = dmatrix(1, 4, 1, 4);
  d1 = x1u - x1l;
  d2 = x2u - x2l;
  bcucof(y, y1, y2, y12, d1, d2, c); /* get the c's */
  /*if (x1u==x1l || x2u==x2l)
          return "Bad input in routine bcuint";*/
  /* !!!!!!!!!!!!!! WARNING: this is a modification
        to the original algo in NR. No guarantee.
            Pb: When x is on a grid coord the search  algo from NR
            returns x1u==x1l which make this algo fail */
  t = (x1u == x1l) ? 0 : (x1 - x1l) / d1;
  u = (x2u == x2l) ? 0 : (x2 - x2l) / d2;
  *ansy = (*ansy2) = (*ansy1) = 0.0;
  for (i = 4; i >= 1; i--) {
    *ansy = t * (*ansy) + ((c[i][4] * u + c[i][3]) * u + c[i][2]) * u + c[i][1];
    *ansy2 = t * (*ansy2) + (3.0 * c[i][4] * u + 2.0 * c[i][3]) * u + c[i][2];
    *ansy1 = u * (*ansy1) + (3.0 * c[4][i] * t + 2.0 * c[3][i]) * t + c[2][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
  free_dmatrix(c, 1, 4, 1, 4);
  return NULL;
}

/*
  !!!!!!  This is not a NR routine. !!!!!
Here we compute the gradient and the cross derivatives
of a tabulated function by centered second order formulae.
y is a matrix of n        ,m. x1 is the value on the 1st axis.
                      x2 is the value on the 2nd axis.
        Input: x1a dvector(1        ,m)
               x2a dvector(1        ,n)
                   ya  dmatrix(1        ,m        ,1        ,n)

    Output: y1a = dy/dx1 dmatrix(1        ,m        ,1        ,n)
                y2a = dy/dx2 dmatrix(1        ,m        ,1        ,n)
                        y12a= d2y/dx1dx2 dmatrix(1        ,m        ,1 ,n)
 */
Err grid_fd_12(double x1a[], double x2a[], double **ya, int m, int n,
               double **y1a, double **y2a, double **y12a) {
  int j, k;
  /*check the grid axis are in ascending order (stricly)*/
  for (j = 2; j <= m; j++)
    if (x1a[j] <= x1a[j - 1])
      return "Wrong values for first axis";
  for (j = 2; j <= n; j++)
    if (x2a[j] <= x2a[j - 1])
      return "Wrong values for second axis";
  /* Centered schemes where possible */
  for (j = 2; j <= m - 1; j++) {
    for (k = 2; k <= n - 1; k++) {
      y1a[j][k] = (ya[j + 1][k] - ya[j - 1][k]) / (x1a[j + 1] - x1a[j - 1]);
      y2a[j][k] = (ya[j][k + 1] - ya[j][k - 1]) / (x2a[k + 1] - x2a[k - 1]);
      y12a[j][k] = (ya[j + 1][k + 1] - ya[j + 1][k - 1] - ya[j - 1][k + 1] +
                    ya[j - 1][k - 1]) /
                   ((x1a[j + 1] - x1a[j - 1]) * (x2a[k + 1] - x2a[k - 1]));
    }
  }
  /* Now on the hedges        , first order upward or downward derivatives */
  for (j = 2; j <= m - 1; j++) {
    y1a[j][1] = (ya[j + 1][1] - ya[j - 1][1]) / (x1a[j + 1] - x1a[j - 1]);
    y1a[j][n] = (ya[j + 1][n] - ya[j - 1][n]) / (x1a[j + 1] - x1a[j - 1]);
    y2a[j][1] = (ya[j][2] - ya[j][1]) / (x2a[2] - x2a[1]);
    y2a[j][n] = (ya[j][n] - ya[j][n - 1]) / (x2a[n] - x2a[n - 1]);
    y12a[j][1] = (ya[j + 1][2] - ya[j + 1][1] - ya[j - 1][2] + ya[j - 1][1]) /
                 ((x1a[j + 1] - x1a[j - 1]) * (x2a[2] - x2a[1]));
    y12a[j][n] =
        (ya[j + 1][n] - ya[j + 1][n - 1] - ya[j - 1][n] + ya[j - 1][n - 1]) /
        ((x1a[j + 1] - x1a[j - 1]) * (x2a[n] - x2a[n - 1]));
  }
  for (k = 2; k <= n - 1; k++) {
    y1a[1][k] = (ya[2][k] - ya[1][k]) / (x1a[2] - x1a[1]);
    y1a[m][k] = (ya[m][k] - ya[m - 1][k]) / (x1a[m] - x1a[m - 1]);
    y2a[1][k] = (ya[1][k + 1] - ya[1][k - 1]) / (x2a[k + 1] - x2a[k - 1]);
    y2a[m][k] = (ya[m][k + 1] - ya[m][k - 1]) / (x2a[k + 1] - x2a[k - 1]);
    y12a[1][k] = (ya[2][k + 1] - ya[1][k + 1] - ya[2][k - 1] + ya[1][k - 1]) /
                 ((x1a[2] - x1a[1]) * (x2a[k + 1] - x2a[k - 1]));
    y12a[m][k] =
        (ya[m][k + 1] - ya[m - 1][k + 1] - ya[m][k - 1] + ya[m - 1][k - 1]) /
        ((x1a[m] - x1a[m - 1]) * (x2a[k + 1] - x2a[k - 1]));
  }

  /*                 and the corners                */
  /* First derivative */
  y1a[1][1] = (ya[2][1] - ya[1][1]) / (x1a[2] - x1a[1]);
  y1a[1][n] = (ya[2][n] - ya[1][n]) / (x1a[2] - x1a[1]);
  y1a[m][1] = (ya[m][1] - ya[m - 1][1]) / (x1a[m] - x1a[m - 1]);
  y1a[m][n] = (ya[m][n] - ya[m - 1][n]) / (x1a[m] - x1a[m - 1]);
  /* Second derivative */
  y2a[1][1] = (ya[1][2] - ya[1][1]) / (x2a[2] - x2a[1]);
  y2a[m][1] = (ya[m][2] - ya[m][1]) / (x2a[2] - x2a[1]);
  y2a[1][n] = (ya[1][n] - ya[1][n - 1]) / (x2a[n] - x2a[n - 1]);
  y2a[m][n] = (ya[m][n] - ya[m][n - 1]) / (x2a[n] - x2a[n - 1]);
  /* Mixed derivative */
  y12a[1][1] = (ya[2][2] - ya[2][1] - ya[1][2] + ya[1][1]) /
               ((x1a[2] - x1a[1]) * (x2a[2] - x2a[1]));
  y12a[m][1] = (ya[m][2] - ya[m - 1][1] - ya[1][2] + ya[1][1]) /
               ((x1a[m] - x1a[m - 1]) * (x2a[2] - x2a[1]));
  y12a[1][n] = (ya[2][n] - ya[2][n - 1] - ya[1][n] + ya[1][n - 1]) /
               ((x1a[2] - x1a[1]) * (x2a[n] - x2a[n - 1]));
  y12a[m][n] = (ya[m][n] - ya[m][n - 1] - ya[m - 1][n] + ya[m - 1][n - 1]) /
               ((x1a[m] - x1a[m - 1]) * (x2a[n] - x2a[n - 1]));
  return NULL;
}
