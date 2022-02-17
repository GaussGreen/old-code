// S.Galluccio: 21 November 2000
// Created source file

#include "math.h"
#include "num_h_allhdr.h"

/*---------------------------------------------------------------------------------------------
Given a vector of n double values passed as a pointer *y  , this function
returns the same function smoothed. The smoothing degree is given by the
parameter nw. Higher nw correspond to a smoother result
-----------------------------------------------------------------------------------------------------*/

void FunctionSmooth(double *y, long int n, long int nw)

{
  long int i, j;
  double wt = 1.0, *Weig, x, *ysm, Sw;

  Weig = dvector(-nw, nw);
  ysm = dvector(1, n);

  wt = 2.0;
  Sw = 0.0;
  for (i = -nw; i <= nw; i++) {
    x = i;
    Weig[i] = 1. / (fabs(x) + 1.) * wt;
    Sw = Sw + Weig[i];
  }
  for (i = nw + 1; i <= n - nw; i++) {
    ysm[i] = 0.0;
    for (j = -nw; j <= nw; j++)
      ysm[i] += 1. / Sw * Weig[j] * y[i + j];
  }
  for (i = 1; i <= nw; i++) {
    ysm[i] = 0.0;
    Sw = 0.0;
    for (j = -i + 1; j <= nw; j++) {
      ysm[i] += Weig[j] * y[i + j];
      Sw = Sw + Weig[j];
    }
    ysm[i] = ysm[i] / Sw;
  }
  for (i = n - nw; i <= n; i++) {
    ysm[i] = 0.0;
    Sw = 0.0;
    for (j = -nw; j <= n - i; j++) {
      ysm[i] += Weig[j] * y[i + j];
      Sw = Sw + Weig[j];
    }
    ysm[i] = ysm[i] / Sw;
  }

  for (i = 1; i <= n; i++)
    y[i] = ysm[i];

  free_dvector(Weig, -nw, nw);
  free_dvector(ysm, 1, n);
}