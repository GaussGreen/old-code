#include "UTALLHDR.H>
#include "math.h"
#include "num_h_allhdr.h"

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// algorithm for Choleski decomposition
///////////////////////////////////////////////////////////////////////////////

void CholDec(void (*CholAlg)(double **, const int, double *),
             double **corr_matr, double **sqrt_matr, const int n)

{

  int i, j;
  double *p = dvector(0, n - 1);

  for (i = 0; i <= n - 1; i++) {
    for (j = i; j <= n - 1; j++) {
      sqrt_matr[i][j] = sqrt_matr[j][i] = corr_matr[i][j];
    }
  }

  (*CholAlg)(sqrt_matr, n, p);

  for (i = 0; i <= n - 1; i++) {
    sqrt_matr[i][i] = p[i];
    for (j = i + 1; j <= n - 1; j++) {
      sqrt_matr[i][j] = 0.0;
    }
  }

  free_dvector(p, 0, n - 1);
}

/////////////////////////////////////////////////////////////////////////////////////////

void CholAlg(double **a, const int n, double *p)

{
  int i, j, k;
  double sum;

  for (i = 0; i <= n - 1; i++) {
    for (j = i; j <= n - 1; j++) {
      for (sum = a[i][j], k = i - 1; k >= 0; k--)
        sum -= a[i][k] * a[j][k];
      if (i == j) {
        if (sum <= 0.0) {
          //					EM_FAIL_TEXT("Fatal: (CholAlg)
          //correlation matrix is not positive definite");
          goto end;
        }
        p[i] = sqrt(sum);
      } else
        a[j][i] = sum / p[i];
    }
  }
end:;
}