/* ===============================================================
   FILENAME : num_f_matrix.c

   PURPOSE:   functions to do several operations on matrixes
   =============================================================== */

#include "num_h_jacobi.h"
#include "utallhdr.h"
#include <math.h"
#include <num_h_gaussj.h"

/* returns transpose (x) */
double **transpose_matrix(double **old_matrix, long old_rl, long old_rh,
                          long old_cl, long old_ch) {
  long new_rl = old_cl, new_rh = old_ch, new_cl = old_rl, new_ch = old_rh;

  long i, j;

  double **new_matrix = dmatrix(new_rl, new_rh, new_cl, new_ch);

  if (!new_matrix)
    return NULL;

  for (i = new_rl; i <= new_rh; i++)
    for (j = new_cl; j <= new_ch; j++)
      new_matrix[i][j] = old_matrix[j][i];

  return new_matrix;
}

/* returns xy */
double **product_matrix(double **x, long rxl, long rxh, long cxl, long cxh,
                        double **y, long ryl, long ryh, long cyl, long cyh) {
  long i, j, k;
  long rows = rxh - rxl + 1, columns = cyh - cyl + 1, common = cxh - cxl + 1;
  double sum;
  double **new_matrix = dmatrix(rxl, rxh, cyl, cyh);

  if (!new_matrix)
    return NULL;

  if (common != ryh - ryl + 1) {
    free_dmatrix(new_matrix, rxl, rxh, cyl, cyh);
    return NULL;
  }

  for (i = 0; i < rows; i++) {
    for (j = 0; j < columns; j++) {
      sum = 0.0;
      for (k = 0; k < common; k++) {
        sum += x[rxl + i][cxl + k] * y[ryl + k][cyl + j];
      }
      new_matrix[rxl + i][cyl + j] = sum;
    }
  }
  return new_matrix;
}

/* returns inv(x) */
double **inverse_matrix(double **x, long rl, long rh) {
  long i, n = rh - rl + 1;
  double **temp = dmatrix(1, n, 1, n), **b = dmatrix(1, n, 1, 1);

  Err err;

  if ((!x) || (!b))
    return NULL;

  for (i = 1; i <= n; i++) {
    memcpy(&(temp[i][1]), &(x[rl + i - 1][rl]), n * sizeof(double));
    b[i][1] = 1.0;
  }

  err = gaussj(temp, n, b, 1);

  if (err) {
    free_dmatrix(temp, 1, n, 1, n);
    free_dmatrix(b, 1, n, 1, 1);
    return NULL;
  }

  free_dmatrix(b, 1, n, 1, 1);

  return temp;
}

/*	Cholesky decomposition */
void nr_choldc(int n, /* 0.. n-1*/
               double **cov, double **chol) {
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      for (sum = cov[i][j], k = i - 1; k >= 0; k--) {
        sum -= cov[i][k] * cov[j][k];
      }
      if (i == j) {
        chol[i][i] = sqrt(sum);
      } else {
        cov[j][i] = chol[j][i] = sum / chol[i][i];
      }
    }
  }
}

/*	Cholesky decomposition with error checking */
Err choldc(int n, double **a, double **chol) {
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      for (sum = a[i][j], k = i - 1; k >= 0; k--)
        sum -= a[i][k] * a[j][k];
      if (i == j) {
        if (sum <= 0.0)
          return serror("Matrix is not positive definite in choldc");
        chol[i][i] = sqrt(sum);
      } else
        a[j][i] = chol[j][i] = sum / chol[i][i];
    }
  }
  return NULL;
}

/*	Inverse Lower Trinagular Matrix */
void inverse_lower_triangular_matrix(double **mat, int n, /*	0..n-1 */
                                     double **res) {
  int i, j, k;
  double sum;

  for (i = n - 1; i >= 0; i--) {
    res[i][i] = 1.0 / mat[i][i];
    for (j = i - 1; j >= 0; j--) {
      sum = 0.0;
      for (k = i; k > j; k--) {
        sum -= res[i][k] * mat[k][j];
      }
      res[i][j] = sum / mat[j][j];
    }
  }
}

/* Algebraic method to make a correlation matrix positive definite */

Err PositiveMatrix(double **Matrix, int n)

{
  Err err = NULL;
  int i, j, k, nrot;
  double *eigen_val, **eigen_vec, **AuxMatrix, Thresh = 0.0, TraceOld = 0.0,
                                               TraceNew = 0.0;

  eigen_val = dvector(0, n - 1);
  eigen_vec = dmatrix(0, n - 1, 0, n - 1);
  AuxMatrix = dmatrix(0, n - 1, 0, n - 1);

  err = jacobi_diagonalisation2(Matrix, n, eigen_val, eigen_vec, &nrot);

  for (i = n - 1; i >= 0; i--) {
    TraceOld += eigen_val[i];

    if (eigen_val[i] < 0.0) {
      eigen_val[i] = Thresh + 1.e-4;
      Thresh += 1.e-4;
    }
  }

  /* corrects for the trace */

  for (i = n - 1; i >= 0; i--) {
    TraceNew += eigen_val[i];
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      AuxMatrix[i][j] = 0.0;

      for (k = 0; k < n; k++) {

        AuxMatrix[i][j] += eigen_vec[i][k] * eigen_vec[j][k] * eigen_val[k] *
                           TraceOld / TraceNew;
      }
    }
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {

      Matrix[i][j] = AuxMatrix[i][j] / sqrt(AuxMatrix[i][i] * AuxMatrix[j][j]);
    }
  }

  free_dvector(eigen_val, 0, n - 1);
  free_dmatrix(eigen_vec, 0, n - 1, 0, n - 1);
  free_dmatrix(AuxMatrix, 0, n - 1, 0, n - 1);

  return NULL;
}