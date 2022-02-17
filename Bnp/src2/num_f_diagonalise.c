/* ========================================================================

  MODULE:	    num_f_diagonalise.c

  LIBRARY:	    NUM_LIB

  FUNCTION:	    diagonalise_symmetric_matrix

  DESCRIPTION:	From Numerical Recipes in C chapter 11
                The most performent way to perform the diagonalisation of
                                a real symmetric matrix of large size  , based
on Householder tridiagonalisation method  , followed by a tridiagonal QL
                                implicit method
                In input  , sym_mat[0..n-1][0..n-1] is a real symmetric matrix
                On output  , eigen_val[0..n-1] returns the eigen values of
sym_mat. eigen_vec[0..n-1][0..n-1] is a matrix that stores all the corresponding
renormalised eigen vectors of sym_mat ( eigen_vec[...][i] corresponds to
eigen_val[i] )

========================================================================== */

#include "num_h_allhdr.h"

static Err sort_eigen_values(double eigen_val[], double **eigen_vec, int n)
/* From NUMERICAL RECIPES IN C
   Given the dimension of the matrix and its eigen_values  ,
   sorts the eigen values in the eigen_val vector  ,
   and replaces the eigen vectors in the right order in **eigen_vec*/
{
  int i, j, k;
  double p;
  Err err = NULL;

  for (i = 0; i < n - 1; i++) {
    p = eigen_val[k = i];
    for (j = i + 1; j < n; j++)
      if (eigen_val[j] >= p)
        p = eigen_val[k = j];
    if (k != i) {
      eigen_val[k] = eigen_val[i];
      eigen_val[i] = p;
      for (j = 0; j < n; j++) {
        p = eigen_vec[j][i];
        eigen_vec[j][i] = eigen_vec[j][k];
        eigen_vec[j][k] = p;
      }
    }
  }
  return err;
}

/* -----------------------------------------------------------------------------
 */

Err diagonalise_symmetric_matrix(double **sym_mat, int n, double eigen_val[],
                                 double **eigen_vec) {
  long i, j;
  double *subdiago;

  Err err = NULL;

  /* Copy sym_mat into eigen_vec */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      eigen_vec[i][j] = sym_mat[i][j];
  }

  /* Performs the Householder tridiagonalisation (the diago is output in eigen
   * val)*/
  subdiago = dvector(0, n - 1);
  err = householder_tridiagonalisation(eigen_vec, n, eigen_val, subdiago);
  if (err)
    return err;

  /* Diagonalise using the QL implicit method */
  err = tridiagonal_QL_implicit(eigen_val, subdiago, n, eigen_vec);
  if (err)
    return err;
  free_dvector(subdiago, 0, n - 1);

  /* Sort the eigen values and eigen vectors */
  err = sort_eigen_values(eigen_val, eigen_vec, n);
  if (err)
    return err;

  /* Return a success string */
  return NULL;

} /* END Err diagonalise_symmetric_matrix (...) */

/* ------------------------------------------------------------------------ */