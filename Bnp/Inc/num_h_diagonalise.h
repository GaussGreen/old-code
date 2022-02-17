/* ========================================================================

  MODULE:	    num_h_diagonalise.h

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

#ifndef NUM_H_DIAGONALISE_H
#define NUM_H_DIAGONALISE_H

Err diagonalise_symmetric_matrix(double **sym_mat, int n, double eigen_val[],
                                 double **eigen_vec);

#endif