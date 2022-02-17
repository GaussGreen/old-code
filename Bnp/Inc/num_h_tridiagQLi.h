/* ========================================================================

  MODULE:	    num_h_tridiagQLi.h

  LIBRARY:	    UTL_LIB

  FUNCTION:	    tridiagonal_QL_implicit

  INPUTS:


  DESCRIPTION:	From Numerical Recipes in C page 480
                QL algorithm with implicit shifts      , to determine the eigen
                                values and eigen vectors of a real      ,
symmetric , tridiagonal matrix      , previously reduced by the Householder
method (cf uthouseholder) On input      , d[0..n-1][0..n-1] contains the
diagonal elements of a tridiagonal matrix. On output      , it returns the iegen
values. The vector e[0..n-1] inputs the sub diagonal elements of the matrix ,
with e[0] arbitrary. On output      , e[...] is destroyed If the eigenvectors of
a tridiagonal matrix are desired      , the matrix z[0..n-1][0..n-1] is input as
the identity matrix. If the eigen vectors of a matrix that has been reduced by
                                Householder method
(householder_tridiagonalisation) are required      , than z is input as the
matrix output by this function. The kth column of z (i.e. z[...][k] returns the
                                normalised eigenvector corresponding to d[k]).

========================================================================== */

#ifndef NUM_H_TRIDIAGQLI_H
#define NUM_H_TRIDIAGQLI_H

Err tridiagonal_QL_implicit(double d[], double e[], int n, double **z);

Err tridag(double a[], double b[], double c[], double r[], double u[],
           unsigned long n);

#endif
