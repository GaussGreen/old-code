/* ========================================================================

  MODULE:	num_h_householder.h

  LIBRARY:	NUM_LIB

  FUNCTION:	householder_tridiagonalisation

  INPUTS:


  DESCRIPTION:	From Numerical Recipes in C page 474
                This function performs the Householder reduction of a real  ,
                                symmetric matrix a[0..n-1][0..n-1]. On output  ,
the matrix a is replaced by the orthogonal matrix effecting the transformation
                                d[0..n-1] returns the diagonal elements of the
tridiagonal matrix  , and e[0..n-1] the off diagonal elements  , with e[0] = 0

========================================================================== */

#ifndef NUM_H_HOUSEHOLDER_H
#define NUM_H_HOUSEHOLDER_H

Err householder_tridiagonalisation(double **a, int n, double d[], double e[]);

#endif