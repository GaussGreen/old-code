
#ifndef __INTALGEBRA_H__
#define __INTALGEBRA_H__

#include "INTSTL.h"
#include "INTUtilities.h"


// Cholesky decomposition of a symmetric positive definite square matrix M = L*L^t
 STLDoubleVectorVector& choldc(STLDoubleVectorVector& L, const STLDoubleVectorVector& M);

// LU decomposition of a non-singular square matrix A
void ludcmp(STLDoubleVectorVector& A, int n, STLIntegerVector& indx, double *d);



#endif
