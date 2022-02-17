
#include	"stdafx.h"
#include "MCInt\INTAlgebra.h"



//----------------------------------------------------------------------------------------------

static STLDoubleVector ___useless_tmp___;

// Cholesky decomposition of a symmetric positive definite square matrix M = L*L^t,
// From NR in C, p. 97
STLDoubleVectorVector& choldc(STLDoubleVectorVector& L, const STLDoubleVectorVector& M)
{
  int n = M.size();
  int i, j, k;
  double sum;

  L.clear();
  L.resize(n, STLDoubleVector(n, 0.0));

  // Check if matrix is nxn square matrix
  for (i = 0; i < n; i++)
    {
    if (n != M[i].size()) ERROR("In choldc: M is not a square matrix. ")
    }

  // Check if M symmetric
  for (i = 0; i < n; i++)
    {
      for (j = i+1; j < n; j++)
        {
          if (M[i][j] != M[j][i]) ERROR("In choldc: M is not symmetric. ")
        }
    }

  for (i = 0; i < n; i++)
    {
      for (j = i; j < n; j++)
        {
          sum = M[i][j];
          for (k = 0; k <= i-1; k++)
            {
              sum -= L[i][k]*L[j][k];
            }
          if (i == j)
            {
              // Matrix not positive definite
              if (sum <= 0.0) ERROR("In choldc: M is not positive definite. ")
              L[i][i] = sqrt(sum);
            }
          else
            {
              L[j][i] = sum/L[i][i];
              L[i][j] = 0.0;
            }
        }
    }

  return L;
}

// A small number
#ifdef TINY
#undef TINY
#endif
#define TINY 1.0e-20 

// LU decomposition of a non-singular square matrix
// From NR in C, pp. 46 + 47
void ludcmp(STLDoubleVectorVector& A, int n, STLIntegerVector& indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;

  // vv stores the implicit scaling of each row.
  STLDoubleVector vv(n);

  // No row interchanges yet.
  *d = 1.0;

  // Loop over rows to get the implicit scaling information
	for (i = 0; i < n; i++)
    { 
      big = 0.0;

	    for (j = 0;j < n; j++)
        {
	        if ((temp = fabs(A[i][j])) > big) big = temp;
        }

      // Matrix singular
      if (big == 0) ERROR("In ludcmp: A is singular. ")

	    // Save the scaling
      vv[i] = 1.0/big; 
    }

  // This is the loop over columns of Crout's method
  for (j = 0; j < n; j++)
    { 
      // This is equation (2.3.12) except for i = j
		  for (i = 0; i < j; i++) 
        {
          sum = A[i][j];
			    for (k = 0; k < i; k++) sum -= A[i][k]*A[k][j];
			    A[i][j] = sum;
			  }

      // Initialize for the search for largest pivot element
      big = 0.0; 
			
      // This is i = j of equation (2.3.12) and i = j+1..N of equation (2.3.13).
      for (i = j; i < n; i++) 
	  { 
		  sum = A[i][j];
		  
		  for (k = 0; k < j; k++) sum -= A[i][k]*A[k][j];
          A[i][j] = sum;

          // Is the figure of merit for the pivot better than the best so far?
		  if ((dum = vv[i]*fabs(sum)) >= big) 
		  {
              big = dum;
			  imax = i;
		  }
	  }

      // Do we need to interchange rows?
      if (j != imax) 
	  {     
		  // Yes,do so...
		  for (k = 0; k < n; k++) 
          { 
			  dum = A[imax][k];
			  A[imax][k] = A[j][k];
			  A[j][k] = dum;
		  }

          // ...and change the parity of d
          *d = -(*d);

          // Also interchange the scale factor
		  vv[imax] = vv[j]; 
	  }

      indx[j] = imax;

      // If the pivot element is zero the matrix is singular (at least to the precision of the
      // algorithm). For some applications on singular matrices, it is desirable to substitute
      // TINY for zero.
      if (A[j][j] == 0.0) A[j][j] = TINY;

      // Now, finally, divide by the pivot element
	  if (j != n) 
	  { 
		  dum = 1.0/(A[j][j]);
		  for (i = j + 1; i < n; i++) A[i][j] *= dum;
	  }

      // Go back for the next column in the reduction
    } 
}


//----------------------------------------------------------------------------------------------

