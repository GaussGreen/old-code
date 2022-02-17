/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/computeexpmatrix.cpp
// Purpose:     Implementation of computation of matrix exp
// Created:     2006/07/06
// RCS-ID:      $Id: computeexpmatrix.cpp,v 1.1 2006/07/06 16:01:12 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @todo We can implement as well Pade approximation, or write it in term of
          ODE solving.

    @todo Compute the whole matrix might not be the most efficient way. 
          Sometimes, we need only the matrix vector.
 */

#include <cmath>

#include "ito33/debug.h"

#include "ito33/numeric/densematrix.h"
#include "ito33/numeric/computeexpmatrix.h"

namespace ito33
{

namespace numeric
{

void ComputeExpMatrix(size_t nNb, double const* const* ppdA, double** ppdExpA)
{
  ASSERT_MSG(nNb < 100, "Your matrix is too big to be efficiently computed!");

  // Find n such that || A / 2^n || < 1. To be safe, use 0.125 instead of 1
  double dNorm = EuclideanNorm(nNb, ppdA);
  size_t n = 0;
  while (dNorm > 0.125)
  {
    dNorm *= 0.5;
    n++;
  }

  // A / 2^n
  DenseMatrix matrixM(nNb);
  double dScale = 1. / pow(2, n);
  for (size_t nI = 0; nI < nNb; nI++)
    for (size_t nJ = 0; nJ < nNb; nJ++)
      matrixM[nI][nJ] = ppdA[nI][nJ] * dScale;

  DenseMatrix matrixTmp(nNb);

  // the taylor series up to 4 degree
  
  // I + A
  for (size_t nI = 0; nI < nNb; nI++)
    for (size_t nJ = 0; nJ < nNb; nJ++)
      ppdExpA[nI][nJ] = matrixM[nI][nJ];

  for (size_t nI = 0; nI < nNb; nI++)
    ppdExpA[nI][nI] += 1;

  // 1 / 2! A^2
  for (size_t nI = 0; nI < nNb; nI++)
    for (size_t nJ = 0; nJ < nNb; nJ++)
    {
      matrixTmp[nI][nJ] = 0;
      for (size_t nK = 0; nK < nNb; nK++)
        matrixTmp[nI][nJ] += matrixM[nI][nK] * matrixM[nK][nJ];

      ppdExpA[nI][nJ] += 0.5 * matrixTmp[nI][nJ];
    }

   // 1 / 3! A^3 + 1 / 4! A^4 = 1 / 3! A^2 ( A + 1 / 4 A^2)
  for (size_t nI = 0; nI < nNb; nI++)
    for (size_t nJ = 0; nJ < nNb; nJ++)
      matrixM[nI][nJ] += 0.25 * matrixTmp[nI][nJ];

  for (size_t nI = 0; nI < nNb; nI++)
    for (size_t nJ = 0; nJ < nNb; nJ++)
      for (size_t nK = 0; nK < nNb; nK++)
        ppdExpA[nI][nJ] += 1. / 6. * matrixTmp[nI][nK] * matrixM[nK][nJ];

  // now exp(A) = exp(M)^{2^n}
  for (size_t m = 0; m < n; m++)
  {
    for (size_t nI = 0; nI < nNb; nI++)
      for (size_t nJ = 0; nJ < nNb; nJ++)
      {
        matrixTmp[nI][nJ] = 0;
        for (size_t nK = 0; nK < nNb; nK++)
          matrixTmp[nI][nJ] += ppdExpA[nI][nK] * ppdExpA[nK][nJ];
      }

    for (size_t nI = 0; nI < nNb; nI++)
      for (size_t nJ = 0; nJ < nNb; nJ++)
        ppdExpA[nI][nJ] = matrixTmp[nI][nJ];
  }
}

} // namespace numeric

} // namespace ito33
