/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/trisparselinearsolver_fixedpoint.cpp
// Purpose:     Fixed point solver class for linear sparse system with tridiagonal part
// Created:     2005/01/16
// RCS-ID:      $Id: trisparselinearsolver_fixedpoint.cpp,v 1.3 2005/06/01 12:53:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// #define STEPDEBUG

#ifdef STEPDEBUG
#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"
#endif

#include <cmath>

#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"
#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/morsematrix.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TriSparseLinearSolverFixedPoint);
}

namespace ito33
{

namespace numeric
{


// Solve the linear sparse system using the fixed point method
void TriSparseLinearSolverFixedPoint::Solve(const TridiagonalMatrix& matrix, 
                                            const MorseMatrix& sparseMatrix,
                                            const double *pdRHS, 
                                            double* pdX)
{
  size_t nNbX = matrix.Dimension();

  size_t nIdx;

  SetMatrixDimension(nNbX);

  m_WorkSolver.EnableTransposeMatrixSolving(m_bIsTranposed);

  #ifdef STEPDEBUG
    std::cout << "  Fixed point Iteration: " << nNbX << " ";
  #endif

  size_t nIterCount;
  for (nIterCount = 0; nIterCount < m_nNbIterMax; nIterCount++)
  {  
    for (nIdx = 0; nIdx < nNbX; nIdx++)
       m_pdIterationPrices[nIdx] = pdX[nIdx];

    sparseMatrix.ProductMatrixVector(pdX, m_pdWorkRHS.Get());

    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = pdRHS[nIdx] - m_pdWorkRHS[nIdx];

    #ifdef STEPDEBUG
      std::cout << "*";
    #endif

    m_WorkSolver.Solve(matrix, m_pdWorkRHS.Get(), pdX);

    bool bConverged = true;
    double dError = 0.0;

    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {
      double scale = (pdX[nIdx] > 1.0) ? pdX[nIdx] : 1.0;
      double dTmpError = fabs(pdX[nIdx] - m_pdIterationPrices[nIdx]) / scale;
      dError = ( dTmpError > dError ) ? dTmpError : dError;

      if ( dError > 1.e-10 )
      {
        bConverged = false;
        break;
      }
    }

    if (bConverged)
      break;
  }

  #ifdef STEPDEBUG
    std::cout << std::endl;
  #endif

  if ( nIterCount == m_nNbIterMax ) // reached the maximum number of iterations
    throw EXCEPTION_MSG
    (
      ITO33_MAX_ITER,
      TRANS("Maximum number of iterations exceeded in fixed point solver.")
    );
}


} // namespace numeric

} // namespace ito33
