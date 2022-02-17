/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/trisparseconstraintsolver_penalty.cpp
// Purpose:     penalty solver class for sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparseconstraintsolver_penalty.cpp,v 1.11 2006/06/12 20:00:33 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

//#define STEPDEBUG

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
#include "ito33/numeric/trisparseconstraintsolver_penalty.h"

#include "ito33/pricing/constraints.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TriSparseConstraintSolverPenalty);
}

namespace ito33
{

namespace numeric
{


// Solve the constrained system using the penalty method
void TriSparseConstraintSolverPenalty::Solve
     (const TridiagonalMatrix& matrix, const MorseMatrix& sparseMatrix,
      const double *pdRHS, 
      const pricing::Constraints* pConstraints,
      double* pdX,
      int* piFlags)
{
  const double dLarge = 1.e9;
  const double dTolerance = 1./dLarge * 10.0;

  size_t nNbX = matrix.Dimension();

  size_t nIdx;

  SetDimension(nNbX);

  const double
    *pdA = matrix.GetA(),
    *pdB = matrix.GetB(),
    *pdC = matrix.GetC();

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();

  // Copy the off-diagonals of the matrix into the work matrix
  // Also initialize the "old" iteration prices. 
  for (nIdx = 0; nIdx < nNbX; nIdx++)
  {
    pdATmp[nIdx] = pdA[nIdx];
    pdCTmp[nIdx] = pdC[nIdx];
    m_pdIterationPrices[nIdx] = 0.0;
  }

  // The penalty algorithm needs a copy of the constraints independent of 
  // the prices.
  Array<double> pdConstraintsArray(nNbX);
  double* pdConstraints = pdConstraintsArray.Get();

  // Use the current flags to get the active constraints
  pConstraints->Get(pdConstraints, piFlags, nNbX);

  // The algorithm has not been restarted
  bool bIsReset = false;

  #ifdef STEPDEBUG
    std::cout << "  Penalty Iteration: " << nNbX << " ";
  #endif

  size_t nIterCount;
  for (nIterCount = 0; nIterCount < m_nNbIterMax; nIterCount++)
  {
    // Start of fixed point iteration (aka incomplete) matrix solve
    sparseMatrix.ProductMatrixVector(pdX, m_pdWorkRHS.Get());

    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = pdRHS[nIdx] - m_pdWorkRHS[nIdx];

    // Apply the penalty method
    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {    
      if (piFlags[nIdx] ) // If the constraint is active ...
      {
        m_pdWorkRHS[nIdx] += dLarge * pdConstraints[nIdx];
        pdBTmp[nIdx] = pdB[nIdx] + dLarge;
      }
      else
      {
        pdBTmp[nIdx] = pdB[nIdx];
      }
    }
   
    #ifdef STEPDEBUG
      std::cout << "*";
    #endif

    m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdX);

    // pdX[i] is used to set piFlag[i], and if piFlag[i] is set, 
    // then pdConstraints[i] is set to the constrained value
    pConstraints->ApplyPenalty(pdX, piFlags, pdConstraints, nNbX);
    
    // Check if the prices have converged.  Cannot simply look at the 
    // iteration flags since the fixed point aspect of the solve
    // may converge after the constraint flags
    for ( nIdx = 0; nIdx < nNbX - 1; nIdx++)
    {
      double dDenom = ( fabs(pdX[nIdx]) > 1.0 ) ? fabs(pdX[nIdx]) : 1.0;
      
      if ( fabs( m_pdIterationPrices[nIdx] - pdX[nIdx] ) / dDenom > dTolerance)
        break;
    }

    // found the optimal solution
    if (nIdx == nNbX - 1)
      break;    

    #ifdef STEPDEBUG
      std::cout << nIdx << " ";
    #endif

    // save the price for comparison when a flag flips
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdIterationPrices[nIdx] = pdX[nIdx];

    // While searching for the free boundary, the penalty method tends to
    // move by only one or two nodes each iteration.  If the actual location
    // has jumped between timesteps (eg. new constraint introduced), it is
    // better to reset the flags, make an explicit step, then let the code 
    // approximate where the new boundary is located. The usual algorithm
    // can then restart with a good initial guess. Allow several iterations
    // to complete so that the fixed point iteration aspect of the solve
    // is hopefully done.
    if ( nIterCount > 15 && bIsReset == false )
    {
#     ifdef STEPDEBUG
      std::cout << "Restarting the iterations ";
#     endif

      bIsReset = true;
      for (nIdx = 0; nIdx < nNbX; nIdx++)
        piFlags[nIdx] = 0;
    }
  }

  #ifdef STEPDEBUG
    std::cout << std::endl;
  #endif

  // Clean-up the flags by applying the constraint with a tolerance in 
  // case the flags are oscillating near the boundary. This has no effect
  // on the price, but can affect Greek computations which use the flags.
  // (see bug 1117)
  pConstraints->Apply(pdX, piFlags, nNbX);

  if (nIterCount == m_nNbIterMax) // reached the maximum number of iterations
    throw EXCEPTION_MSG
    (
      ITO33_MAX_ITER,
      TRANS("Maximum number of iterations exceeded in penalty solver.")
    );

  /*
  // This improves the calibration/accuracy a lot, but probably can be avoided
  // by playing with the flags and prices
  sparseMatrix.ProductMatrixVector(pdX, m_pdWorkRHS.Get());

  for (nIdx = 0; nIdx < nNbX; nIdx++)
    m_pdWorkRHS[nIdx] = pdRHS[nIdx] - m_pdWorkRHS[nIdx];

  // Apply the penalty method
  for (nIdx = 0; nIdx < nNbX; nIdx++)
  {    
    if (piFlags[nIdx] ) // If the constraint is active ...
    {
      m_pdWorkRHS[nIdx] += dLarge * pdX[nIdx];
      pdBTmp[nIdx] = pdB[nIdx] + dLarge;
    }
    else
    {
      pdBTmp[nIdx] = pdB[nIdx];
    }
  }
  
  m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdX);
  */
}


} // namespace pricing

} // namespace ito33
