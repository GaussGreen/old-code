/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/tridiagonalpenaltysolver.cpp
// Purpose:     penalty solver class for tridiagonal system with constraints
// Author:      ZHANG Yunzhi
// Created:     2004/02/10
// RCS-ID:      $Id: tridiagonalpenaltysolver.cpp,v 1.20 2006/06/12 20:00:33 dave Exp $
// Copyright:   (c) 2003- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file numeric/tridiagonalpenaltysolver.cpp
    @brief penalty solver class for tridiagonal system with constraints

    Implementation of penalty solver class for tridiagonal system with
    constraints.
 */

//#define PENALTYSTEPDEBUG
#ifdef PENALTYSTEPDEBUG
#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"
#endif

#include <cmath>

#include "ito33/useexception.h"
#include "ito33/autoptr.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/tridiagonalpenaltysolver.h"

#include "ito33/pricing/constraints.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalPenaltySolver);
}

namespace ito33
{

namespace numeric
{


// Solve the constrained system using the penalty method
void TridiagonalPenaltySolver::Solve(const TridiagonalMatrix &matrix, 
                                     const double *pdRHS, 
                                     const pricing::Constraints *pConstraints,
                                     double *pdX,
                                     int *piFlag)
{
  // Since the original matrix is const, either we can cast the const away,
  // save the diagonal entries before the iterations, and restore them at the
  // end, OR, copy the off-diagonal entries into a work matrix.  The 
  // complexity is the same.  Use the latter approach to be consistent with
  // the frozen method
  //TridiagonalMatrix& WorkMatrix = const_cast<TridiagonalMatrix &>(matrix);

  // Tolerance is roughly the accuracy to which the system will be solved
  const double dLarge = 1.e8;
  const double dTolerance = 1./dLarge;

  size_t nNbX = matrix.Dimension();

  size_t nIdx;

  SetMatrixDimension(nNbX);

  const double
    *pdA = matrix.GetA(),
    *pdB = matrix.GetB(),
    *pdC = matrix.GetC();

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();

  int* piIterationFlag = m_piIterationFlag.Get();

  // The penalty algorithm needs a copy of the constraints independent of 
  // the prices.
  Array<double> pdConstraintsArray(nNbX);
  double* pdConstraints = pdConstraintsArray.Get();

  // Save the position of the flag-flip from the last iteration
  size_t nIdxLastFlip = 0;

  // The algorithm has not been restarted
  bool bIsReset = false;

  // Copy the off-diagonals of the matrix into the work matrix
  // Also initialize the "old" iteration prices. 
  for (nIdx = 0; nIdx < nNbX; nIdx++)
  {
    pdATmp[nIdx] = pdA[nIdx];
    pdCTmp[nIdx] = pdC[nIdx];
    m_pdIterationPrice[nIdx] = 0.0;

    pdConstraints[nIdx] = 0.0;
  }

  // Use the current flags to get the active consraints
  pConstraints->Get(pdConstraints, piFlag, nNbX);

  # ifdef PENALTYSTEPDEBUG
    std::cout << "  Penalty Iteration: " << nNbX << " ";
  # endif

  size_t nIterCount;
  for(nIterCount = 0; 
      nIterCount < m_nNbIterMax; 
      nIterCount++)
  {

    // Apply the penalty method
    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {      
      if (piFlag[nIdx] ) // If the constraint is active ...
      {
        m_pdWorkRHS[nIdx] = pdRHS[nIdx] + dLarge * pdConstraints[nIdx];
        pdBTmp[nIdx] = pdB[nIdx] + dLarge;
      }
      else
      {
        m_pdWorkRHS[nIdx] = pdRHS[nIdx]; // copy rhs
        pdBTmp[nIdx] = pdB[nIdx];
      }
    }
 
    # ifdef PENALTYSTEPDEBUG
      std::cout << "*";
    # endif

    m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdX);

    // Save the old constraints and re-apply using new solution. On exit,
    // piFlag == piIterationFlag, so flipping ptrs is valid
    int* piTmp = piFlag;
    piFlag = piIterationFlag;
    piIterationFlag = piTmp;

    // pdX[i] is used to set piFlag[i], and if piFlag[i] is set, 
    // then pdConstraints[i] is set to the constrained value
    pConstraints->ApplyPenalty(pdX, piFlag, pdConstraints, nNbX);
    
    // Compare the saved constraint values and new constraint values.
    // For some reason (currently unknown and not investigated), the boundary
    // points sometimes flip for no apparent reason. Ignore them.
    for (nIdx = 1; nIdx < nNbX-1; nIdx++)
    {
      if (piIterationFlag[nIdx] != piFlag[nIdx])  
        break;
    }

    // found the optimal solution
    if (nIdx == nNbX-1)
      break;    

    // Check if it is simply numerical round-off that caused a flag to flip.
    // In theory, this check is not needed. In practice, round-off may cause
    // a flag to oscillate.  To avoid some work, only start checking where
    // the previous iteration flipped, or where the current iteration flipped
    // (whichever is less).  Typically, the code breaks immediately.
    size_t nIdx2 = nIdx;
    if (nIdx2 > nIdxLastFlip)
      nIdx2 = nIdxLastFlip;
    for ( ; nIdx2 < nNbX-1; nIdx2++)
    {
      double dDenom = ( fabs(pdX[nIdx2]) > 1.0 ) ? fabs(pdX[nIdx2]) : 1.0;
      if ( fabs( m_pdIterationPrice[nIdx2] - pdX[nIdx2] )/dDenom > dTolerance )
        break;
    }

    // found the optimal solution
    if (nIdx2 == nNbX-1)    
      break;        

    # ifdef PENALTYSTEPDEBUG
      std::cout << nIdx << " ";
    # endif

    // save flip position and the prices for comparison at next iteration
    nIdxLastFlip = nIdx;
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdIterationPrice[nIdx] = pdX[nIdx];
    
    // While searching for the free boundary, this algorithm tends to move 
    // by only one or two nodes each iteration.  If the actual location
    // has jumped between timesteps (eg. new constraint introduced), it is
    // better to reset the flags, make an explicit step, then let the code 
    // approximate where the new boundary is located. The usual algorithm
    // can then restart with a good initial guess
    if ( nIterCount > 6 && bIsReset == false )
    {
#     ifdef PENALTYSTEPDEBUG
      std::cout << "Restarting the iterations ";
#     endif

      bIsReset = true;
      for (nIdx = 0; nIdx < nNbX; nIdx++)
        piFlag[nIdx] = 0;
    }
    
  } // iteration loop

  # ifdef PENALTYSTEPDEBUG
    std::cout << std::endl;
  # endif

  // Clean-up the flags by applying the constraint with a tolerance in 
  // case the flags are oscillating near the boundary. This has no effect
  // on the price, but can affect Greek computations which use the flags.
  // (see bug 1117)
  pConstraints->Apply(pdX, piFlag, nNbX);

  if (nIterCount == m_nNbIterMax) // reached the maximum number of iterations
    throw EXCEPTION_MSG
    (
      ITO33_MAX_ITER,
      TRANS("Maximum number of iterations exceeded in penalty solver.")
    );

}

} // namespace pricing

} // namespace ito33
