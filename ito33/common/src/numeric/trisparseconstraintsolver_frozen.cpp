/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/trisparseconstraintsolver_penalty.cpp
// Purpose:     penalty solver class for sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparseconstraintsolver_frozen.cpp,v 1.13 2006/01/26 15:39:09 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

//#define STEPDEBUG

#define INCOMPLETE_MATRIX_SOLVE

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
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparseconstraintsolver_frozen.h"

#include "ito33/pricing/constraints.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TriSparseConstraintSolverFrozen);
}

namespace ito33
{

namespace numeric
{


TriSparseConstraintSolverFrozen::TriSparseConstraintSolverFrozen
    (size_t nNbSubSystems, size_t nMaxUnknowns, size_t nNbIterMax) 
  :  TriSparseConstraintSolver(nNbSubSystems, nMaxUnknowns, nNbIterMax),
     m_WorkSolver(nMaxUnknowns)
{
}

// Solve the constrained system using the penalty method
void TriSparseConstraintSolverFrozen::Solve
     (const TridiagonalMatrix& matrix, const MorseMatrix& sparseMatrix,
      const double *pdRHS, 
      const pricing::Constraints* pConstraints,
      double* pdX,
      int* piFlags)
{
  size_t nNbX = matrix.Dimension();

#ifdef INCOMPLETE_MATRIX_SOLVE
#else
  size_t nNbS = nNbX / m_nNbSubSystems;
  
  ASSERT(nNbS * m_nNbSubSystems == nNbX);
#endif

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

  for (nIdx = 0; nIdx < nNbX; nIdx++)
  {
    pdBTmp[nIdx] = pdB[nIdx];
    m_pdIterationPrices[nIdx] = - 1.;
  }

  // Use the current flags to constrain the pdX entries, if the flag is set
  pConstraints->Get(pdX, piFlags, nNbX);

  // The algorithm has not been restarted
  bool bIsReset = false;

  #ifdef STEPDEBUG
    std::cout << "  Frozen Iteration: " << nNbX << " ";
  #endif

  size_t nIterCount;
  for (nIterCount = 0; nIterCount < m_nNbIterMax; nIterCount++)
  {
 
    #ifdef STEPDEBUG
      std::cout << "*";
    #endif    
    
#ifdef  INCOMPLETE_MATRIX_SOLVE 
    sparseMatrix.ProductMatrixVector(pdX, m_pdWorkRHS.Get());

    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = pdRHS[nIdx] - m_pdWorkRHS[nIdx]; 
#else
    m_WorkSparseMatrix.Init(sparseMatrix.GetMorseStruct(), 0);
    
    // Get the right hand side 
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = - pdRHS[nIdx];
        
    // A version that use the particular structure of the matrix
    for (size_t nIdxJ = 0; nIdxJ < nNbX; nIdxJ++)
    {
      if ( piFlags[nIdxJ] )
      {
        m_WorkSparseMatrix.SetColumn(nIdxJ);
        sparseMatrix.ProductMatrixVector(pdX, m_pdWorkRHS.Get(), nIdxJ);
      }
      else
        m_WorkSparseMatrix.SetColumn(sparseMatrix, nIdxJ);
    }
   
    
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = - m_pdWorkRHS[nIdx];  
#endif

    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {
      // Not the last unknown
      if ( nIdx != nNbX - 1 )
        if ( piFlags[nIdx + 1] )
        {
          m_pdWorkRHS[nIdx] -= pdC[nIdx] * pdX[nIdx + 1];
          pdCTmp[nIdx] = 0.;
        }
        else
          pdCTmp[nIdx] = pdC[nIdx];

      // Not the first unknown
      if ( nIdx != 0 ) 
        if ( piFlags[nIdx - 1] )
        {
          m_pdWorkRHS[nIdx] -= pdA[nIdx] * pdX[nIdx - 1];
          pdATmp[nIdx] = 0.; 
        }
        else
          pdATmp[nIdx] = pdA[nIdx];
    }
   
    // todo: a complete solving might be avoided to speed up the overall
    // performance, as with penalty solver, but a balance needs to be found
    // A possible choice is to iterate only a fixed number, or to require a 
    // lower precision
#ifdef INCOMPLETE_MATRIX_SOLVE
   m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdX);
#else
    m_pWorkSolver->Solve(m_WorkMatrix, m_WorkSparseMatrix, 
                         m_pdWorkRHS.Get(), pdX);
#endif

    // copy the current constraint flags for later comparison
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_piIterationFlags[nIdx] = piFlags[nIdx]; 

    // pdX[i] is used to set piFlag[i], and if constraint[i] is active, 
    // then pdX[i] is also reset to the constrained value
    pConstraints->Apply(pdX, piFlags, nNbX);
    
    // Compare the saved constraint values and new constraint values
#ifdef INCOMPLETE_MATRIX_SOLVE   
    // Check price error
    for ( nIdx = 0; nIdx < nNbX-1; nIdx++)
    {
      double dDenom = ( fabs(pdX[nIdx]) > 1.0 ) ? fabs(pdX[nIdx]) : 1.0;
      if ( fabs( m_pdIterationPrices[nIdx] - pdX[nIdx] )/dDenom > 1.e-12 )
        break;
    }

    if(piFlags[nNbX - 1] != 0 && piFlags[nNbX - 2] == 0)
      piFlags[nNbX - 1] = 0;
    if(piFlags[0] != 0 && piFlags[1] == 0)
      piFlags[0] = 0;

    // found the optimal solution
    if (nIdx == nNbX-1)
      break; 
#else
    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {
      if (   m_piIterationFlags[nIdx] != piFlags[nIdx]
          // sometimes we are in the this situation: at the two consecutive
          // iterations, the flags at a point are different but the price
          // values are the same (exactly the constraint value). We can 
          // consider this point as the free boundary and say that the solver
          // converges at this point.
          && m_pdIterationPrices[nIdx] != pdX[nIdx]  
         )
      {
        break;
      }
    }
 
    // The following situation may happen and cause non-convergence. So we     
    // have to correct it artificially for each sub system
    for (size_t nIdxSS = 0; nIdxSS < m_nNbSubSystems; nIdxSS++)
    {
      int* piFlagsTmp = piFlags + nIdxSS * nNbS;
      
      if ( piFlagsTmp[nNbS - 1] != 0 && piFlagsTmp[nNbS - 2] == 0)
        piFlagsTmp[nNbS - 1] = 0;
      if ( piFlagsTmp[0] != 0 && piFlagsTmp[1] == 0 )
        piFlagsTmp[0] = 0;
    }
    
    if (nIdx == nNbX)
      break; 
#endif

    #ifdef STEPDEBUG
      std::cout << nIdx << " ";
    #endif

    // save the price for comparison when a flag flips
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdIterationPrices[nIdx] = pdX[nIdx];  

    // While searching for the free boundary, the frozen method tends to
    // move by only one or two nodes each iteration.  If the actual location
    // has jumped between timesteps (eg. new constraint introduced), it is
    // better to reset the flags, make an explicit step, then let the code 
    // approximate where the new boundary is located. The usual algorithm
    // can then restart with a good initial guess. In case a fixed point
    // interation (incomplete solve) is being used, allow several iterations
    // to pass so that the fixed point aspect of the solve is done
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

  // reached the maximum number of iterations
  if ( nIterCount == m_nNbIterMax ) 
    throw EXCEPTION_MSG
    (
      ITO33_MAX_ITER,
      TRANS("Maximum number of iterations exceeded in frozen solver.")
    );
}


} // namespace numeric

} // namespace ito33
