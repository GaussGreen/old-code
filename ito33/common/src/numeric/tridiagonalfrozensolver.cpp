/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/tridiagonalfrozensolver.cpp
// Purpose:     frozen solver class for tridiagonal system with constraints
// Author:      ZHANG Yunzhi
// Created:     2004/02/10
// RCS-ID:      $Id: tridiagonalfrozensolver.cpp,v 1.18 2006/05/24 14:11:21 nabil Exp $
// Copyright:   (c) 2003- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file numeric/tridiagonalfrozensolver.cpp
    @brief frozen solver class for tridiagonal system with constraints

    implementation of frozen solver class for tridiagonal system with
    constraints.
 */

//#define FROZENSTEPDEBUG
#ifdef FROZENSTEPDEBUG
  #include "ito33/beforestd.h"
  #include <iostream>
  #include "ito33/afterstd.h"
#endif

#include <math.h>

#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"
#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalfrozensolver.h"

#include "ito33/pricing/constraints.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalFrozenSolver);
}

namespace ito33
{

namespace numeric
{


// Solve the constrained system using the frozen method
void TridiagonalFrozenSolver::Solve(const TridiagonalMatrix &matrix, 
                                    const double *pdRHS, 
                                    const pricing::Constraints *pConstraints,
                                    double *pdX,
                                    int *piFlag)
{
  size_t nNbX = matrix.Dimension();
  size_t nIdx;

  SetMatrixDimension( nNbX );

  const double
    *pdA = matrix.GetA(),
    *pdB = matrix.GetB(),
    *pdC = matrix.GetC();

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();
 
  const double dTolerance = 1.e-7;
  
  // The algorithm has not been restarted
  bool bIsReset = false;

  // Use the existing values in piFlag to set the values in pdX
  // ie. pdX[n] = constrained value iff piFlag[n] != 0
  pConstraints->Get(pdX, piFlag, nNbX);

  // the diagonal part of the tridiagonal matrix at each iteration doesn't change
  for (nIdx = 0; nIdx < nNbX; nIdx++)
  {
    pdBTmp[nIdx] = pdB[nIdx]; 
    m_pdIterationPrice[nIdx] = -1;
  }

# ifdef FROZENSTEPDEBUG
    std::cout << "  Frozen Iteration: " << " ";
# endif


  size_t nCounter;
  for(nCounter = 0;
      nCounter < m_nNbIterMax;
      nCounter++)
  {

#   ifdef FROZENSTEPDEBUG
      std::cout << "*";
#   endif

    // Get the right hand side 
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdWorkRHS[nIdx] = pdRHS[nIdx];
        
    // Loop over the rows of the matrix. 
    // First row is special case
    if (piFlag[1])
    {
      m_pdWorkRHS[0] -= pdC[0] * pdX[1];
      pdCTmp[0] = 0.0;
    }
    else
      pdCTmp[0] = pdC[0];

    // Internal rows
    for (nIdx = 1; nIdx < nNbX - 1; nIdx++)
    {
      if (piFlag[nIdx + 1])
      {
        m_pdWorkRHS[nIdx] -= pdC[nIdx] * pdX[nIdx + 1];
        pdCTmp[nIdx] = 0.;
      }
      else
        pdCTmp[nIdx] = pdC[nIdx];

      if (piFlag[nIdx - 1])
      {
        m_pdWorkRHS[nIdx] -= pdA[nIdx] * pdX[nIdx - 1];
        pdATmp[nIdx] = 0.; 
      }
      else
        pdATmp[nIdx] = pdA[nIdx];
    }

    // last row is special case
    if (piFlag[nNbX - 2])
    {
      m_pdWorkRHS[nNbX - 1] -= pdA[nNbX - 1] * pdX[nNbX - 2];
      pdATmp[nNbX - 1] = 0.0;
    }
    else
      pdATmp[nNbX - 1] = pdA[nNbX - 1];


    // now solve the updated system
    m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdX);

    // copy the current constraint flags for later comparison
    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {
      m_piIterationFlag[nIdx] = piFlag[nIdx]; 
    }

    // update the constraint flags, and prices, if needed
    pConstraints->Apply(pdX, piFlag, nNbX);

    // Compare the saved constraint values and new constraint values
    for (nIdx = 0; nIdx < nNbX; nIdx++)
    {
      if (m_piIterationFlag[nIdx] != piFlag[nIdx]
        // sometimes we are in the this situation: at the two consecutive
        // iterations, the flags at a point are different but the price
        // values are the same (exactly the constraint value). We can 
        // consider this point as the free boundary and say that the solver
        // converges at this point.
        && m_pdIterationPrice[nIdx] != pdX[nIdx]  
        )      
        break;
    }

    /*
      The following situation may happen and causes non-convergence. So we
      have to correct it artificially.
      */
    
    if(piFlag[nNbX - 1] != 0 && piFlag[nNbX - 2] == 0)
      piFlag[nNbX - 1] = 0;
    if(piFlag[0] != 0 && piFlag[1] == 0)
      piFlag[0] = 0;


    // found the optimal solution
    if (nIdx == nNbX)
      break;    

    // Check if it is simply numerical round-off that caused a flag to flip
    for ( ; nIdx < nNbX-1; nIdx++)
    {
      double dDenom = ( fabs(pdX[nIdx]) > 1.0 ) ? fabs(pdX[nIdx]) : 1.0;
      if ( fabs( m_pdIterationPrice[nIdx] - pdX[nIdx] )/dDenom > dTolerance )
        break;
    }

    // found the optimal solution
    if (nIdx == nNbX-1)
      break; 

#   ifdef FROZENSTEPDEBUG
    std::cout << nIdx << " ";
#   endif

    // save the price for possibly free boundary check explained above.
    for (nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdIterationPrice[nIdx] = pdX[nIdx];

    // While searching for the free boundary, this algorithm tends to move 
    // by only one or two nodes each iteration.  If the actual location
    // has jumped between timesteps (eg. new constraint introduced), it is
    // better to reset the flags, make an explicit step, then let the code 
    // approximate where the new boundary is located. The usual algorithm
    // can then restart with a good initial guess
    if ( nCounter > 6 && bIsReset == false )
    {
#     ifdef FROZENSTEPDEBUG
      std::cout << "Restarting the iterations ";
#     endif

      bIsReset = true;
      for (nIdx = 0; nIdx < nNbX; nIdx++)
        piFlag[nIdx] = 0;
    }

  }

# ifdef FROZENSTEPDEBUG
    std::cout << std::endl;
# endif
    
  ASSERT_MSG( !(piFlag[nNbX - 1] != 0 && piFlag[nNbX - 2] == 0)
                && !(piFlag[0] != 0 && piFlag[1] == 0),
              "We are in trouble: we got constraints only at boundary point "
              "but not at the points around them. We probably have some "
              "mesh problem here.");

  if (nCounter == m_nNbIterMax) // reached the maximum number of iterations
    throw EXCEPTION_MSG
    (
      ITO33_MAX_ITER,
      TRANS("Maximum number of iterations exceeded")
    );

}


} // namespace pricing

} // namespace ito33
