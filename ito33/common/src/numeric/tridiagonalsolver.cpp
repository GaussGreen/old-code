/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/tridiagonalsolver.cpp
// Purpose:     tridiagonal matrix solver class
// Author:      ICARE
// Created:     2003/11/20
// RCS-ID:      $Id: tridiagonalsolver.cpp,v 1.13 2005/06/01 12:53:41 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/exception.h"
#include "ito33/numeric/exception.h"
#include "ito33/error.h"
#include "ito33/gettext.h"
#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/numeric/tridiagonalsolver.h"

extern const ito33::Error ITO33_DIV0;

using namespace ito33::numeric;

using ito33::numeric::Exception;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalSolver);
}

// Solve the system matrix * pdX = pdRHS where matrix is a triadiagonal matrix.
// Return values in pdX, which must be initialized by the caller.
void TridiagonalSolver::Solve
     (const TridiagonalMatrix &matrix, const double *pdRHS, double *pdX) const
{
  // the real number of unknowns of the tridiagonal system
  size_t nNbUnknowns = matrix.Dimension();

  // make sure the size of the work array is big enough 
  if ( m_nMaxUnknowns < nNbUnknowns )
  {
    m_nMaxUnknowns = nNbUnknowns;
    m_pdWorkArray = Array<double>(m_nMaxUnknowns);
  }

  size_t nIdx;
  double dTmp;

  // pdA is the lower diagonal, pdB the main diagonal, pdC the upper diagonal 
  const double *pdB = matrix.GetB();
  const double *pdA, *pdC;
  
  if ( m_bIsTranposed )
  {
    pdA = matrix.GetC() - 1; 
    pdC = matrix.GetA() + 1;
  }
  else
  {
    pdA = matrix.GetA();
    pdC = matrix.GetC();
  }

  dTmp = pdB[0];

  // check a divide by zero on the first step of the algorithm
  if (dTmp == 0.0)
    throw EXCEPTION_MSG
          (
            ITO33_DIV0,
            TRANS("Division by zero in TridiagonalSolver")
          );

  pdX[0] = pdRHS[0] / dTmp;

  // decomposition and forward solve
  for (nIdx = 1; nIdx < nNbUnknowns; nIdx++) 
  {
    m_pdWorkArray[nIdx] = pdC[nIdx - 1] / dTmp;
    dTmp = pdB[nIdx] - pdA[nIdx] * m_pdWorkArray[nIdx];
    
    if (dTmp == 0.0)
      throw EXCEPTION_MSG
            (
              ITO33_DIV0,
              TRANS("Division by zero in TridiagonalSolver")
            ); 
    
    pdX[nIdx] = (pdRHS[nIdx]- pdA[nIdx] * pdX[nIdx - 1]) / dTmp;
  }

  // Backsolve
  for (nIdx = nNbUnknowns - 2; nIdx != size_t(-1); nIdx--)
    pdX[nIdx] -= m_pdWorkArray[nIdx + 1] * pdX[nIdx + 1];
}
