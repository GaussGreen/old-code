/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalconstraintsolver.cpp
// Purpose:     base solver class for tridiagonal system with constraints
// Author:      ZHANG Yunzhi
// Created:     2004/02/10
// RCS-ID:      $Id: tridiagonalconstraintsolver.cpp,v 1.6 2005/12/27 16:11:47 wang Exp $
// Copyright:   (c) 2003- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalfrozensolver.cpp
    @brief base solver class for tridiagonal system with constraints

    impelementation of base solver class for tridiagonal system with
    constraints.
 */

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/tridiagonalconstraintsolver.h"

using namespace ito33;
using namespace ito33::numeric;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalConstraintSolver);
}

void TridiagonalConstraintSolver::SolveGreek(const TridiagonalMatrix &matrix,
                                             const double *pdRHS, 
                                             const int *piFlag,
                                             double *pdGreeks)
{
  size_t nNbX = matrix.Dimension();
  size_t nIdx;

  Resize( nNbX );

  m_WorkMatrix = matrix;

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();

  for(nIdx = 0; nIdx < nNbX; nIdx++)
  {
    if(piFlag[nIdx]) // in constraints case
    {
      pdATmp[nIdx] = 0;
      pdBTmp[nIdx] = 1;
      pdCTmp[nIdx] = 0;

      m_pdWorkRHS[nIdx] = 0;
    }
    else
      m_pdWorkRHS[nIdx] = pdRHS[nIdx];
  }

  // now solve the updated system
  m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdGreeks);
  
}

void TridiagonalConstraintSolver::SolveGreekWithSpecialConstraint
                                  ( const TridiagonalMatrix &matrix,
                                    const double *pdRHS, 
                                    const int *piFlag,
                                    const double* pdConstraintGreeks,
                                    double *pdGreeks)
{
  size_t nNbX = matrix.Dimension();
  size_t nIdx;

  Resize( nNbX );

  m_WorkMatrix = matrix;

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();

  for(nIdx = 0; nIdx < nNbX; nIdx++)
  {
    if(piFlag[nIdx]) // in constraints case
    {
      pdATmp[nIdx] = 0;
      pdBTmp[nIdx] = 1;
      pdCTmp[nIdx] = 0;

      m_pdWorkRHS[nIdx] = pdConstraintGreeks[nIdx];
    }
    else
      m_pdWorkRHS[nIdx] = pdRHS[nIdx];
  }

  // now solve the updated system
  m_WorkSolver.Solve(m_WorkMatrix, m_pdWorkRHS.Get(), pdGreeks);
  
}
