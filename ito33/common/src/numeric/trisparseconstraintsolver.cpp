/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/trisparseconstraintsolver.cpp
// Purpose:     penalty solver class for non linear sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparseconstraintsolver.cpp,v 1.8 2005/06/29 15:45:11 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparseconstraintsolver.h"

using namespace ito33;
using namespace ito33::numeric;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TriSparseConstraintSolver);
}

TriSparseConstraintSolver::TriSparseConstraintSolver
    (size_t nNbSubSystems, size_t nMaxUnknowns, size_t nNbIterMax)
  :  m_nNbSubSystems(nNbSubSystems),
     m_nMaxUnknowns(nMaxUnknowns),
     m_WorkMatrix(nMaxUnknowns),                  
     m_pdWorkRHS(nMaxUnknowns),
     m_piIterationFlags(nMaxUnknowns),
     m_pdIterationPrices(nMaxUnknowns),
     m_nNbIterMax(nNbIterMax),
     m_pWorkSolver( new TriSparseLinearSolverFixedPoint(nMaxUnknowns) )
{
  ASSERT_MSG(nNbIterMax,
    "The max iterations number of matrix solver must not be zero.");
}

void TriSparseConstraintSolver::SolveGreek
     (const TridiagonalMatrix& matrix,
      const MorseMatrix& sparseMatrix,
      const double* pdRHS,  
      const int* piFlags,
      double* pdGreeks)
{
  // Basic idea: If the constraint is active at node i, then the price of 
  // the contract is fixed at node i.  Any sensitivity (derivative with 
  // respect to model param), is then zero at node i.  So, if we are
  // solving A x = RHS, then for all constrained nodes i, we set
  //   A(i,i) = 1
  //   A(i,j) = 0  for all j != i
  //   RHS(i) = 0
  // For all unconstrained nodes, we leave A and RHS unchanged.
  // 
  // Note that this argument does NOT apply for delta, gamma and theta.

  size_t nNbX = matrix.Dimension();

  SetDimension( nNbX );

  // When a constraint is active, the tridiagonal matrix will be set to
  // unit diagonal, and the sparse matrix is cleared.  However, we do
  // not want to change the matrices passed in (in fact, they are const),
  // so use the work matrices instead
  m_WorkMatrix = matrix;
  m_WorkSparseMatrix.Init(sparseMatrix.GetMorseStruct(), 0);

  double
    *pdATmp = m_WorkMatrix.GetA(),
    *pdBTmp = m_WorkMatrix.GetB(),
    *pdCTmp = m_WorkMatrix.GetC();

  for (size_t nIdx = 0; nIdx < nNbX; nIdx++)
  {
    if(piFlags[nIdx]) // in constraints case
    {
      pdATmp[nIdx] = 0;
      pdBTmp[nIdx] = 1;
      pdCTmp[nIdx] = 0;

      m_pdWorkRHS[nIdx] = 0;

      // The sparse matrix was already initialized to zero
    }
    else
    {
      m_pdWorkRHS[nIdx] = pdRHS[nIdx];

      m_WorkSparseMatrix.SetRow(sparseMatrix, nIdx);

      // The tridiagonal matrix was already copied
    }
  }

  // now solve the updated system
  m_pWorkSolver->Solve(m_WorkMatrix, m_WorkSparseMatrix, m_pdWorkRHS.Get(), 
                       pdGreeks);
}
