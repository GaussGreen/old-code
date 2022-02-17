/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/trisparselinearsolver_gmres.cpp
// Purpose:     Fixed point solver class for linear sparse system with tridiagonal part
// Created:     2005/01/16
// RCS-ID:      $Id: trisparselinearsolver_gmres.cpp,v 1.4 2005/06/01 13:20:08 wang Exp $
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
#include "ito33/autoptr.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/morsematrix.h"
#include "ito33/numeric/tridiagonalpreconditioner.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TriSparseLinearSolverGmres);
}

namespace ito33
{

namespace numeric
{


// todo: optmize this function
CVecteurDouble TriSparseLinearSolverGmres::operator*(CVecteurDouble &X) const
{
  CVecteurDouble Prod(X.size()); 

  if ( m_bIsTranposed )
  {
    m_pSparseMatrix->ProductTransposeMatrixVector( X.get(), Prod.get() );
    m_pTridiagonalMatrix->ProductTransposeMatrixVector
                          (X.get(), Prod.get(), true);
  }
  else
  {
    m_pSparseMatrix->ProductMatrixVector( X.get(), Prod.get() );
    m_pTridiagonalMatrix->ProductMatrixVector(X.get(), Prod.get(), true);
  }

  return Prod;
}

// Solve the linear sparse system using the GMRES method
void TriSparseLinearSolverGmres::Solve(const TridiagonalMatrix& matrix, 
                                       const MorseMatrix& sparseMatrix,
                                       const double *pdRHS, 
                                       double* pdX)
{
  m_pTridiagonalMatrix = &matrix;

  m_pSparseMatrix = &sparseMatrix;

  m_pPreconditioner = AutoPtr<TridiagonalPreconditioner>
                      ( new TridiagonalPreconditioner(matrix) );

  m_pPreconditioner->EnableTransposeMatrixSolving(m_bIsTranposed);

  m_WorkSolver.Solve(*this, *m_pPreconditioner, pdRHS, pdX);

  #ifdef STEPDEBUG

  for (size_t nIdx = 0; nIdx < m_WorkSolver.GetNbIters(); nIdx++)
    std::cout << '*';

  std::cout << std::endl;
  
  #endif
}


} // namespace numeric

} // namespace ito33
