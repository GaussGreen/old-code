/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/tridiagonalsolver.cpp
// Purpose:     Preconditioner using a tridiagonal matrix
// Author:      Wang
// Created:     2004/04/27
// RCS-ID:      $Id: tridiagonalpreconditioner.cpp,v 1.4 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalpreconditioner.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalPreconditioner);

}


using ito33::numeric::TridiagonalPreconditioner;

// Solve the system matrix * pdX = pdRHS where matrix is a triadiagonal matrix.
// Return values in pdX, which must be initialized by the caller.
void TridiagonalPreconditioner::Solve(const double *pdRHS, double *pdX) const
{
  m_solver.Solve(m_matrix, pdRHS, pdX);
}



CVecteurDouble TridiagonalPreconditioner::Solve
               (const CVecteurDouble &pdRHS) const
{
  CVecteurDouble pdX( pdRHS.dim() );

  m_solver.Solve(m_matrix, &pdRHS[0], &pdX[0]);

  return pdX;
}
