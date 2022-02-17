/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalpreconditioner.h
// Purpose:     A preconditioner using a tridiagonal matrix
// Author:      Wang
// Created:     2004/04/27
// RCS-ID:      $Id: tridiagonalpreconditioner.h,v 1.3 2005/06/01 12:51:08 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalpreconditioner.h
    @brief A preconditioner using a tridiagonal matrix
 */

#ifndef _ITO33_NUMERIC_TRIDIAGONALPRECONDITIONER_H_
#define _ITO33_NUMERIC_TRIDIAGONALPRECONDITIONER_H_

#include "ito33/common.h"
#include "ito33/array.h"

#include "ito33/mv/mvvd.h"

#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"

namespace ito33
{

namespace numeric
{


/** 
   A preconditioner using a tridiagonal matrix

   It uses a tridiagonal solver to implement the method Solve.
 */
class TridiagonalPreconditioner
{ 

public: 
  
  /**
     Ctor stores a reference to a tridiagonal matrix,
     and creates a tridiagonal solver using the matrix.

     @param matrix the tridiagonal matrix to be used as preconditioner
   */
  TridiagonalPreconditioner(const TridiagonalMatrix & matrix) 
    : m_matrix(matrix), m_solver(matrix.Dimension())
  {
  }

  // default destructor is ok
 
  void EnableTransposeMatrixSolving(bool bIsTranposed = true)
  {
    m_solver.EnableTransposeMatrixSolving(bIsTranposed);
  }

  /** 
     Solve a tridiagonal system.

     @param pdRHS the right hand side array of the tirdiagonal system
     @param pdX the array of the solution
   */ 
  void Solve(const double *pdRHS, double *pdX) const;

  /**
     Solve method required by the caller of the preconitioner
     as for example GMRES

     @param pdRHS the right hand side array of the system
     @return the array of the solution

     @todo The actual implementation has big performance issue
           as the resulted solution is copied several times. A possible
           solution is to change the prototype required by gmres:
           the solution vector can be passed by reference as parameter.
           Note that in this case, attention to memory allocation!
   */
  CVecteurDouble Solve(const CVecteurDouble &pdRHS) const;


private:

  /// The tridiagonal preconitioner matrix 
  const TridiagonalMatrix &m_matrix;

  /// Work array used by the solving algorithm
  TridiagonalSolver m_solver;
 
  // Cannot copy or assign
  NO_COPY_CLASS(TridiagonalPreconditioner);

}; // class TridiagonalPreconditioner


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_TRIDIAGONALPRECONDITIONER_H_

