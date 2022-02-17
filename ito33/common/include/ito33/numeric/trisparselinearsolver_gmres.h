/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/trisparselinearsolver_gmres.h
// Purpose:     GMRES solver class for linear sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparselinearsolver_gmres.h,v 1.2 2005/03/30 11:23:36 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/trisparselinearsolver_gmres.h
    @brief GMRES solver class for linear sparse system 
 */

#ifndef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_GMRES_H_
#define _ITO33_NUMERIC_TRISPARSELINEARSOLVER_GMRES_H_

#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ito33/numeric/gmressolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

namespace ito33
{

namespace numeric
{

  class MorseMatrix;
  class TridiagonalMatrix;
  class TridiagonalPreconditioner;

/// Fixed point solver class for linear sparse system with tridiagonal part
class TriSparseLinearSolverGmres : public TriSparseLinearSolver
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
   */
  TriSparseLinearSolverGmres
      (size_t nMaxUnknowns = 10, size_t nNbIterMax = 40) 
     : TriSparseLinearSolver(nMaxUnknowns, nNbIterMax),
       m_WorkSolver(nMaxUnknowns, 35)
  {
  }
 
  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param sparseMatrix a sparse matrix 
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pdX (output) the array of the solution
     
     the input constraint flag values serve as initial guess.
   */
  void Solve(const TridiagonalMatrix& matrix,
             const MorseMatrix& sparseMatrix,
             const double* pdRHS, 
             double* pdX);

  CVecteurDouble operator*(CVecteurDouble &X) const;


protected:

  /// The tridiagonal matrix
  const TridiagonalMatrix* m_pTridiagonalMatrix;

  /// The sparse matrix
  const MorseMatrix* m_pSparseMatrix;

  /// The GMRES solver
  GmresSolver m_WorkSolver;

  /// The tridiagonal preconditioner
  AutoPtr<TridiagonalPreconditioner> m_pPreconditioner;


private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TriSparseLinearSolverGmres);

}; // class TriSparseLinearSolverGmres


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_GMRES_H_
