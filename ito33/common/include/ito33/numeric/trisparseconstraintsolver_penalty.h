/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/trisparseconstraintsolver_penalty.h
// Purpose:     penalty solver class for non linear sparse system 
// Created:     2005/01/15
// RCS-ID:      $Id: trisparseconstraintsolver_penalty.h,v 1.4 2005/05/19 09:21:40 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/trisparseconstraintsolver_penalty.h
    @brief penalty solver class for non linear psarse system 
 */

#ifndef _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_PENALTY_H_
#define _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_PENALTY_H_

#include "ito33/numeric/trisparseconstraintsolver.h"
#include "ito33/common.h"

namespace ito33
{

namespace pricing
{
  class Constraints;
}

namespace numeric
{

/// penalty solver class for tridiagonal system with constraints
class TriSparseConstraintSolverPenalty: public TriSparseConstraintSolver
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nNbSubsystems The number of sub system that the whole system
                          is composed of.
     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
   */
  TriSparseConstraintSolverPenalty
     (size_t nNbSubSystems, size_t nMaxUnknowns = 10, size_t nNbIterMax = 40) 
    : TriSparseConstraintSolver(nNbSubSystems, nMaxUnknowns, nNbIterMax),
      m_WorkSolver(nMaxUnknowns)
  {
  }
 
  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param sparseMatrix a sparse matrix 
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pConstraints the constraints for the solution pdX
     @param pdX (output) the array of the solution
     @param piFlags (input & output) the contraints flag for the solution
     
     the input constraint flag values serve as initial guess.
   */
  void Solve(const TridiagonalMatrix& matrix,
             const MorseMatrix& sparseMatrix,
             const double* pdRHS, 
             const pricing::Constraints* pConstraints,
             double* pdX,
             int* piFlags);


private:
 
  /// The normal tridiagonal matrix solver
  TridiagonalSolver m_WorkSolver;

  // Cannot copy or assign
  NO_COPY_CLASS(TriSparseConstraintSolverPenalty);

}; // class TriSparseConstraintSolverPenalty


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_PENALTY_H_
