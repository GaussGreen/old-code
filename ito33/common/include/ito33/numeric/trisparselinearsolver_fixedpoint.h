/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/trisparselinearsolver_fixedpoint.h
// Purpose:     Fixed point solver class for linear sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparselinearsolver_fixedpoint.h,v 1.3 2005/03/30 11:23:36 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/trisparselinearsolver_fixedpoint.h
    @brief Fixed point solver class for linear sparse system 
 */

#ifndef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_FIXEDPOINT_H_
#define _ITO33_NUMERIC_TRISPARSELINEARSOLVER_FIXEDPOINT_H_

#include "ito33/common.h"
#include "ito33/array.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

namespace ito33
{

namespace numeric
{

  class MorseMatrix;

/// Fixed point solver class for linear sparse system with tridiagonal part
class TriSparseLinearSolverFixedPoint : public TriSparseLinearSolver
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
   */
  TriSparseLinearSolverFixedPoint
      (size_t nMaxUnknowns = 0, size_t nNbIterMax = 40) 
     : TriSparseLinearSolver(nMaxUnknowns, nNbIterMax),
       m_WorkSolver(nMaxUnknowns),
       m_pdIterationPrices(nMaxUnknowns),
       m_pdWorkRHS(nMaxUnknowns)
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


protected:

  /**
    Let the solver work on a system whose dimension is nDimension.

    When nDimension is greater the the maximum size of the work arrays.
    This function does re-allocate memory for them.
    
    @param nDimension The dimension of the system to be solved.
   */
  void SetMatrixDimension(size_t nDimension)
  {
    if (nDimension > m_nMaxUnknowns)
    {
      m_nMaxUnknowns = nDimension;

      m_pdWorkRHS = Array<double>(nDimension);

      m_pdIterationPrices = Array<double>(nDimension);
    }
  }

  /// The normal tridiagonal matrix solver
  TridiagonalSolver m_WorkSolver;
  
  /// The work array for right hand side
  Array<double> m_pdWorkRHS;

  /// helper for iterative method
  Array<double> m_pdIterationPrices;


private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TriSparseLinearSolverFixedPoint);

}; // class TriSparseLinearSolverFixedPoint


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_FIXEDPOINT_H_
