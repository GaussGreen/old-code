/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/trisparselinearsolver.h
// Purpose:     Base solver class for linear sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparselinearsolver.h,v 1.4 2005/06/01 12:51:08 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/trisparselinearsolver.h
    @brief Base solver class for linear sparse system 
 */

#ifndef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_H_
#define _ITO33_NUMERIC_TRISPARSELINEARSOLVER_H_

#include "ito33/common.h"

namespace ito33
{

namespace numeric
{

  class TridiagonalMatrix;
  class MorseMatrix;

/// Fixed point solver class for linear sparse system with tridiagonal part
class TriSparseLinearSolver
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
   */
  TriSparseLinearSolver(size_t nMaxUnknowns, size_t nNbIterMax) 
                      : m_nMaxUnknowns(nMaxUnknowns),
                        m_nNbIterMax(nNbIterMax),
                        m_bIsTranposed(false)
  {
  }

  /// virtual dtor for base class
  virtual ~TriSparseLinearSolver() { }
 
  void EnableTransposeMatrixSolving(bool bIsTranposed = true)
  {
    m_bIsTranposed = bIsTranposed;
  }

  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param sparseMatrix a sparse matrix 
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pdX (output) the array of the solution
     
     the input constraint flag values serve as initial guess.
   */
  virtual void Solve(const TridiagonalMatrix& matrix,
                     const MorseMatrix& sparseMatrix,
                     const double* pdRHS, 
                     double* pdX) = 0;


protected:

  /// The maximum number of unknowns, used to allocate memory
  size_t m_nMaxUnknowns;
  
  /// The maximum number of iterations
  size_t m_nNbIterMax;

  /// Boolean indicates if we should soving M^T X = B instead of MX = B
  bool m_bIsTranposed;

}; // class TriSparseLinearSolver


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRISPARSELINEARSOLVER_H_
