/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/trisparseconstraintsolver.h
// Purpose:     Solver class for non linear sparse system with tridiagonal part
// Created:     2005/01/15
// RCS-ID:      $Id: trisparseconstraintsolver.h,v 1.7 2005/06/01 14:25:10 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/trisparseconstraintsolver.h
    @brief Solver class for non linear sparse system with tridiagonal part
 */

#ifndef _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_H_
#define _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_H_

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/morsematrix.h"

namespace ito33
{

namespace pricing
{
  class Constraints;
}

namespace numeric
{

  class TriSparseLinearSolver;

/// Solver class for non linear sparse system with tridiagonal part
class TriSparseConstraintSolver 
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
     @param nNbSubsystems The number of sub system that the whole system
                          is composed of.
     REQUIRE: nNbIterMax must not be zero.
   */
  TriSparseConstraintSolver(size_t nNbSubSystems, 
                            size_t nMaxUnknowns, size_t nNbIterMax);

  /// virtual dtor for base class
  virtual ~TriSparseConstraintSolver() {}
 
  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param sparseMatrix a sparse matrix 
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pConstraints the constraints for the solution pdX
     @param pdX (output) the array of the solution
     @param piFlags (output) the contraints flag for the solution
   */
  virtual void Solve(const TridiagonalMatrix& matrix, 
                     const MorseMatrix& sparseMatrix,
                     const double* pdRHS, 
                     const pricing::Constraints* pConstraints,
                     double* pdX,
                     int* piFlags) = 0;


  /** 
     Solve a tridiagonal system for greeks, with known constraints flag.

     In fact, this function works for greeks type variables. That is, when
     the constraints has been applied to the price at a point, (the 
     constraints flag of the price is different from 0), this variable value
     at this point is systematically 0. "Fugit" is one of the exemples.

     @param matrix a tridiagonal matrix
     @param sparseMatrix a sparse matrix 
     @param pdRHS the right hand side array of the tirdiagonal system
     @param piFlags the contraints flag for the price
     @param pdGreeks (output) the array of the solution
   */
  virtual void SolveGreek(const TridiagonalMatrix& matrix, 
                          const MorseMatrix& sparseMatrix,
                          const double *pdRHS, 
                          const int* piFlags,
                          double* pdGreeks);


protected:

  /// Really change the size of the solver
  void SetDimension(size_t nSize)
  {
    m_WorkMatrix.SetDimension(nSize);
    
    if (m_nMaxUnknowns < nSize)
    {
      m_nMaxUnknowns = nSize; 

      // Resize for m_pWorkSolver and/or m_WorkSparseMatrix is done in 
      // Solve or SolveGreeks
      
      m_pdWorkRHS = Array<double>(m_nMaxUnknowns);

      m_piIterationFlags = Array<int>(m_nMaxUnknowns); 

      m_pdIterationPrices = Array<double>(m_nMaxUnknowns); 
    }
  }

  /// The number of sub system that the whole system is composed of.
  size_t m_nNbSubSystems;

  /// The size of the tridiagonal matrices to be solved
  size_t m_nMaxUnknowns;

  /// The tridiagonal work matrix for the solver
  TridiagonalMatrix m_WorkMatrix;

  /// The work array for right hand side
  Array<double> m_pdWorkRHS;

  /// The max number of iterations 
  const size_t m_nNbIterMax;

  // Intermediate flags
  Array<int> m_piIterationFlags;

  /// Intermediate iteration prices
  Array<double> m_pdIterationPrices;

  /// The Greek solver needs a normal linear solver, as do some derived classes
  AutoPtr<TriSparseLinearSolver> m_pWorkSolver;

  /// The sparse work matrix for the solver
  MorseMatrix m_WorkSparseMatrix;

private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TriSparseConstraintSolver);

}; // class TriSparseConstraintSolver 


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRISPARSECONSTRAINTSOLVER_H_
