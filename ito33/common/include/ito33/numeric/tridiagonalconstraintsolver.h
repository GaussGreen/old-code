/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalconstraintsolver.h
// Purpose:     solver class for tridiagonal system with constraints
// Author:      ICARE
// Created:     2004/02/10
// RCS-ID:      $Id: tridiagonalconstraintsolver.h,v 1.7 2005/12/27 16:12:07 wang Exp $
// Copyright:   (c) 2003- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalconstraintsolver.h
    @brief solver class for tridiagonal system with constraints
 */

#ifndef _ITO33_NUMERIC_TRIDIAGONALCONSTRAINTSOLVER_H_
#define _ITO33_NUMERIC_TRIDIAGONALCONSTRAINTSOLVER_H_

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/common.h"

namespace ito33
{

namespace pricing
{
  class Constraints;
}

namespace numeric
{

/** 
   solver class for tridiagonal system with constraints
 */
class TridiagonalConstraintSolver 
{

public: 
  
  /**
    Ctor creates a work array. 

    @param nMaxUnknowns The suggested maximum number for the unknowns
    @param nNbIterMax The max iterations number of solver

    REQUIRE: nNbIterMax must not be zero.
   */
  TridiagonalConstraintSolver(size_t nMaxUnknowns, 
                              size_t nNbIterMax) :
      m_nMaxUnknowns(nMaxUnknowns),
      m_WorkSolver(nMaxUnknowns),
      m_WorkMatrix(nMaxUnknowns + 1), // the ctor doesn't accept 0. Any way
                                      // the real dimension will be reset 
      m_pdWorkRHS(nMaxUnknowns),
      m_nNbIterMax(nNbIterMax)
  {
    ASSERT_MSG(nNbIterMax,
      "The max iterations number of matrix solver must not be zero.");
  }

  /// dtor
  virtual ~TridiagonalConstraintSolver() {}
 
  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pConstraints the constraints for the solution pdX
     @param pdX (output) the array of the solution
     @param piFlag (output) the contraints flag for the solution
   */
  virtual void Solve(const TridiagonalMatrix &matrix, 
                     const double *pdRHS, 
                     const pricing::Constraints *pConstraints,
                     double *pdX,
                     int *piFlag) = 0;
  
  /** 
     Solve a tridiagonal system for greeks, with known constraints flag.

     In fact, this function works for greeks type variables. That is, when
     the constraints has been applied to the price at a point, (the 
     constraints flag of the price is different from 0), this variable value
     at this point is systematically 0. "Fugit" is one of the exemples.

     @param matrix a tridiagonal matrix
     @param pdRHS the right hand side array of the tirdiagonal system
     @param piFlag the contraints flag for the price
     @param pdGreeks (output) the array of the solution
   */
  virtual void SolveGreek(const TridiagonalMatrix &matrix, 
                          const double *pdRHS, 
                          const int *piFlag,
                          double *pdGreeks);

  /** 
     Solve a tridiagonal system for greeks, with known constraints flag,
     but with a special constraint.

     In fact, this function works for greeks type variables. That is, when
     the constraints has been applied to the price at a point, (the 
     constraints flag of the price is different from 0), this variable value
     at this point depends on the greek value of the constraint at this point. 
     "Vega for CBOption" is one of the exemples.

     @param matrix a tridiagonal matrix
     @param pdRHS the right hand side array of the tirdiagonal system
     @param piFlag the contraints flag for the price
     @param pdConstraintGreeks the array of the greeks of the constraint
     @param pdGreeks (output) the array of the solution

     @todo{The implementation of this function is just so similar to the 
           original function that it may be a good idea to find a way to share 
           common part of the implementation. A possible solution is to let the 
           SolveGreek function take a callback function that will return the 
           greek value at constrainted point.}
   */
  virtual void 
  SolveGreekWithSpecialConstraint(const TridiagonalMatrix &matrix, 
                                   const double *pdRHS, 
                                   const int *piFlag,
                                   const double *pdConstraintGreeks,
                                   double *pdGreeks); 

protected:

  /**
    Really change the size of the solver
    */
  void Resize(size_t nSize)
  {
    m_nMaxUnknowns = nSize;

    m_pdWorkRHS = Array<double>(nSize);

    // we don't need to resize m_WorkSolver here. In fact, it is done 
    // automaticly Solver::solve() function
  }

  /// The size of the tridiagonal matrices to be solved
  size_t m_nMaxUnknowns;

  /// The normal tridiagonal matrix solver
  TridiagonalSolver m_WorkSolver;

  /// The work matrix for the solver
  TridiagonalMatrix m_WorkMatrix;

  /// The work array for right hand side
  Array<double> m_pdWorkRHS;

  /// The max number of iterations 
  const size_t m_nNbIterMax;

private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TridiagonalConstraintSolver);

}; // class TridiagonalConstraintSolver 


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRIDIAGONALCONSTRAINTSOLVER_H_
