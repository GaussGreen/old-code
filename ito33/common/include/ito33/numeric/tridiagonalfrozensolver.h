/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalfrozensolver.h
// Purpose:     frozen solver class for tridiagonal system with constraints
// Author:      ZHANG Yunzhi
// Created:     2004/02/10
// RCS-ID:      $Id: tridiagonalfrozensolver.h,v 1.8 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalfrozensolver.h
    @brief frozen solver class for tridiagonal system with constraints
 */

#ifndef _ITO33_NUMERIC_TRIDIAGONALFROZENSOLVER_H_
#define _ITO33_NUMERIC_TRIDIAGONALFROZENSOLVER_H_

#include "ito33/numeric/tridiagonalconstraintsolver.h"
#include "ito33/common.h"

namespace ito33
{

namespace pricing
{
  class Constraints;
}

namespace numeric
{

/// frozen solver class for tridiagonal system with constraints
class TridiagonalFrozenSolver: public TridiagonalConstraintSolver
{

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
     @param nNbIterMax The max iterations number of solver
   */
  TridiagonalFrozenSolver(size_t nMaxUnknowns = 1, 
                          size_t nNbIterMax = 40) :
     TridiagonalConstraintSolver(nMaxUnknowns, nNbIterMax),
     m_piIterationFlag(nMaxUnknowns),
     m_pdIterationPrice(nMaxUnknowns)
  {
  }

  /// dtor
  virtual ~TridiagonalFrozenSolver() {}
 
  /** 
     Solve a tridiagonal system with constraint

     @param matrix a tridiagonal matrix
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pConstraints the constraints for the solution pdX
     @param pdX (output) the array of the solution
     @param piFlag (output) the contraints flag for the solution
   */
  void Solve(const TridiagonalMatrix &matrix, 
             const double *pdRHS, 
             const pricing::Constraints *pConstraints,
             double *pdX,
             int *piFlag);


protected:

  /**
    Let the solver work on a system whose dimension is nDimension.

    When nDimension is greater the the maximum size of the work arrays.
    This function does re-allocate memory for them.
    
    @param nDimension The dimension of the system to be solved.
   */
  void SetMatrixDimension(size_t nDimension)
  {
    m_WorkMatrix.SetDimension(nDimension);

    if (nDimension > m_nMaxUnknowns)
    {
      TridiagonalConstraintSolver::Resize(nDimension);

      m_piIterationFlag = Array<int>(nDimension);
      m_pdIterationPrice = Array<double>(nDimension);
    }
  }

  /// Flags helper for the frozen method
  Array<int> m_piIterationFlag;

  /// helper for the frozen method
  Array<double> m_pdIterationPrice;

private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TridiagonalFrozenSolver);

}; // class TridiagonalFrozenSolver


} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_TRIDIAGONALFROZENSOLVER_H_
