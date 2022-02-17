/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalsolver.h
// Purpose:     tridiagonal matrix solver class
// Author:      ICARE
// Created:     2003/11/20
// RCS-ID:      $Id: tridiagonalsolver.h,v 1.12 2005/06/01 12:51:08 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalsolver.h
    @brief Tridiagonal matrix solver class

    Tridiagonal matrix solver class
 */

#ifndef _ITO33_NUMERIC_TRIDIAGONALSOLVER_H_
#define _ITO33_NUMERIC_TRIDIAGONALSOLVER_H_

#include "ito33/common.h"
#include "ito33/array.h"

#include "ito33/numeric/tridiagonalmatrix.h"


namespace ito33
{

namespace numeric
{


/** 
   A tridiagonal solver. 

   It stores the work array needed by the solver to avoid reallocation.
 */
class TridiagonalSolver
{ 

public: 
  
  /**
     Ctor creates a work array. 

     @param nMaxUnknowns The suggested maximum number for the unknowns
   */
  TridiagonalSolver(size_t nMaxUnknowns = 10)
                  : m_nMaxUnknowns(nMaxUnknowns),
                    m_bIsTranposed(false)
  {
    // Allocate here, if the size of the problem is bigger, should reallocate
    m_pdWorkArray = Array<double>(m_nMaxUnknowns);
  }

  // default destructor is ok
 
  void EnableTransposeMatrixSolving(bool bIsTranposed = true)
  {
    m_bIsTranposed = bIsTranposed;
  }

  /** 
     Solve a tridiagonal system.

     @param matrix a tridiagonal matrix
     @param pdRHS the right hand side array of the tirdiagonal system
     @param pdX the array of the solution
   */
  void Solve(const TridiagonalMatrix &matrix, 
             const double *pdRHS, 
             double *pdX) const;


protected:

  /// The size of the tridiagonal matrices to be solved
  mutable size_t m_nMaxUnknowns;

  /// Work array used by the solving algorithm
  mutable Array<double> m_pdWorkArray;

  /// Boolean indicates if we should soving M^T X = B instead of MX = B
  bool m_bIsTranposed;


private:
 
  // Cannot copy or assign
  NO_COPY_CLASS(TridiagonalSolver);

}; // class TridiagonalSolver


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_TRIDIAGONALSOLVER_H_

