/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/gmressolver.h
// Purpose:     Generalized Minimum Residual method (GMRES) solver
// Author:      Wang
// Created:     2004/04/27
// RCS-ID:      $Id: gmressolver.h,v 1.3 2005/01/17 13:42:13 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_GMRESSOLVER_H_
#define _ITO33_NUMERIC_GMRESSOLVER_H_

#include "ito33/mv/mvmd.h"
#include "ito33/numeric/gmres.h"

namespace ito33
{

namespace numeric
{


/**
   Wrapp the gmres function into a class.

   Useful for a serie of problems of same kind.
 */
class GmresSolver
{ 

public: 

  /**
     Ctor initializes some common parameters

     @param nNbUnknowns the number of unknowns of the linear system
     @param nRestart Restart number of the gmres solver
     @param nMaxIteration Maximum number of iterations
     @param dTolerance the tolerance required for the solving
   */
  GmresSolver(size_t nNbUnknowns, size_t nRestart = 35, 
              size_t nMaxIteration = 1000, 
              double dTolerance = 1.e-12)
            : m_nNbUnknowns(nNbUnknowns),
              m_nRestart(nRestart),
              m_nMaxIteration(nMaxIteration),
              m_dTolerance(dTolerance),
              m_Hessenberg(m_nRestart + 1, m_nRestart),
              m_nNbIters(0)
  {
  }

  // Default dtor is ok
 
  template <class Matrix, class Preconditioner>
  void Solve(const Matrix &matrix, const Preconditioner &preconditioner,
             const double *pdB, double *pdX) const
  {
    CVecteurDouble 
      B(const_cast<double *>(pdB), m_nNbUnknowns, MV_Vector_::ref),
      X(pdX, m_nNbUnknowns, MV_Vector_::ref);

    size_t
      nRestart = m_nRestart,
      nMaxIteration = m_nMaxIteration;
    double
      dTolerance = m_dTolerance;

    m_nNbIters = 0;

    GMRES(matrix, X, B, preconditioner, m_Hessenberg, 
          nRestart, nMaxIteration, dTolerance);

    m_nNbIters = nMaxIteration;
  }  

  size_t GetNbIters() const { return m_nNbIters; }
  

private:

  size_t m_nNbUnknowns;
  
  size_t m_nRestart;
  
  size_t m_nMaxIteration;
  
  double m_dTolerance;

  mutable size_t m_nNbIters;

  mutable MATRIX_double m_Hessenberg;

}; // class GmresSolver


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_GMRESSOLVER_H_
