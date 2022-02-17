///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/bisecnewton.h                             //
// Purpose:          mix of  newton method and secant one                    // 
// Author:           Laurence                                                //
// Created:          29/09/2003                                              //
// RCS-ID:           $Id: bisecnewton.h,v 1.12 2006/03/23 09:27:14 yann Exp $                                                    //
// Copyright (C) Trilemma LLP  2000 -                                        //
///////////////////////////////////////////////////////////////////////////////

/**
   @file  ito33/numeric/bisecnewton.h

   @brief Hybrid method between Newton-Raphson and bisection as in Numerical 
          Recipes p366. 
  
   It is rather a quick method, but it needs the derivative of the function
 */
#ifndef _ITO33_NUMERIC_BISECNEWTON_H_
#define _ITO33_NUMERIC_BISECNEWTON_H_

#include <cmath>

#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

/**
   Declaration of error codes
*/
//@{

extern const ito33::Error ITO33_DIFF_PARAM;
extern const ito33::Error ITO33_BAD_PARAM;
extern const ito33::Error ITO33_MAX_ITER;
extern const ito33::Error ITO33_NEG_TOL;

//@}

namespace ito33
{

namespace numeric
{


class BisecNewton
{
public:

  class Exception : public ito33::numeric::Exception
  {
  public:
    /**
        Ctor for the Exception object.

        Use the standard EXCEPTION macro to create Exception objects, this
        frees you from having to type __FILE__, __LINE__ and __FUNCTION__
     */
    Exception(int errorCode,
         const char *message,
         const char *filename,
         size_t line,
         const char *function)
       : ito33::numeric::Exception(errorCode, message, filename, line, function)
      {
      }
  };

  /**
     Ctor  
     
     @param dTolerance is optionnal, and gives the tolerance expected

     @param nNbMaxIter is optional, and gives the maximum number of iterations
   */
  BisecNewton(double dTolerance = 1e-8, size_t nNbMaxIter = 1000)
            : m_dTolerance(dTolerance), m_nNbMaxIter(nNbMaxIter),
              m_dInitialFuncValue1(0), m_dInitialFuncValue2(0), m_nNbIter(0),
              m_bInitialGuessSet(false)
  {
    if (dTolerance <= 0.0)
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in secant solver")
            );
  }

  /// Sets initial guess
  void SetInitialGuess(double dInitialGuess)
  {
    m_bInitialGuessSet = true;
    m_dInitialGuess = dInitialGuess;
  }

  /// clear initial guess
  void ClearInitialGuess()
  {
    m_bInitialGuessSet = false;
  }

  /**
     @param func is a template function that must take 3 parameters, two 
            of them as references. The function is the one to find the zero for 
            first variable parameter given. The second (and referenced) 
            parameter must be the result of the function at that point, 
            the third (and also referenced), its derivative

     @param dInitial1 is a boundary of the interval in which we expect
                     to find the solution

     @param dInitial2 is the other boundary
   */
  template<typename T>
  double operator()(T& func, double dInitial1, double dInitial2)
  {   
    //Reset some values to flush out result of a possible previous call
    m_dInitialFuncValue1 = 0;
    m_dInitialFuncValue2 = 0;
    m_nNbIter = 0;

    double 
      dInf,
      dSup,
      dFuncInf,        
      dFuncSup,         
      dDerivInf,        
      dDerivSup,        
      dRoot,
      dFuncRoot,        
      dDerivRoot,       
      dStep;

    m_nNbIter = 0;

    if (dInitial1 < dInitial2) 
    {
      dInf = dInitial1;
      dSup = dInitial2; 
    }
    else
    {
      dInf = dInitial2; 
      dSup = dInitial1; 
    }

    double dDistance = fabs(dSup - dInf);

    func(dInf, dFuncInf, dDerivInf);
    
    // store the value
    m_dInitialFuncValue1 = dFuncInf;

    if (fabs(dFuncInf) < m_dTolerance)
      return dInf;

    // store the value
    func(dSup, dFuncSup, dDerivSup);

    m_dInitialFuncValue2 = dFuncSup;

    if (fabs(dFuncSup) < m_dTolerance)
      return dSup;

    if (dDistance < m_dTolerance ) 
      throw EXCEPTION_MSG
            (
              ITO33_DIFF_PARAM,
              TRANS("Initial values are equal and are not root in BisecNewton")
            );

    if (dFuncInf * dFuncSup > 0.)
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Initial values do not bracket root in BisecNewton")
            );

    if(m_bInitialGuessSet)
      dRoot = m_dInitialGuess;
    else
      dRoot =  0.5 * (dInf + dSup);
    
    func(dRoot, dFuncRoot, dDerivRoot);

    if (fabs(dFuncRoot) < m_dTolerance)
      return dRoot;

    while (fabs(dFuncRoot) > m_dTolerance)
    {
      ++m_nNbIter;

      if (m_nNbIter > m_nNbMaxIter)
        throw EXCEPTION_MSG
              (
                ITO33_MAX_ITER,
                TRANS("Iteration number greater than max in BisecNewton")
              );

      double dTgInf = (dRoot - dInf) * dDerivRoot - dFuncRoot;
      double dTgSup = (dRoot - dSup) * dDerivRoot - dFuncRoot;
      double dParam1 = fabs(dFuncRoot * 2.);
      double dParam2 = fabs(dDistance * dDerivRoot);

      if 
        (   dTgInf * dTgSup > 0.       // tangent nullifies out of interval
        || dParam1 > dParam2          // bisection more efficient than Newton  
        || fabs(dDerivRoot) < 1e-15   )     // newton could explode
      {
        // bisection method
        dStep = .5 * dDistance;
        dRoot = dInf + dStep;
      }
      else
      {   
        // Newton method
        dStep = dFuncRoot / dDerivRoot;
        dRoot -= dStep;
      }

      func(dRoot, dFuncRoot, dDerivRoot);

      if (dFuncRoot * dFuncInf < 0.)
        dSup = dRoot;
      else
        dInf = dRoot;
      
      dDistance = dSup - dInf;  
    }
       
    return dRoot;
  }

  double GetInitialFuncValue1() const { return m_dInitialFuncValue1; }

  double GetInitialFuncValue2() const { return m_dInitialFuncValue2; }

  size_t GetNbIter() const { return m_nNbIter; }


private: 

  double m_dTolerance;
  size_t m_nNbMaxIter;
 
  double m_dInitialFuncValue1;
  double m_dInitialFuncValue2; 
  
  size_t m_nNbIter;

  bool m_bInitialGuessSet;
  double m_dInitialGuess;
};


}  // namespace numeric 

}  // namespace ito33


#endif // #ifndef _ITO33_NUMERIC_BISECNEWTON_H_
