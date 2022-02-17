/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/nonlinearsolver.h
// Purpose:     non linear equation 1D solvers
// Author:      Pedro Ferreira
// Created:     10.08.03
// RCS-ID:      $Id: nonlinearsolver.h,v 1.23 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/nonlinearsolver.h
    @brief Solvers for non linear equations in one variable
    This header defines several solvers for non linear equations
    in one variable.
 */

#ifndef _ITO33_NUMERIC_NONLINEARSOLVER_H_
#define _ITO33_NUMERIC_NONLINEARSOLVER_H_

#include <cmath>

#include "ito33/debug.h"
#include "ito33/error.h"
#include "ito33/gettext.h"
#include "ito33/numeric/exception.h"

/// Declaration of error codes
extern const ito33::Error ITO33_DIFF_PARAM, ITO33_BAD_PARAM,
                          ITO33_MAX_ITER, ITO33_UNEXPECTED,
                          ITO33_NEG_TOL, ITO33_DIV0;

namespace ito33
{
namespace numeric
{

// TODO: move the following several functions to some common numeric header

/// default tolerance for various algorithms
const double DEFAULT_TOLERANCE = 1e-5;

/**
    Check whether the number if 0 up to the given precision.

    @param dValue the value to compare with 0
    @param dTolerance if the abs value of the number is smaller than this
           value it is considered 0
    
    @return true if the number is 0 up to given tolerance
 */
inline bool IsZero(double dValue, double dTolerance = DEFAULT_TOLERANCE)
{
  return fabs(dValue) < dTolerance;
}

/**
    Check whether the given values are of the same sign.

    If one of the values is null this function returns false, see
    HaveSameSignOrZero() if you don't want this.

    @param value1 first value (positive, negative or zero)
    @param value2 second value
   
    @return true if both values are positive or negative, false if they have
            different signs
 */
inline bool HaveSameSign(double value1, double value2)
{
  // alternatively we could do "value1*value2 > 0." but I'm not sure if this
  // is really faster -- to be profiled if anybody cares
  return (value1 > 0.) == (value2 > 0.);
}

/**
    Check whether the given values are of the same sign.

    If one of the values is null this function returns true, use
    HaveSameSign() if you don't want this.

    @param value1 first value (positive, negative or zero)
    @param value2 second value
    
    @return true if both values are positive or negative or one of them is 0,
            false if they have different signs
 */
inline bool HaveSameSignOrZero(double value1, double value2)
{
  // alternatively we could do "value1*value2 >= 0." but I'm not sure if this
  // is really faster -- to be profiled if anybody cares
  return (value1 > 0.) == (value2 > 0.) || IsZero(value1) || IsZero(value2);
}

/** 
   Implement the Secant solver for one dimensional non linear problem.
  
   @todo Add the reference of this method 
 */
class Secant
{
public:
  
  /**
     Constructor for the secant method.

     Constructs a new Secant algorithm. Checks all the values.

     @param dTolerance the tolerance for the result
     @param nNbMaxIter the maximum number of iterations
   */
  Secant(double dTolerance = DEFAULT_TOLERANCE, size_t nNbMaxIter = 1000)
       : m_dTolerance(dTolerance), m_nNbMaxIter(nNbMaxIter),
         m_nNbIter(0)
  {
    if (dTolerance <= 0.0)
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in secant solver")
            );
  }

  /**
     Get the number of iterations.

     @return the number of iterations
   */
  size_t GetNbIter() const { return m_nNbIter; }

  /**
     Find the root of a one variable function using the secant method.

     Finds the zero of a real one variable function by using the secant algorithm.
     By using this algorithm there is no need to have the derivative of the
     function.
     
     @param func the function to find the root
     @param dInitial1 the first initial guess for the zero
     @param dInitial2 the second initial guess for the zero
     
     @return The root of the function
   */
  template<typename F>
  double operator()(F& func, double dInitial1, double dInitial2);


protected:
  
  double m_dTolerance;
  
  size_t m_nNbMaxIter;
  
  size_t  m_nNbIter;

  NO_COPY_CLASS(Secant);

};

/// Implementation of operator ()
template<typename F>
double Secant::operator()(F& func, double dInitial1, double dInitial2)
{
  double dTnm1 = dInitial2;
  double dTnm2 = dInitial1;

  for ( m_nNbIter = 0; m_nNbIter < m_nNbMaxIter; ++m_nNbIter )
  {
    double dResult = func(dTnm1);
    if ( fabs(dResult) <= m_dTolerance )
    {
      // we achieved the wanted precision
      return dTnm1;
    }

    double dDiff = dResult - func(dTnm2);
    if ( fabs(dDiff) < m_dTolerance )
    {
      throw EXCEPTION_MSG
            (
              ITO33_DIV0,
              TRANS("Divide by zero in secant solver")
            );
    }

    double dTn = dTnm1 - dResult * (dTnm1 - dTnm2) / dDiff;
    dTnm2 = dTnm1;
    dTnm1 = dTn;
  }

  throw EXCEPTION_MSG
        (
          ITO33_MAX_ITER,
          TRANS("Maximum number of iterations exceeded in secant solver")
        );
}


/** 
   Implement the Regula Falsi solver for one dimensional non linear problem.
  
   @todo Add the reference of this method 
 */
class RegulaFalsi
{
public:

  /**
     Constructor for the Regula Falsi method.

     Constructs a new Secant algorithm. Checks all the values.

     @param dTolerance the tolerance for the result
     @param nNbMaxIter the maximum number of iterations
     @param nTreshold maximum number of iterations before
            the function f(a) or F(b) on [a,b] starts to get divided
            by 2. Not that at this point the tolerance is divided
            by 10.

      If the function f(x) is strongly curved on the starting interval,
      the approach to the root can be one-sided in the sense that one end
      of the interval remains fixed rather than the endpoint that is discarded
      oscillating from one end to the other. This can be partially suppressed
      by artificially decreasing the value of f(x) at the stagnant end of the
      interval to f(x)/2 at each iteration. This does not alter the sign of 
      f(x) at that endpoint, but does tend to accelerate the convergence 
      towards the root by inducing a more equal distribution of the interval
      end which is dropped. 

   */
  RegulaFalsi(double dTolerance = DEFAULT_TOLERANCE, size_t nNbMaxIter = 1000,
    size_t nTreshold = 15)
            : m_dTolerance(dTolerance),
              m_nNbMaxIter(nNbMaxIter),
              m_nNbIter(0),
              m_nTreshold(nTreshold)

  {
    if (dTolerance <= 0.0)
      throw EXCEPTION_MSG
           (
            ITO33_NEG_TOL,
            TRANS("Negative tolerance in regula falsi solver")
           );
  }

  /**
     Given a point find the matching one which would bracket the root
     between these two points supposing that the function is monotonuous.

     This function may be used to choose the second inital value for
     RegulaFalsi ctor. Notice that this can only work for monotonuous
     function but that we have no way to check for this here so care should
     be taken in the caller code to ensure that this is indeed the case.

     @param func the function to find the root for
     @param dInitial1 the "first" initial point (which may just as well
                       be the second one)

     @return the second initial point for RegulaFalsi ctor
   */
  template<typename F>
  static double FindMatchingRootBracket(F& func, double dInitial1);

  /**
     Get the number of iteration.

     @return number of iterations
  */
  size_t GetNbIter() const { return m_nNbIter; }

  /**
     Find the root of a one variable function using the RegulaFalsi method.

     Finds the zero of a real one variable function by using the RegulaFalsi algorithm.
     By using this algorithm there is no need to have the derivative of the
     function.
     
     @param func the function to find the root
     @param dInitial1 the first initial guess for the zero
     @param dInitial2 the second initial guess for the zero
     
     @return The root of the function
   */
  template<typename F>
  double operator()(F& func, double dInitial1, double dInitial2);


protected:

  double m_dTolerance;

  size_t m_nNbMaxIter;

  size_t m_nNbIter;

  size_t m_nTreshold;

  NO_COPY_CLASS(RegulaFalsi);

};


template<typename F>
double RegulaFalsi::operator()(F& func, double dInitial1, double dInitial2)
{
  double
    dTnm1, dTnm2,
    dDistance,
    dFunc1, dFunc2,
    dStep,
    dRoot, dFRoot;

  double dTolerance = m_dTolerance;

  m_nNbIter = 0;

  if (dInitial1 < dInitial2)
  {
    dTnm1 = dInitial1;
    dTnm2 = dInitial2;
  }
  else
  {
    dTnm1 = dInitial2;
    dTnm2 = dInitial1;
  }

  dDistance = fabs(dTnm2 - dTnm1);

  dFunc1 = func(dTnm1);

  if (fabs(dFunc1) < dTolerance)
    return dTnm1;

  dFunc2 = func(dTnm2);

  if (fabs(dFunc2) < dTolerance)
    return dTnm2;

  if (dDistance < dTolerance )
    throw EXCEPTION_MSG
          (
            ITO33_DIFF_PARAM,
            TRANS("Initial values are equal and are not root in RegulaFalsi")
          );

  if (dFunc1 * dFunc2 > 0.)
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Initial values do not bracket root in RegulaFalsi")
          );

  for ( ; m_nNbIter < m_nNbMaxIter; ++m_nNbIter )
  {


    // dFunc1 and dFunc2 are always opposed in sign, so can't be equal
    dStep = dFunc1 / ( dFunc2 - dFunc1 ); 
    dRoot = dTnm1 - dDistance * dStep;

    if (dRoot < dTnm1 || dRoot > dTnm2)
      throw EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              TRANS("RegulaFalsi solver fail to find the root")
            );

    dFRoot = func(dRoot);

    // treshold has been reached decrease tolerance
    // to force a few more iterations
    if ( m_nNbIter == m_nTreshold )
      dTolerance *= .1;
        
    if ( fabs(dFRoot) <= dTolerance )
    {
      // we achieved the wanted precision
      return dRoot;
    }

    if (dFRoot * dFunc1 < 0.)
    {
      dTnm2 = dRoot;
      dFunc2 = dFRoot;

      if ( m_nNbIter >= m_nTreshold )
        dFunc1 *= .5;
    }
    else
    {
      dTnm1 = dRoot;
      dFunc1 = dFRoot;

      if ( m_nNbIter >= m_nTreshold)
        dFunc2 *= .5;
    }


    dDistance = fabs(dTnm2 - dTnm1);
  }

  throw EXCEPTION_MSG
        (
          ITO33_MAX_ITER,
          TRANS("Maximum number of iterations exceeded in RegulaFalsi solver")
        );
}

template <typename F>
/* static */
double RegulaFalsi::FindMatchingRootBracket(F& func, double dInitial1)
{
  // first find out in which direction we should be looking for the second
  // point
  const double dValue1 = func(dInitial1);
  if ( IsZero(dValue1) )
  {
    // return the same point, the algorithm will notice when it is called with
    // 2 identical points and return the same value as root which is exactly
    // what we want
    return dInitial1;
  }

  double 
    inc = HaveSameSign(func(dInitial1 + 1.) - dValue1, dValue1) ? -1. : 1.;

  // now advance in ever increasing steps until we get to the other side of
  // the root or the max number of iterations is spent without achieving this
  // goal
  const size_t MAX_ITER = 100; // Yet another magic number?

  size_t nIter = 0;

  double dInitial2 = dInitial1 + inc;
  while ( HaveSameSignOrZero(dValue1, func(dInitial2)) )
  {
    inc *= 2.;
    dInitial2 += inc;

    if ( ++nIter > MAX_ITER )
    {
      throw EXCEPTION_MSG
            (
              ITO33_MAX_ITER,
              TRANS("Maximum number of iterations exceeded in"
                    "RegulaFalsi::FindMatchingRootBracket()")
            );
    }
  }

  return dInitial2;
}

}  // namespace numeric

}  // namespace ito33

#endif // _ITO33_NUMERIC_NONLINEARSOLVER_H_

