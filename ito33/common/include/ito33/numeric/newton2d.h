///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/newton2d.h                                //
// Purpose:          Newton2d method                                         // 
// Author:           Ito33                                                   //
// Created:          2004/12/16                                              //
// RCS-ID:           $Id: newton2d.h,v 1.6 2005/12/02 14:41:00 wang Exp $                                                  //
// Copyright (C) Trilemma LLP  2000 -                                        //
///////////////////////////////////////////////////////////////////////////////

/**
   @file  ito33/numeric/newton2d.h

   @brief Newton-Raphson for two variables
  
   It is rather a quick method, 
   but it needs the derivative of the function
 */
#ifndef _ITO33_NUMERIC_NEWTON2D_H_
#define _ITO33_NUMERIC_NEWTON2D_H_

#include "ito33/beforestd.h"
#ifdef NEWTON2DDEBUG
  #include <iostream>
#endif
#include <cmath>
#include "ito33/afterstd.h"



#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/numericerror.h"

extern const ito33::Error ITO33_NEG_TOL;

namespace ito33
{

namespace numeric
{

class Newton2D
{
public:

  /**
     Ctor  
     
     @param dTolerance is optionnal, and gives the tolerance expected

     @param nNbMaxIter is optional, and gives the maximum number of iterations

     @param dParam1LowerBound first parameter lower bound
     @param dParam1UpperBound first parameter upper bound
     @param dParam2LowerBound second parameter lower bound
     @param dParam2UpperBound second parameter lower bound

     //when there are no boundaries set the value to a
     //large number i.e. 1.e10;
   */
  Newton2D(double dParam1LowerBound,
         double dParam1UpperBound,
         double dParam2LowerBound,
         double dParam2UpperBound,
         double dTolerance = 1e-6, 
         size_t nNbMaxIter = 15)
            : m_dParam1LowerBound(dParam1LowerBound),
              m_dParam1UpperBound(dParam1UpperBound),
              m_dParam2LowerBound(dParam2LowerBound),
              m_dParam2UpperBound(dParam2UpperBound),
              m_dTolerance(dTolerance),
              m_dInvTolerance(1./m_dTolerance),
              m_nNbMaxIter(nNbMaxIter), 
              m_nNbIter(0)             
  {
    
    if ( dTolerance <= 0.0 )
      throw EXCEPTION_MSG
            (
              ITO33_NEG_TOL,
              TRANS("Negative tolerance in Newton solver.")
            );
            
  }

  /**
     @param func is a template function that must take 
     dParam1,dParam2,dPrice1,dPrice2,dF,dPrice1Error,dPrice2Error
     7 parameters:
          dParam1 current value of first variable
          dParam2 current value of second variable
          &dPrice1 first cds price
          &dPrice2 second cds price
          &dF      value of the objective function

     @return numerical error
   */
  template<typename T>
  ito33::numeric::NumericError operator()(T& func, double &dParam1, 
                                          double &dParam2)
  {
    //reset values  
    m_nNbIter = 0;
        
    double dPrice1       = 0.;
    double dPrice2       = 0.;
    double dDeriv1Param1 = 0.;    
    double dDeriv2Param1 = 0.; 
    double dDeriv1Param2 = 0.; 
    double dDeriv2Param2 = 0.;
    double dF            = 1.;
    double dParam1Range  = m_dParam1UpperBound - m_dParam1LowerBound;
    double dParam2Range  = m_dParam2UpperBound - m_dParam2LowerBound;
    double dZeroTol = 1.e-5;

    // These are the parameter updates at each iteration. They are used to 
    // determine finite difference shifts. Typically, one would use the 
    // square root of machine epsilon for the shift (assume magnitude 1).
    // However, since we compute derivatives by finite differences, we
    // are more interested in the actual function obtained after the
    // parameter update. In other words, if the parameter is going to be 
    // shifted by 0.1, we are interested in the value F(x+0.1) instead
    // of the actual derivative of F(x) at x. Of course, we don't know
    // what update is going to be computed. It is hoped that a scaled value
    // of the previous update is a good proxy.  To start, use a scaled
    // value of the initial parameter value. Also, limit the updates
    // so the shifts cannot obtain outrageous values.
    double dUpdate1      = fabs(dParam1) * 1.e-3;
    double dUpdate2      = fabs(dParam2) * 1.e-3;
    double dMaxUpdate1   = dParam1Range * 0.2;
    double dMaxUpdate2   = dParam2Range * 0.2;
    
    if (dUpdate1 < 1.e-7)
      dUpdate1 = 1.e-7;

    if (dUpdate2 < 1.e-7)
      dUpdate2 = 1.e-7;

    // In case the bounds are set to the same value, add a small value
    dParam1Range += 1.e-14;
    dParam2Range += 1.e-14;

    // Make sure we can do forward differencing
    if ( dParam1 > m_dParam1UpperBound - m_dTolerance )
      dParam1 = m_dParam1UpperBound - m_dTolerance;

    if ( dParam2 > m_dParam2UpperBound - m_dTolerance)
      dParam2 = m_dParam2UpperBound - m_dTolerance;


    // Compute the objective function, and init all the pricing data
    // with the initial guess
    func(dParam1,dParam2,dPrice1,dPrice2,dF);

    while ( dF  > m_dTolerance )
    { 
      m_nNbIter++;

      // Computing the regular price and objective function is
      // done before the while loop, and at the end of the while loop

      double dShift;
      //shift param1
      {
        dShift =  0.25 * fabs(dUpdate1);
        if (dShift < 1.e-7)
          dShift = 1.e-7;
        if ( dParam1 + dShift > m_dParam1UpperBound )
          dShift = m_dTolerance;

        func(dParam1 + dShift,dParam2,dDeriv1Param1,dDeriv2Param1,dF);

        dDeriv1Param1 = (dDeriv1Param1 - dPrice1)/dShift;
        dDeriv2Param1 = (dDeriv2Param1 - dPrice2)/dShift;

        if (fabs(dDeriv1Param1) < dZeroTol)
          dDeriv1Param1 = 0.0;

        if (fabs(dDeriv2Param1) < dZeroTol)
          dDeriv2Param1 = 0.0;
      }

      //shift param2
      {
        dShift = 0.25 * fabs(dUpdate2);
        if (dShift < 1.e-7)
          dShift = 1.e-7;
        if ( dParam2 + dShift > m_dParam2UpperBound )
          dShift = m_dTolerance;

        func(dParam1,dParam2 + dShift,dDeriv1Param2,dDeriv2Param2,dF);

        dDeriv1Param2 = (dDeriv1Param2 - dPrice1)/dShift;
        dDeriv2Param2 = (dDeriv2Param2 - dPrice2)/dShift;

        if (fabs(dDeriv1Param2) < dZeroTol)
          dDeriv1Param2 = 0.0;

        if (fabs(dDeriv2Param2) < dZeroTol)
          dDeriv2Param2 = 0.0;
      }

      // Check if one of the parameters does not affect the price, or if 
      // one of the prices is unaffected by the parameters. Also
      // avoids a divide by zero.
      if ( fabs(dDeriv1Param2) < dZeroTol && fabs(dDeriv2Param2) < dZeroTol)
      {
        dUpdate1 = -dPrice1/dDeriv1Param1;
        dUpdate2 = 0.0;
      }
      else if (fabs(dDeriv1Param1) < dZeroTol && fabs(dDeriv2Param1) < dZeroTol)
      {
        dUpdate2 = -dPrice2/dDeriv2Param2;
        dUpdate1 = 0.0; 
      }
      else if (fabs(dDeriv1Param1) < dZeroTol && fabs(dDeriv1Param2) < dZeroTol)
      {
        dUpdate1 = -dPrice2/dDeriv2Param1;
        dUpdate2 = 0.0;
      }
      else if (fabs(dDeriv2Param1) < dZeroTol && fabs(dDeriv2Param2) < dZeroTol)
      {
        dUpdate1 = -dPrice1/dDeriv1Param1;
        dUpdate2 = 0.0;
      }
      else
      {
        double dDet = dDeriv1Param1*dDeriv2Param2 
                    - dDeriv1Param2*dDeriv2Param1 ;
     
        dDet = 1. / dDet;

        double dParam1New = dDeriv1Param2*dPrice2
                          - dDeriv2Param2*dPrice1;
        double dParam2New = dDeriv2Param1*dPrice1
                          - dDeriv1Param1*dPrice2;

        dUpdate1 = dParam1New * dDet;        
        dUpdate2 = dParam2New * dDet;
      }

      // Limit the updates.  Rationale: If the previous guess was so bad 
      // that the parameter estimate dramatically changes, the derivative
      // estimation was probably wrong. Make sure we don't over-compensate.
      if ( dUpdate1 > dMaxUpdate1)
        dUpdate1 = dMaxUpdate1;
      else if ( -dUpdate1 > dMaxUpdate1)
        dUpdate1 = -dMaxUpdate1;

      dParam1 += dUpdate1;

      if ( dUpdate2 > dMaxUpdate2)
        dUpdate2 = dMaxUpdate2;
      else if ( -dUpdate2 > dMaxUpdate2)
        dUpdate2 = -dMaxUpdate2;

      dParam2 += dUpdate2;


      //TODO: add linesearch method

      // check value constraints. Don't let values get too high
      // since we (forward) shift to get derivatives
      if ( dParam1 > m_dParam1UpperBound - m_dTolerance )
        dParam1 = m_dParam1UpperBound - m_dTolerance;

      if ( dParam1 < m_dParam1LowerBound )
        dParam1 = m_dParam1LowerBound;

      if ( dParam2 > m_dParam2UpperBound - m_dTolerance)
        dParam2 = m_dParam2UpperBound - m_dTolerance;
 
      if ( dParam2 < m_dParam2LowerBound )
        dParam2 = m_dParam2LowerBound;
 
 
      //compute the objective function with the new value
      func(dParam1,dParam2,dPrice1,dPrice2,dF);

      # ifdef NEWTON2DDEBUG
        std::cout << std::endl << " Iteration #" <<  m_nNbIter << std::endl;

        std::cout << "dParam1 = " << dParam1 << std::endl;
        std::cout << "dParam2 = " << dParam2 << std::endl;

        std::cout << "Current Price Deriv1 = " << dPrice1 << std::endl;
        std::cout << "Current Price Deriv2 = " << dPrice2 << std::endl;

        std::cout << "dF = " << dF << std::endl;
      # endif

      //bound the maximum number of 
      //NewtonRaphson iterations
      if ( m_nNbIter > m_nNbMaxIter  )
        return ITO33_TOO_MANY_ITERATION;

    } //end while

    return ITO33_NO_ERROR;
  } //end newton

  /**
    Get the number on Newton's iteration
  */
  size_t GetNbIter() const 
  { 
    return m_nNbIter; 
  }


private: 

  ///Tolerance for convergence
  double m_dTolerance;
  
  ///Inverse tolerance
  double m_dInvTolerance;
  
  ///Maximum Number of Iteration
  size_t m_nNbMaxIter;
  
  ///Current number of Iteration
  size_t m_nNbIter;

  ///Bounds for the different parameters
  double m_dParam1LowerBound;
  double m_dParam1UpperBound;
  double m_dParam2LowerBound;
  double m_dParam2UpperBound;

};


}  // namespace numeric 

}  // namespace ito33


#endif // #ifndef _ITO33_NUMERIC_NEWTON2D_H_
