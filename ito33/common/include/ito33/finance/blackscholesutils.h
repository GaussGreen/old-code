///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/blackscholesutils.h                       //
// Purpose:          Useful functions for the resolution of black-scholes PDE// 
// Author:           ZHANG Yunzhi                                            //
// Created:          june 2000                                               //
// RCS-ID:           $Id: blackscholesutils.h,v 1.9 2006/03/24 16:20:06 yann Exp $                                                    //
// Copyright (C) Trilemma LLP  2000 -                                               //
///////////////////////////////////////////////////////////////////////////////

/**
  
  @file ito33/finance/blackscholesutils.h

  @brief This file contains the declaration of several functions useful in
  black-scholes PDE. Formulae can be found in Wilmott "Derivatives", p91-113

 */

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------


#ifndef _ITO33_FINANCE_BLACKSCHOLESUTILS_H_
#define _ITO33_FINANCE_BLACKSCHOLESUTILS_H_

namespace ito33
{
namespace finance
{

/**
   Compute the call price from the BS solution given by Wilmott 
         "Derivatives" p 97.

  @param dSpot is the Spot price of the option ("S" in Wilmot)
  @param dStrike is the Strike price of the option ("E" in Wilmott)
  @param dVolatility is the Volatility ("\f$ \sigma \f$ in Wilmott)
  @param dMaturity is the Relative Maturity ("T-t" in BS formulae)
  @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
  @param dForeignRate is the Foreign Continuous rate "D" in 
          Wilmott
 */
double CalculateCallPrice(double dSpot, double dStrike, double dVolatility,
                          double dMaturity, 
                          double dRate, double dForeignRate);

//--------------------------------------------------------------------------

/**

   Compute the put price from the BS solution given by Wilmott
         "Derivatives" p 99.

  @param dMaturity is the Relative Maturity ("T-t" in BS formulae)
  @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
  @param dSpot is the Spot price of the option ("S" in Wilmott)
  @param dStrike is the Strike price of the option ("E" in Wilmott)
  @param dVolatility is the Volatility ("\f$ \sigma \f$ in Wilmott)
  @param dForeignRate is the Foreign Continuous rate "D" in 
  Wilmott

 */
double CalculatePutPrice(double dSpot, double dStrike, double dVolatility, 
                         double dMaturity, double dRate, 
                         double dForeignRate);

//--------------------------------------------------------------------------

/**

   @brief {Calculate the Vega of a Call option in the BS equation from the 
   formula proposed by Wilmott in "Derivatives" p 107. The Vega is the 
   sensitivity of the option price to the volatility :
   \f$ Vega = \frac{\partial V}{\partial \sigma}\f$.}

   @param dSpot is the Spot price of the option ("S" in Wilmott)
   @param dStrike is the Strike price of the option ("E" in Wilmott)
   @param dMaturity is the Relative Maturity ("T-t" in BS formulae)
   @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
   @param dForeignRate is the Foreign Continuous rate "D" in 
   Wilmott
   @param dVolatility is the Volatility ("\f$ \sigma \f$ in Wilmott)
  
 */
double CalculateVegaCall(double dSpot, double dStrike, double dMaturity, 
                         double dRate, double dForeignRate, 
                         double dVolatility);

//--------------------------------------------------------------------------

/**

   @brief {Calculate for a call option, the volatility that matches a given dPrice
   thanks to a Newton method. See Wilmmott "Derivatives" p 110
   It also returns the Vega, witch corresponds to the calculated volatility. If the
   probleme does not have any solution, dPrice = 0 and dVega = 0 are returned.}


   @param dSpot is the Spot price of the option ("S" in Wilmott)
   @param dStrike is the Strike price of the option ("E" in Wilmott)
    @param dPrice is the given price ("V" in Wilmott)
   @param dMaturity is the Relative Maturity ("T-t" in BS formulae)
   @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
   @param dForeignRate is the Foreign Continuous rate "D" in 
   Wilmott
  
   

 */
double CalculateImpliedVolatility(double dSpot, double dStrike, double dPrice, 
                                  double dMaturity, double dRate, 
                                  double dForeignRate);

//--------------------------------------------------------------------------

/**

   @brief {Calculate the Delta of a call/put option of the BS equation from
   the formula given by Wilmott "Derivatives" p 103. The delta is the 
   sensitivity of the option to the underlying :
   \f$ \Delta = \frac{\partial V}{\partial S} \f$.}

   @param iCallPut is equal to 1 when the option is a Call, 0 when it is a put
   @param dSpot is the Spot price of the option ("S" in Wilmott)
   @param dStrike is the Strike price of the option ("E" in Wilmott)
   @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
   @param dForeignRate is the Foreign Continuous rate "D" in 
   Wilmott
   @param dVolatility is the Volatility ("\f$ \sigma \f$ in Wilmott)
   @param dMaturity is the Relative Maturity ("T-t" in BS formulae)

 */
double CalculateDelta(int iCallPut, double dSpot, double dStrike, double dRate,
                      double dForeignRate, double dVolatility, 
                      double dMaturity);


//--------------------------------------------------------------------------


/**

   Get the Strike that corresponds to the Delta of a given Call/put 
   thanks to a Newton method.

   @param dMaturity is the Relative Maturity ("T-t" in BS formulae)
   @param dSpot is the Spot price of the option ("S" in Wilmott)
   @param dDelta is the given Delta if delta < 0, the option is a call
   @param dVolatility is the Volatility ("\f$ \sigma \f$ in Wilmott)
   @param dRate is the Domestical Continuous Rate ("r" in Wilmott)
   @param dForeignRate is the Foreign Continuous rate "D" in 
   Wilmott

 */
double GetStrikeFromDelta(double dMaturity, double dSpot, double dDelta, 
                          double dVolatility, double dRate, 
                          double dForeignRate);


} // namespace finance
} // namespace ito33


#endif // _ITO33_FINANCE_BLACKSCHOLESUTILS_H_
