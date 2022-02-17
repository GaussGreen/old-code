///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/blackscholes.h                            //
// Purpose:          Usual tools in black-scholes PDE                        // 
// Author:           ZHANG Yunzhi                                            //
// Created:          ??                                                      //
// RCS-ID:           $Id: blackscholes.h,v 1.16 2006/01/09 15:43:42 wang Exp $                                                    //
// Copyright (C) Trilemma LLP  2000 -                                               //
///////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/blackscholes.h

   @brief This file contains a "BlackScholes" class. Its members are the usual
   Black-Scholes parameters : Spot, Strike, Rate, ForeignRate, Maturity (in
   fact relative maturity) and Volatility, for a Call or a Put. It also contains
   member functions that call functions defined in blackscholesutils.
 */

#ifndef _ITO33_FINANCE_BLACHSCHOLES_H_
#define _ITO33_FINANCE_BLACHSCHOLES_H_

#include "ito33/numeric/normaldistribution.h"
#include "ito33/finance/blackscholesutils.h"

namespace ito33
{
namespace finance
{

/**
   Contains usual Black-Scholes quantities for call-put options. 
   Member functions give price, vega and delta.
 */
class BlackScholes
{
protected:
  
  double
    m_dSpot,       ///< The spot price
    m_dMaturity,   ///< The relative maturity (maturity - initial time)
    m_dRate,       ///< Local interest Rate
    m_dForeignRate,///< Foreign interest rate
    m_dStrike,     ///< Strike price
    m_dVolatility; ///< Volatility
    
  int
    m_iCallPut;    ///< Indicator of option type : call or put

public:

  /// brief Default constructor
  BlackScholes() { }
 
  /// virtual destructor  
  virtual ~BlackScholes() { }

  /**
     Overloaded constructor where all member variables are initialised.
    
     @param iCallPut If 1, the option is a Call, a put for 0.
     @param dSpot is the Spot
     @param dStrike is the strike price
     @param dRate is the local interest rate
     @param dForeignRate is the foreign interest rate
     @param dVolatility is the volatility
     @param dMaturity is the reltive maturity (maturity - initial time)
   */
  BlackScholes(int iCallPut, double dSpot, double dStrike, double dRate, 
               double dForeignRate, double dVolatility, double dMaturity)
             : m_iCallPut(iCallPut),
               m_dSpot(dSpot), m_dStrike(dStrike), 
               m_dRate(dRate), m_dForeignRate(dForeignRate),
               m_dVolatility(dVolatility), 
               m_dMaturity(dMaturity) { }

  /**
     Assigns m_dSpot with:
     @param dSpot a spot
   */
  void SetSpotSharePrice(double dSpot) { m_dSpot = dSpot; }

  /**
     Assigns m_dStrike with :
     @param dStrike a strike
   */
  void SetStrike(double dStrike) {m_dStrike = dStrike;}

  /**
     Assigns m_dRate with:
     @param dRate an interest rate
   */
  void SetRate(double dRate) {m_dRate = dRate;}

  /**
     Assigns m_dForeignRate with:
     @param dForeignRate a foreign interest rate
   */
  void SetForeignRate(double dForeignRate) {m_dForeignRate = dForeignRate;}

  /**
     Assigns m_dMaturity with:
     @param dMaturity a relative maturity
   */
  void SetMaturity(double dMaturity) {m_dMaturity = dMaturity;}

  /**
     Assigns m_dVolatility with:
     @param dVolatility a volatility
   */
  void SetVolatility(double dVolatility) {m_dVolatility = dVolatility;}

  /**
     Assigns m_iCallPut:
     @param iCallPut an integer (1 or 0)
   */
  void SetCallPut(int iCallPut) {m_iCallPut = iCallPut;}
  
  /**
     Calls some functions in blackscholesutils.h to compute the price,
     given the member variables.

     @return the price 
   */
  double ComputePrice()
  {
    if (m_iCallPut == 1)
      return CalculateCallPrice(m_dSpot, m_dStrike, m_dVolatility, 
                                m_dMaturity,m_dRate, m_dForeignRate);
    else
      return CalculatePutPrice(m_dSpot, m_dStrike, m_dVolatility, 
                               m_dMaturity, m_dRate, m_dForeignRate);
  }

  /**
     Calls CalculateVegaCall() in blackscholesutils.h to compute the Vega,
     given the member variables.

     @return the Vega 
   */
  double ComputeVega()
  {
    return CalculateVegaCall(m_dSpot, m_dStrike, m_dMaturity, 
                             m_dRate, m_dForeignRate, m_dVolatility);
  }

  /**
     Calls  CalculateDelta() in blackscholesutils.h to compute the delta,
     given the member variables.

     @return the delta 
   */
  double ComputeDelta()
  {
    return CalculateDelta(m_iCallPut, m_dSpot, m_dStrike, m_dRate, 
                          m_dForeignRate, m_dVolatility, m_dMaturity);
  }

};


} // namespace finance

} // namespace ito33


#endif // #ifndef _ITO33_FINANCE_BLACHSCHOLES_H_
