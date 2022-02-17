///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/exoticoption/curran.h   
// Purpose:          Curran's formula to price Asian Option
//                   very common for practitionner, details are
//                   in Risk
// Author:           ITO 33 Canada                                       
// Created:          April 27, 2005                                                  
// RCS-ID:           $Id: curran.h,v 1.5 2006/03/24 16:12:52 yann Exp $           
// Copyright (C) 2005 - Trilemma LLP  
///////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/exoticoption/curran.h

   @brief Curran Formula to compute discretely observed
          Asian option using some approximation. This formula
          is very often used on the street for simple cases.
          "Beyond Average Intelligence", Michael Curran (1992),
          Risk vol. 5, number. 10, pp. 60. 
          (or over the rainbow/RISK Publications)
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_CURRAN_H_
#define _ITO33_FINANCE_EXOTICOPTION_CURRAN_H_

namespace ito33
{

namespace finance
{

/**
  Price of a fixed strike asian call option
  using curran's formula. No hazard rate

  @param dSpot          spot price
  @param dStrike        strike
  @param dMaturityTime  maturity time from today in years
  @param nNbAveraging   number of averaging points
  @param dStartTime     time of the first averaging points
  @param dVolatility    volatility
  @param dInterestRate  rate Interest rate continuous rate
  @param dDividend      continuous dividends
  @param dA             current average
  @param nNbIntoAveraging number of points into averaging

  @return call price

  @noexport
*/
double CurranCall(double dSpot, 
                  double dStrike, 
                  double dMaturityTime,
                  size_t nNbAveraging,
                  double dStartTime, 
                  double dVolatility, 
                  double dInterestRate, 
                  double dDividend,
                  double dA=0,
                  size_t nNbIntoAveraging=0);


/**
  Price of a fixed strike asian put option

  @param dSpot          spot price
  @param dStrike        strike
  @param dMaturityTime  maturity time from today in years
  @param nNbAveraging   number of averaging points
  @param dStartTime     time of the first averaging points
  @param dVolatility    volatility
  @param dInterestRate  rate Interest rate continuous rate
  @param dDividend      continuous dividends
  @param dA             current average
  @param nNbIntoAveraging number of points into averaging

  @return put price

  @noexport
*/
double CurranPut(double dSpot, 
                  double dStrike, 
                  double dMaturityTime,
                  size_t nNbAveraging,
                  double dStartTime, 
                  double dVolatility, 
                  double dInterestRate, 
                  double dDividend,
                  double dA=0,
                  size_t nNbIntoAveraging=0);

/**
  Price of a fixed strike asian  option

  @param nIsCall      1 is put -1 if call
  @param dSpot          spot price
  @param dStrike        strike
  @param dMaturityTime  maturity time from today in years
  @param nNbAveraging   number of averaging points
  @param dStartTime     time of the first averaging points
  @param dVolatility    volatility
  @param dInterestRate  rate Interest rate continuous rate
  @param dDividend      continuous dividends
  @param dA             current average
  @param nNbIntoAveraging number of points into averaging

  @return put price

  @noexport
*/
double PriceCurran(int nIsCall,
                   double dSpot, 
                   double dStrike, 
                   double dMaturityTime,
                   size_t nNbAveraging,
                   double dStartTime, 
                   double dVolatility, 
                   double dInterestRate, 
                   double dDividend,
                   double dA,
                   size_t nNbIntoAveraging);
                   
/**
  Check that the parameters are valid

  @param dSpot          spot price
  @param dStrike        strike
  @param dMaturityTime  maturity time from today in years
  @param nNbAveraging   number of averaging points
  @param dStartTime     time of the first averaging points
  @param dVolatility    volatility
  @param dInterestRate  rate Interest rate continuous rate
  @param dDividend      continuous dividends
  @param dA             current average
  @param nNbIntoAveraging number of points into averaging
*/
void CheckParam(double dSpot, 
                double dStrike, 
                double dMaturityTime,
                size_t nNbAveraging,
                double dStartTime, 
                double dVolatility, 
                double dInterestRate, 
                double dDividend,
                double dA,
                size_t nNbIntoAveraging);

} // namespace finance

} // namespace ito33


#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_CURRAN_H_
