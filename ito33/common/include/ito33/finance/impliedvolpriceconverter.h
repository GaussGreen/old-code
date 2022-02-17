///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/impliedvolpriceconverter.h  
// Purpose:          functions to get price from implied vol and vice versa
// Author:           ITO 33 Canada
// Created:          2005/07/18
// RCS-ID:           $Id: impliedvolpriceconverter.h,v 1.3 2006/08/15 20:50:30 wang Exp $                                                    
// Copyright:        (C) 2005- Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
  
  @file ito33/finance/impliedvolpriceconverter.h

  @brief Given an option with implied volatility
  return the price

 */

#ifndef _ITO33_FINANCE_IMPLIEDVOLPRICECONVERTER_H_
#define _ITO33_FINANCE_IMPLIEDVOLPRICECONVERTER_H_

#include "ito33/finance/option.h"

namespace ito33
{

namespace finance
{

/**
    Gets price from implied volatility.

    @param option Option contract
    @param dImpliedVol the Black-Scholes vol
    @return price corresponding to given Black-Scholes vol
 */
double GetPriceFromImpliedVol(const Option& option, double dImpliedVol);

/**
    Gets vega from implied volatility.

    @param option Option contract
    @param dImpliedVol the Black-Scholes vol
    @return vega corresponding to given Black-Scholes implied vol
 */
double GetVegaFromImpliedVol(const Option& option, double dImpliedVol);

/**
    Gets implied volatility from price.

    @param option Option contract
    @param dPrice market price

    @return implied volatility
 */
double GetImpliedVolFromPrice(const Option& option, double dPrice);

} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_IMPLIEDVOLPRICECONVERTER_H_
