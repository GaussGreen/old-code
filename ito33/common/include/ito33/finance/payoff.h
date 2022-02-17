/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/payoff.h
// Purpose:     payoff interface 
// Author:      WANG Xuewen
// Created:     2003/09/25
// RCS-ID:      $Id: payoff.h,v 1.8 2004/10/28 14:45:24 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/payoff.h
    @brief The declaration of the payoff class.   
 */

#ifndef _ITO33_FINANCE_PAYOFF_H_
#define _ITO33_FINANCE_PAYOFF_H_

#include <cstddef> // for size_t

namespace ito33 
{

namespace finance 
{


/// Payoff is the general interface for different kinds of payoff
class Payoff 
{
public:

  /// dummy virtual destructor 
  virtual ~Payoff() { }
  
  /**
     Get the payoff values at an array of spots.

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  virtual void Get(const double *pdS, double *pdPrices, size_t nNbS) const = 0;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PAYOFF_H_

