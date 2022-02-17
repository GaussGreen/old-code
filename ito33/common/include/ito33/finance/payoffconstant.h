/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/payoffconstant.h
// Purpose:     Payoff class defined by a constant value
// Author:      David
// Created:     2004/02/13
// RCS-ID:      $Id: payoffconstant.h,v 1.4 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/payoffconstant.h
    @brief Payoff class defined by a constant value
 */

#ifndef _ITO33_FINANCE_PAYOFFCONSTANT_H_
#define _ITO33_FINANCE_PAYOFFCONSTANT_H_

#include "ito33/finance/payoff.h"

namespace ito33 
{

namespace finance 
{

/**
  Payoff defined by a constant value

  This payoff is primarily used for testing the rate modifications, but may 
  also be useful for pricing bonds.
*/
class PayoffConstant : public Payoff
{
public:

  /**
    constructor, set the constant value
  */
  PayoffConstant(double dPayoffValue)
    : m_dPayoffValue(dPayoffValue)
  {
  }
 
  /// dummy virtual destructor
  virtual ~PayoffConstant() { }

  /**
     Get the payoff values at an array of spots. The first parameter
     (the spots) is not required.
    
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double* /*pdS*/, double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
      pdPrices[nIdx] = m_dPayoffValue;
  }

protected:

  /// The constant value
  double m_dPayoffValue;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PAYOFFCONSTANT_H_
