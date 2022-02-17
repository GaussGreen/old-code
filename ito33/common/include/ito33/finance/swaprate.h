/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/swaprate.h
// Purpose:     SwapRate class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: swaprate.h,v 1.6 2006/08/22 22:54:42 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/swaprate.h
    @brief The declaration of SwapRate class.

    Class for SwapRate.
 */

#ifndef ITO33_FINANCE_SWAPRATE_H_
#define ITO33_FINANCE_SWAPRATE_H_

#include "ito33/finance/yieldcurveleg.h"
#include "ito33/finance/frequency.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    class for swap rate, the time unit of holding period must be "Month" or
    "Year".

    @nocreate
 */
class ITO33_DLLDECL SwapRate : public YieldCurveLeg
{
public:

  /**
      Creates a swap rate.

      Note that the time unit of holding period must be "Month" or "Year".

      @noexport
   */
  SwapRate( double dRate,
            size_t nMaturityDuration,
            TimeUnit maturityUnit,
            Frequency paymentFrequency
            );

  /// Gets number of months of the maturity term.
  size_t GetNumberMonths() const
  {
    return m_nNumberMonths;
  }

  /// Gets the payment frequency.
  Frequency GetPaymentFrequency() const
  {
    return m_paymentFrequency;
  }

protected:

  /// frequency
  Frequency m_paymentFrequency;

  /// number of months
  size_t m_nNumberMonths;
};


} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_SWAPRATE_H_
