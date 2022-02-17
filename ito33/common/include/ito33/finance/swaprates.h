/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/swaprates.h
// Purpose:     SwapRates class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: swaprates.h,v 1.3 2006/08/22 21:39:57 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/swaprates.h
    @brief The declaration of SwapRates class.

    Class for SwapRates.
 */

#ifndef ITO33_FINANCE_SWAPRATES_H_
#define ITO33_FINANCE_SWAPRATES_H_

#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/finance/swaprate.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{


/**
    Swap rate Terms.

    @iterator GetAll
 */
class ITO33_DLLDECL SwapRates
{
public:

  /// creates cash rate terms with basis 30/360
  SwapRates() : m_dccSwap(Date::DayCountConvention_30360) {}

  /// The day count convention used for cash rate
  void SetBasis(Date::DayCountConvention dcc)
  {
    Validate(dcc);

    m_dccSwap = dcc;
  }
  
  /**
      Defines a swap rate for a holding period and payment frequency.
      The time unit for the holding period must be "Month" or "Year".
   */
  void AddLeg(double dRate,
                   size_t nMaturityDuration,
                   TimeUnit maturityUnit,
                   Frequency paymentFrequency
                   )
  {
    m_pSwapRates.push_back(SwapRate(dRate,
                                    nMaturityDuration,
                                    maturityUnit,
                                    paymentFrequency));
  }

  /**
      A version with a more consistent signature with others.

      @noexport
   */
  void AddLeg(size_t nMaturityDuration,
              TimeUnit maturityUnit,
              Frequency paymentFrequency,
              double dRate)
  {
    m_pSwapRates.push_back(SwapRate(dRate,
                                    nMaturityDuration,
                                    maturityUnit,
                                    paymentFrequency));
  }

  /// The day count convention used for cash rate
  Date::DayCountConvention GetBasis() const
  {
    return m_dccSwap;
  }

  /**
      Gets all SwapRate objects.

      @noexport COM
   */
  const std::vector<SwapRate>& GetAll() const
  {
    return m_pSwapRates;
  }

private:

  /// Swap basis
  Date::DayCountConvention m_dccSwap;

  /// Swap rates
  std::vector<SwapRate> m_pSwapRates;
};



} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_SWAPRATES_H_

