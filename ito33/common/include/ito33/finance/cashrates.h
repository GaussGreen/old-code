/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashrates.h
// Purpose:     CashRates class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: cashrates.h,v 1.3 2006/08/22 21:39:57 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashrates.h
    @brief The declaration of CashRates class.
 */

#ifndef ITO33_FINANCE_CASHRATES_H_
#define ITO33_FINANCE_CASHRATES_H_

#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/finance/cashrate.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Cash rate Terms.

    @iterator GetAll
 */
class ITO33_DLLDECL CashRates
{
public:

  /// creates cash rate terms with basis 30/360
  CashRates() : m_dccCash(Date::DayCountConvention_30360) {}

  /// The day count convention used for cash rate
  void SetBasis(Date::DayCountConvention dcc)
  {
    Validate(dcc);

    m_dccCash = dcc;
  }
  
  /// defines cash rate for a holding period
  void AddLeg(double dRate,
              size_t nMaturityDuration, TimeUnit maturityUnit)
  {
    m_pCashRates.push_back(CashRate(dRate, nMaturityDuration, maturityUnit));
  }

  /**
      A version with a more consistent signature with others.

      @noexport
   */
  void AddLeg(size_t nMaturityDuration,
              TimeUnit maturityUnit,
              double dRate)
  {
    m_pCashRates.push_back(CashRate(dRate, nMaturityDuration, maturityUnit));
  }

  /// The day count convention used for cash rate
  Date::DayCountConvention GetBasis() const
  {
    return m_dccCash;
  }

  /**
      Gets all CashRate objects.

      @noexport COM
   */
  const std::vector<CashRate>& GetAll() const
  {
    return m_pCashRates;
  }

private:

  /// cash basis
  Date::DayCountConvention m_dccCash;

  /// cash rates
  std::vector<CashRate> m_pCashRates;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef ITO33_FINANCE_CASHRATES_H_
