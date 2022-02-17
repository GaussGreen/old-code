/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/zerorate.h
// Purpose:     ZeroRate class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: zerorate.h,v 1.2 2005/04/28 15:05:18 cosmin Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/zerorate.h
    @brief The declaration of ZeroRate class.
 */

#ifndef ITO33_FINANCE_ZERORATE_H_
#define ITO33_FINANCE_ZERORATE_H_

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
  ZeroRate class. It represents the zero coupon rate for the given
  actual number of days.

  @nocreate
  */
class ITO33_DLLDECL ZeroRate
{
public:

  /// default constructor doesn't nothing
  ZeroRate() {}

  /// defines a zero rate for given actual number of days
  ZeroRate(double dRate, size_t nActualNumberOfDays)
    : m_dRate(dRate), m_nDuration(nActualNumberOfDays)
  { 
  }

  /// Gets rate value
  double GetRate() const
  {
    return m_dRate;
  }

  /// Gets the duration being the actual number of days till maturity
  size_t GetDuration() const
  {
    return m_nDuration;
  }

#ifndef __CPP2ANY__
  /// compare two ZeroRates
  bool operator<(const ZeroRate& other) const
    { return m_nDuration < other.m_nDuration; }
  /// compare two ZeroRates
  bool operator>(const ZeroRate& other) const
    { return m_nDuration > other.m_nDuration; }
  /// compare two ZeroRates
  bool operator<=(const ZeroRate& other) const
    { return m_nDuration <= other.m_nDuration; }
  /// compare two ZeroRates
  bool operator>=(const ZeroRate& other) const
    { return m_nDuration >= other.m_nDuration; }
#endif
  
private:

  /// rate value
  double m_dRate;

  /// actual number of days till maturity
  size_t m_nDuration;
};

} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_ZERORATE_H_

