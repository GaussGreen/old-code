/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurveleg.h
// Purpose:     YieldCurveLeg class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: yieldcurveleg.h,v 1.2 2006/08/21 14:52:27 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurveleg.h
    @brief The declaration of YieldCurveLeg class.
 */

#ifndef ITO33_FINANCE_YIELDCURVELEG_H_
#define ITO33_FINANCE_YIELDCURVELEG_H_

#include "ito33/date.h"

#include "ito33/finance/timeunit.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    yield curve leg
 */
class ITO33_DLLDECL YieldCurveLeg
{
public:

  /// creates a yield curve leg
  YieldCurveLeg( double dRate,
                size_t nMaturityDuration,
                TimeUnit maturityUnit
                );

  /**
      The interest rate. This is can be interbank cash rate, swap rate or
      annually compounded rate etc. with respect to the type of the yield
      curve.

      @return interest rate
   */
  double GetRate() const
  {
    return m_dRate;
  }

  /**
      This is maturity expressed as a number of TimeUnits.

      @return the first part of the maturity expressed as 'Number of TimeUnit'
   */
  size_t GetMaturityDuration() const
  {
    return m_nMaturityDuration;
  }

  /**
      The unit in which the maturity is expressed.

      @return 'TimeUnit' part of the maturity expressed as 'Number of TimeUnit'.
   */
  TimeUnit GetMaturityUnit() const
  {
    return m_maturityUnit;
  }


protected:
  
  /// interest rate
  double m_dRate;

  /// duration part of maturity
  size_t m_nMaturityDuration;

  /// time unit part of maturity
  TimeUnit m_maturityUnit;

};

} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_YIELDCURVELEG_H_

