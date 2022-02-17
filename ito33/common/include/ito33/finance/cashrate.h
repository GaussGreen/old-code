/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashrate.h
// Purpose:     CashRate class
// Author:      ZHANG Yunzhi
// Created:     2004/04/11
// RCS-ID:      $Id: cashrate.h,v 1.5 2006/08/22 22:54:42 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashrate.h
    @brief The declaration of CashRate class.
 */

#ifndef ITO33_FINANCE_CASHRATE_H_
#define ITO33_FINANCE_CASHRATE_H_

#include "ito33/finance/yieldcurveleg.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    CashRate class.

    In general, it concerns the maturities ranging from 1 day to 1 year.
    A cash rate r at maturity t means that the interest of 1 dollar will
    be (1 + r) * t.

    @nocreate
 */
class ITO33_DLLDECL CashRate: public YieldCurveLeg
{
public:

  /**
      Create a cash rate leg.

      For example, a 1-month cash rate of 2.5%.

      @param the cash rate with the given maturity. In the example, this
             is 0.025.
      @param nMaturityDuration the maturity as number of TimeUnit. In the
             example, this is 1.
      @param maturityUnit the TimeUnit in which the maturity is expressed.
             In the example, this is TimeUnit_Month.
             
      @noexport
   */
  CashRate( double dRate,
            size_t nMaturityDuration,
            TimeUnit maturityUnit
            )
          : YieldCurveLeg(dRate, nMaturityDuration, maturityUnit)
  {
  }
};

} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_CASHRATE_H_
