/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/slidingmaturity.h
// Purpose:     Convert maturity in months to date for a sliding instruemnt
// Created:     2005/07/20
// RCS-ID:      $Id: slidingmaturity.h,v 1.1 2005/07/20 17:19:39 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/date.h"

extern const ito33::Error ITO33_SLIDINGCDS_MATURITY;

namespace ito33
{

namespace finance
{

/**
   Determines the maturity of a sliding instrument.

   @param today Reference date
   @param months rough maturity as a number of months
 */
Date GetSlidingMaturity(Date referenceDate, size_t months)
{
  CHECK_COND(months > 0, ITO33_SLIDINGCDS_MATURITY);

  // The year of today
  size_t nYear = referenceDate.GetYear();

  // Find the next regular date: March 20, June 20, Sep 20, Dec 20
  // What if today is already a regular date? convention?
  Date regularDate;

  if ( referenceDate >= Date(nYear, Date::Dec, 20) )
    regularDate = Date(nYear + 1, Date::Mar, 20);
  else
    for (size_t n = 0; n < 4; n++)
      if ( Date(nYear, Date::Month(3 * ( n + 1)), 20) > referenceDate )
      {
        regularDate = Date(nYear, Date::Month(3 * ( n + 1)), 20);

        break;
      }

  // If there is too few days between the reference date and the current 
  // regular date, we may need to use the regular date next to the current one
  // as first payment date
  // For now, this is not implemented because of missing of specif
  // This brings up also the problem of maturity: Is it from the current 
  // regular date or the first payment date?
  Date firstPaymentDate = regularDate;

  Date maturityDate = firstPaymentDate;
  maturityDate.AddMonths(months);

  return maturityDate;
}


} // namespace finance

} // namespace ito33
