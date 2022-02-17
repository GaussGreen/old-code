/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/rateutils.h
// Purpose:     Some useful function for the rates
// Author:      Nabil
// Created:     2005/06/03
// RCS-ID:      $Id: rateutils.h,v 1.4 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_FINANCE_RATEUTILS_H_
#define _ITO33_FINANCE_RATEUTILS_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/yieldcurve.h"

namespace ito33
{

namespace finance
{

/**
   Actually this is local function called by ComputeFloatingRates.
   Computes the forward reference (ex: LIBOR) rates for a given dates.

   For exmple, for floating coupons,at a given payment date t_i, the reference 
   rate used to compute the payment amount is the forward rate r(t_ref, t_i),
   where t_ref = t_{i-1} - iFixingDelay. 

   @param previousDate the previous date befor the first one of pDates.
   @param pDates the dates
   @param iFixingDelay the fixing delay
   @param pYieldCurve the yield curve

   @return the forward reference rates
*/
std::vector<double>
ComputeReferenceRates(Date previousDate, 
                      const std::vector<Date>& pDates,
                      int iFixingDelay,
                      const shared_ptr<YieldCurve>& pYieldcurve);

/**
   Computes floating rates.
 
   For exmple, for floating coupons,at a given payment date t_i, the reference 
   rate used to compute the payment amount is the forward rate r(t_ref, t_i),
   where t_ref = t_{i-1} - iFixingDelay.
   The coupon rate is:
   c(t_i) = Multiplier * r(t_ref, t_i) + Margin, this value is capped and 
   floored.
    
   @param previousDate the previous date befor the first one of pDates.
   @param pDates the dates
   @param dcc day count convention of the cash rate (libor)
   @param iFixingDelay the fixing delay
   @param dMargin the margin
   @param dMultiplier the multiplier
   @param dCap the cap for the coupon rate
   @param dFloor the floor for the coupon rate
   @param pYieldCurve the yield curve

   @return the coupon rates computed.   
 */
std::vector<double>
ComputeFloatingRates(Date previousDate, 
                   const std::vector<Date>& pDates,
                   Date::DayCountConvention dcc,
                   int iFixingDelay,
                   double dMargin,
                   double dMultiplier,
                   double dCap,
                   double dFloor,
                   const shared_ptr<YieldCurve>& pYieldcurve);

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_RATEUTILS_H_
