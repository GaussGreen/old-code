/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/rateutils.cpp
// Purpose:     Some useful function for the rates
// Author:      Nabil
// Created:     2005/06/03
// RCS-ID:      $Id: rateutils.cpp,v 1.4 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/arraycheckers.h"
#include "ito33/exception.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/finance/rateutils.h"

extern const ito33::Error ITO33_BAD_DATA;

namespace ito33
{

namespace finance
{

std::vector<double>
ComputeReferenceRates(Date previousDate, 
                      const std::vector<Date>& pDates,
                      Date::DayCountConvention dcc,
                      int iFixingDelay,
                      const shared_ptr<YieldCurve>& pYieldcurve)
{
  Date referenceDate = previousDate.AddDays(-iFixingDelay);

  CHECK_COND_MSG( referenceDate < pDates[0],
                  ITO33_BAD_DATA,
                  "Reference date after the one of the first unknown rate." );
   
  CHECK_COND_MSG( referenceDate >= pYieldcurve->GetReferenceDate(),
                  ITO33_BAD_DATA,
                  "Reference date used to compute the first unknown rate is "
                  "before the reference date of the yield curve.");
  
  CheckIncreasingOrder
  (
    pDates, 
    "pDates must be a non empty array of increasing dates."
  );
  
  size_t
    nIdxRate,
    nNbRates = pDates.size();
  
  if( iFixingDelay <= 0 )
  {
    for( nIdxRate = 1; nIdxRate < nNbRates; ++nIdxRate )
    {
      CHECK_COND_MSG
      ( 
        Date::DaysDiff( pDates[nIdxRate - 1], pDates[nIdxRate] ) > 
        - iFixingDelay,
        ITO33_BAD_DATA,
        "FixingDelay must be greater than the difference (in days) between "
        "the dates of any two consecutive unknown rates."
      );      
    }
  }
  
  std::vector<double> pdReferenceRates;

  double 
    dReferenceDate,
    dRate,
    dDate;
  
  for(nIdxRate = 0; nIdxRate < nNbRates; ++nIdxRate)
  {    
    dReferenceDate = GetDoubleFrom( referenceDate );

    dDate = GetDoubleFrom( pDates[nIdxRate] );
    
    dRate = ( 1. / pYieldcurve->GetForwardDiscountFactor(dReferenceDate, dDate) 
      - 1. ) / Date::YearsDiff(referenceDate, pDates[nIdxRate], dcc);

    pdReferenceRates.push_back(dRate);

    referenceDate = pDates[nIdxRate];
    referenceDate.AddDays(-iFixingDelay);
  }

  return pdReferenceRates;
}

std::vector<double>
ComputeFloatingRates(Date previousDate, 
                   const std::vector<Date>& pDates,
                   Date::DayCountConvention dcc,
                   int iFixingDelay,
                   double dMargin,
                   double dMultiplier,
                   double dCap,
                   double dFloor,
                   const shared_ptr<YieldCurve>& pYieldcurve)
{
  
  std::vector<double> paymentRates;

  //Compute reference rates.
  std::vector<double> 
    pdReferenceRates = ComputeReferenceRates(previousDate, pDates, 
                                             dcc, iFixingDelay, pYieldcurve);

  //Compute coupons
  size_t nNbCoupons = pdReferenceRates.size();
  
  double dRate;

  for( size_t nIdxCoupon = 0; nIdxCoupon < nNbCoupons; ++nIdxCoupon)
  {
    dRate = pdReferenceRates[nIdxCoupon];

    dRate = dMultiplier * dRate + dMargin;
    
    if ( dRate < dFloor )
      dRate = dFloor;

    if ( dRate > dCap )
      dRate = dCap;

    //This rate obtained is annually, then:
    dRate *= Date::YearsDiff(previousDate, pDates[nIdxCoupon], dcc);
    
    paymentRates.push_back(dRate);

    previousDate = pDates[nIdxCoupon];
  }

  return paymentRates;

}

} // namespace finance

} // namespace ito33
