///////////////////////////////////////////////////////////////////////////
// Name:        pricing/conversionprovisions.cpp
// Purpose:     Abstract base class for conversion provisions
// Author:      Laurence
// Created:     2004/03/12
// RCS-ID:      $Id: conversionprovisions.cpp,v 1.17 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////



#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/constants.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/conversionprovisions.h"

namespace ito33
{
  //Implementation of the AutoPtrDeleter for ConversionProvisions class.
  ITO33_IMPLEMENT_AUTOPTR(pricing::ConversionProvisions);
}

namespace ito33
{

  using namespace finance;

namespace pricing
{


// begin <= time < end (+)  begin < time <= end(-)
size_t ConversionProvisions::GetIndexContinousConversionAt
       (double dTime, bool bPlus)
{
  size_t
    nIdx;

  if ( bPlus )
    for(nIdx = 0; nIdx < m_nNbConversions; nIdx++)
    {
      if( GetStartTime(nIdx) != GetEndTime(nIdx) 
          && !ito33::numeric::IsBefore(dTime, GetStartTime(nIdx)) 
          && ito33::numeric::IsBefore(dTime, GetEndTime(nIdx)) )
        break;
    }
  else
    for(nIdx = 0; nIdx < m_nNbConversions; nIdx++)
    {
      if( GetStartTime(nIdx) != GetEndTime(nIdx) 
          && ito33::numeric::IsAfter(dTime, GetStartTime(nIdx)) 
          && !ito33::numeric::IsAfter(dTime, GetEndTime(nIdx)) )
        break;
    }
  
  if(nIdx >= m_nNbConversions)
    nIdx = INVALIDINDEX;

  return nIdx;
}

double ConversionProvisions::GetTriggerLevel() const
{
  size_t nIdxConversion = m_pParams->GetIndexConversion();

  ASSERT_MSG(nIdxConversion < INVALIDINDEX, "get trigger but not conversion");

  return m_pdTriggerRates[nIdxConversion] * m_pParams->GetConversionPrice();
}

/*
  Implementation can be optimized if we compute only the gross parity for spots
  above trigger
*/
bool ConversionProvisions::GetConversionConstraintValues
     (const double* pdS, size_t nNbS, 
      const double* pdNewSharePrices, double* pdValues) const
{
  size_t nIdxConversion = m_pParams->GetIndexConversion();
  if ( nIdxConversion == INVALIDINDEX )
    return false;

  if ( m_pbIsActive[nIdxConversion] == false )
    return false;
  
  GetGrossParities(pdS, nNbS, pdNewSharePrices, pdValues);      
  
  double dTriggerLevel = GetTriggerLevel();
  
  size_t nIdxS;
  for (nIdxS = 0; 
       nIdxS < nNbS && numeric::LessThanTrigger(pdS[nIdxS], dTriggerLevel); 
       nIdxS++)
    pdValues[nIdxS] = 0;
 
  // Add accured and coupon for spot above trigger
  double dTmp = - m_pdCashs[m_pParams->GetIndexConversion()];
 
  if (m_bKeepAccrued)
    dTmp += m_pParams->GetAccruedInterest();

  if (!m_bForfeitCoupon)
    dTmp += m_pParams->GetCouponAmount();  
  
  for ( ; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] += dTmp;

  return true;
}

bool ConversionProvisions::GetForcedConversionValues
     (const double* pdS, size_t nNbS, const double* pdNewSharePrices,
      double* pdValues) const
{
  if ( m_pParams->GetIndexConversion() == INVALIDINDEX )
    return false;

  double dTmp = - m_pdCashs[m_pParams->GetIndexConversion()];
  
  if (m_bKeepAccrued) 
    dTmp += m_pParams->GetAccruedInterest();
  
  if (!m_bForfeitCoupon)
    dTmp += m_pParams->GetCouponAmount();

  GetGrossParities(pdS, nNbS, pdNewSharePrices, pdValues);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] += dTmp;

  return true;
}

void ConversionProvisions::RemoveTriggers()
{
  for (size_t nIdx = 0; nIdx < m_nNbConversions; nIdx++)
  {
    m_pCoCoTypes[nIdx] = finance::CoCoType_Max;
    m_pdTriggerRates[nIdx] = 0.;
  }
}

void ConversionProvisions::DeactivateDiscretePeriods()
{
  // Only deactivate quarterly coco type periods.  All others
  // should be included in the base path of path dependent pricing.
  // This function should only be called for path dependent pricing.
  for (size_t nIdx = 0; nIdx < m_nNbConversions; nIdx++)
  {
    if ( m_pCoCoTypes[nIdx]==CoCoType_CheckQuarterlyAndConvertDuringNextQuarter
      || m_pCoCoTypes[nIdx]==CoCoType_CheckQuarterlyAndConvertAsOfCheckDate)
    {
      m_pbIsActive[nIdx] = false;
    } //end if    
  } //end loop
} // ConversionProvisions::DeactivateDiscretePeriods


bool ConversionProvisions::HasPathDepCoCo(double dValuationTime, 
                                          double dMaturityTime)
{
  // Scan the periods for a coco clause
  for (size_t nIdx = 0; nIdx < m_nNbConversions; nIdx++)
    if ( m_pbIsActive[nIdx] )
    {
      CoCoType cocoType = GetPeriodType(nIdx);
      double dStartTime = GetStartTime(nIdx);
      double dEndTime   = GetEndTime(nIdx);
     
      // Make sure this conversion window is within the pricing dates
      if ( numeric::IsAfter(dStartTime, dMaturityTime)
        || numeric::IsBefore(dEndTime, dValuationTime) )
        continue;

      if (   cocoType != CoCoType_Max 
          && cocoType != CoCoType_CheckAnyTimeAndConvertOnCheckDate)
      {
        return true;
      }
    } // loop over conversion periods

  return false;
}


} // namespace pricing

} // namespace ito33
