///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cbconversions.cpp
// Purpose:     convertible bond conversions class
// Author:      Laurence
// Created:     2004/03/15
// RCS-ID:      $Id: cbconversions.cpp,v 1.72 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/triggerincurrencyof.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"
#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/pricing/cashflows.h"
#include "ito33/pricing/cbconversions.h"
#include "ito33/pricing/cbcalls.h"

#include "ito33/pricing/cblikeparams.h"

namespace ito33
{

namespace pricing
{
   using namespace finance;


CBConversions::CBConversions
               ( const shared_ptr<ConversionSchedule>& pConversions,
                 Date valuationDate ) 
               : ConversionProvisions()
{
  if (!pConversions)
    return;

  const ConversionSchedule::Elements& conversions = pConversions->GetAll();
  if(conversions.empty())
    return;

  m_bKeepAccrued = pConversions->GetKeepAccrued();
  m_bForfeitCoupon = pConversions->GetForfeitCoupon();

  size_t nNbConversions = conversions.size();
  
  // Assume no windows will be broken up. If they are, the vectors
  // will need to resize themselves when push_back is called. 
  m_pdStartTimes.reserve(nNbConversions);
  m_pdEndTimes.reserve(nNbConversions);
  m_pdRatios.reserve(nNbConversions);
  m_pdTriggerRates.reserve(nNbConversions);
  m_pCoCoTypes.reserve(nNbConversions);
  m_pdCashs.reserve(nNbConversions);
  m_pbIsActive.reserve(nNbConversions);
  
  // Needed for coco
  bool bIsLastTriggerMet = pConversions->GetIsLastTriggerConditionMet();

  // Loop over the windows, copying data and breaking up windows as needed
  ConversionSchedule::Elements::const_iterator iter = conversions.begin();
  for ( ; iter != conversions.end(); ++iter)
  {
    const ConversionPeriod & period = *(*iter); 

    // Get the information about this window
    Date startDate      = period.GetStartDate();
    Date endDate        = period.GetEndDate();    
    double dTriggerRate = 0.0;
    double dExtremeRate = 0.0;
    double dRateChange  = 0.0;
    CoCoType cocoType   = CoCoType_Max;

    if ( period.HasCoCo() )
    {
      dTriggerRate = period.GetTrigger();
      dExtremeRate = period.GetExtremeTrigger();
      dRateChange  = period.GetChangeRate();
      cocoType     = period.GetCoCoType();
    }

    // if the window is too long, break it into quarters
    if ( IsWindowTooLong(startDate, endDate, cocoType, dTriggerRate, 
                         dExtremeRate, dRateChange, bIsLastTriggerMet) )
    {
      BreakupWindow(startDate, endDate, period.GetRatio(), period.GetCash(),
                    cocoType, dTriggerRate, dExtremeRate, dRateChange, 
                    bIsLastTriggerMet, valuationDate);
    }
    else
    {
      m_pdStartTimes.push_back( GetDoubleFrom(startDate) );
      m_pdEndTimes.push_back( GetDoubleFrom( endDate ) );    
      
      m_pdRatios.push_back( period.GetRatio() );
      m_pdCashs.push_back( period.GetCash() );
      
      m_pbIsActive.push_back( true );

      // Set as a normal period if the last trigger condition was met,
      // this period contains the valuation date, and it is a CoCo period
      // that allows conversion after the trigger is checked.
      if ( bIsLastTriggerMet == true &&
           valuationDate >= startDate && valuationDate < endDate && 
          (cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter || 
           cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate ||
           cocoType == CoCoType_CheckAnyTimeAndConvertAsOfCheckDate) )
      {
        dTriggerRate = 0.0;
        cocoType = CoCoType_Max;
      }

      // Set the coco features
      m_pCoCoTypes.push_back( cocoType );
      m_pdTriggerRates.push_back( dTriggerRate );

    } // if the period was too long and had to be broken

  } // loop over windows

  // If windows were broken up, the number of conversion windows
  // will be greater than the number of periods in the schedule
  m_nNbConversions = m_pdStartTimes.size();
}

double CBConversions::GetConversionPrice
       (TriggerAsPercentageOf triggerAsPercentageOf) const
{
  /*
  ASSERT_MSG(m_pParams->GetIndexConversion() != INVALIDINDEX,
       "GetConversionPrice must be called when index of conversion is valid");
  */

  // fall back to the first conversion ratio
  size_t nIdxConversion = m_pParams->GetIndexConversion();

  if (nIdxConversion == INVALIDINDEX)
    nIdxConversion = 0;

  double dConversionPrice;
  double dConversionRatio = GetRatio(nIdxConversion);

  if ( triggerAsPercentageOf == TriggerAsPercentageOf_IssuePrice )
    dConversionPrice = m_pParams->GetCBLike().GetIssuePrice()
                     / dConversionRatio;
  else if ( triggerAsPercentageOf == TriggerAsPercentageOf_Claim)
    dConversionPrice = m_pParams->GetClaim() / dConversionRatio;
  else
    dConversionPrice = (m_pParams->GetClaim() - m_pParams->GetAccruedInterest()) 
                     / dConversionRatio;
  
  if ( m_pParams->GetCBLike().IsCrossCurrency() )
  {
    if ( m_triggerInCurrencyOf == TriggerInCurrencyOf_Underlying )
      dConversionPrice /= m_dFixedFXRate;    
    else 
      dConversionPrice /= m_pParams->GetFXRate();  
  }

  return dConversionPrice; 
}

double CBConversions::GetConversionPrice(double dTime, 
  TriggerAsPercentageOf triggerAsPercentageOf, bool bPlus) const
{
  size_t nIdxConversion = m_pParams->GetConversions()
                        ->GetIndexContinousConversionAt(dTime, bPlus);

  // fall back to the first conversion ratio
  if (nIdxConversion == INVALIDINDEX)
    nIdxConversion = 0;

  double dConversionRatio = GetRatio(nIdxConversion);
  double dConversionPrice;

  if ( triggerAsPercentageOf == TriggerAsPercentageOf_IssuePrice )
    dConversionPrice = m_pParams->GetCBLike().GetIssuePrice()
                     / dConversionRatio;  
  else if ( triggerAsPercentageOf == TriggerAsPercentageOf_Claim)   
     dConversionPrice = 
        m_pParams->GetCBLike().GetClaim(dTime, bPlus) / dConversionRatio;  
  else
    dConversionPrice = (  m_pParams->GetCBLike().GetClaim(dTime, bPlus)
                         - m_pParams->GetCashFlows()
                                    ->GetAccruedInterest(dTime, bPlus)) 
                     / dConversionRatio; 

  if ( m_pParams->GetCBLike().IsCrossCurrency() )
  {
    if ( m_triggerInCurrencyOf == TriggerInCurrencyOf_Underlying )
      dConversionPrice /= m_dFixedFXRate;    
    else 
      dConversionPrice /= m_pParams->GetFXRate(dTime);  
  }
  
  return dConversionPrice;
}

bool CBConversions::GetGrossParities
     (const double* /* pdS */, size_t nNbS, const double* pdNewSharePrices, 
      double* pdValues) const
{
  if ( m_pParams->GetIndexConversion() == INVALIDINDEX )
    return false;
  
  double dRatio = m_pdRatios[m_pParams->GetIndexConversion()];
  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    dRatio *= m_pParams->GetFXRate();
  
  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] = dRatio * pdNewSharePrices[nIdxS];  
    
  return true;
}

void CBConversions::ComputeRoots(double dTime,
  size_t &nNbRoots, numeric::mesh::Root *pRoots, bool bPlus )
{
  // maybe roots should be calculated in Calls? (Zhang 08/17/2004)  
  CBCalls* pCalls = (CBCalls*)m_pParams->GetCalls();

  // note: the index is setup before calling this function
  size_t
    nIdxCall = m_pParams->GetIndexCall(),
    nIdxConversion = m_pParams->GetIndexConversion();

  nNbRoots = 1;
   
  TriggerAsPercentageOf of = pCalls->GetTriggerAsPercentageOf();

  double dTriggerLevel = pCalls->GetTriggerRate(nIdxCall)
                       * GetConversionPrice(dTime, of, bPlus);  

  // in case of a virtual call, call constraint is off, but we still want to
  // add call trigger as root.
  if ( !pCalls->IsActive() && dTriggerLevel > 0.) 
  {
    pRoots[0].m_dRoot     = dTriggerLevel;
    pRoots[0].m_iRootType = numeric::mesh::RootType_Real;

    return;
  }

  double dConvAccrued;

  if (m_bKeepAccrued)
    dConvAccrued = m_pParams->GetCashFlows()->GetAccruedInterest(dTime, bPlus);
  else
    dConvAccrued = 0.;
  
  double dRatio = m_pdRatios[nIdxConversion];

  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    dRatio *= m_pParams->GetFXRate(dTime);
  
  pRoots[0].m_dRoot = (   pCalls->GetStrikeWithoutCoupon(dTime, bPlus)  
                        + m_pdCashs[nIdxConversion] - dConvAccrued ) / dRatio;
  
  if (pRoots[0].m_dRoot < dTriggerLevel)
    pRoots[0].m_dRoot = dTriggerLevel;

  pRoots[0].m_iRootType = numeric::mesh::RootType_Real;
}

bool CBConversions::IsWindowTooLong(Date startDate,
                                    Date endDate,
                                    CoCoType cocoType,
                                    double dTriggerRate,
                                    double dExtremeRate,
                                    double dRateChange,
                                    bool bIsLastTriggerMet)
{

  if (cocoType == CoCoType_Max)
    return false;

  // If the last trigger condition was met, "convert as of check date"
  // coco periods allow conversion. No need to divide.
  if ( bIsLastTriggerMet == true &&
       (cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate ||
        cocoType == CoCoType_CheckAnyTimeAndConvertAsOfCheckDate) )
    return false;

  // Assume that any quarterly coco type needs to be broken.  For non-quarterly
  // periods, break up if the length is greater than one quarter and the trigger
  // rate changes
  Date tmpEndDate = startDate;
  tmpEndDate.AddMonths(3);
  if ( cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate
    || cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter
    || (   dRateChange != 0.0 && dExtremeRate != dTriggerRate 
        && tmpEndDate < endDate) )
  {
    return true;
  }
  else
  {
    return false;
  }
}

void CBConversions::BreakupWindow(Date windowStartDate,
                                  Date finalEndDate,
                                  double dRatio,
                                  double dCash,
                                  CoCoType cocoType,
                                  double dTriggerRate,
                                  double dExtremeRate,
                                  double dRateChange,
                                  bool bIsLastTriggerMet,                                   
                                  Date valuationDate)
{

  // Break into quarterly periods (3 months)
  Date windowEndDate  = windowStartDate;
  windowEndDate.AddMonths(3);

  // Make sure we are not extending the period (this function
  // is always called for certain coco windows, regardless of the
  // window length)
  if (windowEndDate > finalEndDate)
    windowEndDate = finalEndDate;

  for ( ; ; )
  {
    // Store the current data
    m_pdStartTimes.push_back( GetDoubleFrom(windowStartDate) );
    m_pdEndTimes.push_back( GetDoubleFrom(windowEndDate) );    
    m_pdRatios.push_back( dRatio );
    m_pdCashs.push_back( dCash );
    m_pbIsActive.push_back( true );

    if ( bIsLastTriggerMet == true &&
         cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter &&
         valuationDate >= windowStartDate &&
         valuationDate < windowEndDate )
    {
      m_pdTriggerRates.push_back( 0.0 );
      m_pCoCoTypes.push_back( CoCoType_Max );
    }
    else
    {
      m_pdTriggerRates.push_back( dTriggerRate );
      m_pCoCoTypes.push_back( cocoType );
    }

    // Check if this was the last window
    if (windowEndDate >= finalEndDate)
      break;

    // Update for the next window
    windowStartDate = windowEndDate;
    windowEndDate.AddMonths(3);
    if (windowEndDate > finalEndDate)
      windowEndDate = finalEndDate;

    // The rate change is an annual rate, and we are compounding
    // 4 times a year (since periods are one quarter)
    dTriggerRate *= (1.0 + dRateChange / 4.0);
    if (dRateChange > 0.0 && dTriggerRate > dExtremeRate)
      dTriggerRate = dExtremeRate;
    else if (dRateChange < 0.0 && dTriggerRate < dExtremeRate)
      dTriggerRate = dExtremeRate;

  } // for loop creating new small windows
}


} // namesapce pricing

} // namespace ito33
