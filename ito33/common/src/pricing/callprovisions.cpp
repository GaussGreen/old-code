///////////////////////////////////////////////////////////////////////////
// Name:        pricing/callprovisions.cpp
// Purpose:     Abstract base class for call provisions
// Author:      Nabil
// Created:     2004/03/11
// RCS-ID:      $Id: callprovisions.cpp,v 1.36 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/bondlike/call.h"

#include "ito33/pricing/callprovisions.h"
#include "ito33/pricing/cblikeparams.h"

namespace ito33
{
  //Implementation of the AutoPtrDeleter for CallProvisions class.
  ITO33_IMPLEMENT_AUTOPTR(pricing::CallProvisions);
}

namespace ito33
{

using numeric::AreTimesEqual;
  
namespace pricing
{


void CallProvisions::GetCommonCallData(const finance::Call& call)
{
  m_bKeepAccrued     = call.GetKeepAccrued();
  m_bForfeitCoupon   = call.GetForfeitCoupon();
  m_bHasNoticePeriod = call.HasNoticePeriod();
  m_dNoticePeriod    = call.GetNoticePeriod() * ONEDAY;
  m_nTriggerPeriod   = call.GetTriggerPeriod();
  m_nTriggerHistory  = call.GetTriggerHistory();
  m_triggerAsPercentageOf = call.GetTriggerAsPercentageOf();

  // Cap the trigger history in the pricing level. Otherwise, 
  // the path dependent code tries to return an invalid path
  if (m_nTriggerHistory > m_nTriggerPeriod)
    m_nTriggerHistory = m_nTriggerPeriod;
  
  if (call.HasMakeWhole()) // no make-whole
  {
     m_iMakeWholeType = call.GetMakeWholeType();
     
     if(m_iMakeWholeType == finance::MakeWholeType_Premium)
       m_dMakeWholePremium = call.GetMakeWholePremium();
     else
       m_bDiscountCouponForMakeWhole = call.IsPVCouponMakeWhole();
  }
}

// begin <= time < end (+)  begin < time <= end(-)
size_t CallProvisions::GetIndexContinousCallAt(double dTime, bool bPlus)
{
  size_t nIdx;

  if ( bPlus )
    for (nIdx = 0; nIdx < m_nNbCalls; nIdx++)
    {
      if( GetStartTime(nIdx) != GetEndTime(nIdx) 
          && !ito33::numeric::IsBefore(dTime, GetStartTime(nIdx)) 
          && ito33::numeric::IsBefore(dTime, GetEndTime(nIdx)) )
        break;
    }
  else
    for(nIdx = 0; nIdx < m_nNbCalls; nIdx++)
    {
      if( GetStartTime(nIdx) != GetEndTime(nIdx) 
          && ito33::numeric::IsAfter(dTime, GetStartTime(nIdx)) 
          && !ito33::numeric::IsAfter(dTime, GetEndTime(nIdx)) )
        break;
    }

  if ( nIdx >= m_nNbCalls )
    nIdx = INVALIDINDEX;

  return nIdx;
}

double CallProvisions::GetTriggerLevel() const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT_MSG(nIdxCall < INVALIDINDEX, "get trigger but not call");

  return m_pdTriggerRates[nIdxCall] * m_pParams->GetConversionPriceForCall();
}

bool CallProvisions::GetCallConstraintValues
      ( const double* pdS, size_t nNbS, double* pdValues,
        size_t& nIdxStartConversion, const double* pdNewSharePrices,
        bool bNoticed) const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT_MSG(m_pbIsActives[nIdxCall], 
             "Calling GetCallConstraintValues for desactivated calls");

  if ( nIdxCall == INVALIDINDEX )
    return false;
   
  nIdxStartConversion = INVALIDINDEX;

  /*
    If there is no call notice, compute internally the call strikes.
    Otherwise, use the external values computed by solving the equation on 
    call notice.
  */
  if ( !bNoticed )
    GetCallStrikes(pdS, nNbS, pdNewSharePrices, pdValues); 

  size_t nIdxS;

  double dTriggerLevel = GetTriggerLevel();
  
  // No max constraints for spot smaller than the trigger
  for (nIdxS = 0;
       nIdxS < nNbS  && numeric::LessThanTrigger(pdS[nIdxS], dTriggerLevel);
       nIdxS++)
  {
    pdValues[nIdxS] = 1.e12;
  }
  
  // record the index where the issuer can call back 
  size_t nIdxSAboveTrigger = nIdxS;

  /*
    Take the max of call strikes and forced conversion values (if any).
    The holder of the security can convert regardless the conversion trigger
    when the issuer calls back.
  */
  Array<double> pdForcedConversionValues(nNbS);

  bool 
    bHasConversion = m_pParams->GetConversions()->GetForcedConversionValues
                     ( pdS, nNbS, pdNewSharePrices,
                       pdForcedConversionValues.Get() );

  if (bHasConversion)
  {
    for (nIdxS = nIdxSAboveTrigger;
         nIdxS < nNbS && pdValues[nIdxS] > pdForcedConversionValues[nIdxS];
         nIdxS++)
     ;
    nIdxStartConversion = nIdxS; // this is output
    for (; nIdxS < nNbS; nIdxS++)
    {
      // after this point conversion is forced
      pdValues[nIdxS] = pdForcedConversionValues[nIdxS];  
    }
  }
  
  /*
    Add Make whole premium for soft call. The holder of the security gets 
    the make whole premium whenever the issuer makes a soft call.
  */
  if (dTriggerLevel > 0 && HasMakeWhole() )
  {
    double dPremium = m_pParams->GetCurrentMakeWholePositivePremium();
    for (nIdxS = nIdxSAboveTrigger; nIdxS < nNbS; nIdxS++)
      pdValues[nIdxS] += dPremium;
  }
  
  return true;
}

void CallProvisions::DeactivateSoftCall()
{
  size_t nIdxCall;

  for ( nIdxCall = 0; nIdxCall <  m_nNbCalls; nIdxCall++)
  {   
    double dStartTime   = GetStartTime(nIdxCall);
    double dEndTime     = GetEndTime(nIdxCall);
    double dTriggerRate = GetTriggerRate(nIdxCall);

    if ( !AreTimesEqual(dStartTime, dEndTime) && dTriggerRate > 0 )
    {
      m_pbIsActives[nIdxCall] = false;
    }
  } 
}

bool CallProvisions::IsActive() const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  if (nIdxCall == INVALIDINDEX)
    return false;

  return m_pbIsActives[nIdxCall];
} // CallProvisions::IsActive() const


void CallProvisions::DeactivateNoticePeriod()
{
  m_bHasNoticePeriod = false;
  m_dNoticePeriod = 0.0;
} //CallProvisions::DeactivateNoticePeriod()

void CallProvisions::DeactivateTriggerPeriod()
{
 m_nTriggerPeriod = 0;
} // CallProvisions::DeactivateTriggerPeriod()

void CallProvisions::DeactivateAllCalls()
{
  DeactivateNoticePeriod();
  DeactivateTriggerPeriod();

  for (size_t nIdxCall = 0; nIdxCall <  m_nNbCalls; nIdxCall++)
    m_pbIsActives[nIdxCall] = false;
 
}//CallProvisions::DeactivateAllCalls()

bool CallProvisions::HasPathDepCall(double dValuationTime, 
                                    double dStoppingTime) const
{
  if ( m_nTriggerPeriod > 0 )
  {
    /*
      Loop over the periods to make sure at least one of them
      applies during the pricing period
    */
    for (size_t nIdx = 0 ; nIdx < m_nNbCalls ; nIdx++)
      if (m_pbIsActives[nIdx])
      {
        double dStartTime = m_pdStartTimes[nIdx];
        double dEndTime   = m_pdEndTimes[nIdx];

        // Make sure this conversion window overlaps the pricing dates.
        // Use IsEqualOrXXX to ensure overlap, and not just matching
        // at valuation or stopping time
        if (   numeric::IsEqualOrAfter(dStartTime, dStoppingTime)
            || numeric::IsEqualOrBefore(dEndTime, dValuationTime) )
          continue;

        // filter out hard call periods and monodates
        if (   m_pdTriggerRates[nIdx] > 0 
            && !numeric::AreTimesEqual(dStartTime, dEndTime))
        {
          return true;
        }
      }
  }

  return false;
}


} // namespace pricing

} // namespace ito33
