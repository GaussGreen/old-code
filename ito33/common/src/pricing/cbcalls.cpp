///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cbcalls.cpp
// Purpose:     convertible bond calls class
// Author:      Nabil
// Created:     2004/03/11
// RCS-ID:      $Id: cbcalls.cpp,v 1.75 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/dateutils.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/pricing/cashflows.h"
#include "ito33/pricing/cbcalls.h"

#include "ito33/pricing/cblikeparams.h"

namespace ito33
{

using namespace finance;

namespace pricing
{
  
CBCalls::CBCalls(const shared_ptr<finance::CallSchedule>& pCall) 
               : CallProvisions()
{
  if ( !pCall )
  {
    return;
  }

  const CallSchedule::Elements& calls = pCall->GetAll();

  if ( calls.empty() ) 
  {
    return;
  }

  GetCommonCallData(*pCall);
  
  size_t nNbCalls = calls.size();
  
  m_pdStartTimes.resize( nNbCalls);
  m_pdEndTimes.resize(nNbCalls);
  m_pdStrikeRates.resize(nNbCalls);
  m_pdTriggerRates.resize(nNbCalls);
  m_pbIsActives.resize(nNbCalls);
  m_pdYieldToCall.resize(nNbCalls);

  double dShiftTimeValue = 0.0;

  if ( m_bHasNoticePeriod )
  {
    dShiftTimeValue = - m_dNoticePeriod;
  }
  
  CallSchedule::Elements::const_iterator iter = calls.begin();
  size_t nIdxCall = 0;

  for ( ; iter != calls.end(); ++iter, ++nIdxCall)
  {
    const CallPeriod& callPeriod = *(*iter);

    m_pdStartTimes[nIdxCall] = GetDoubleFrom( callPeriod.GetStartDate() )
                             + dShiftTimeValue;
    
    m_pdEndTimes[nIdxCall] = GetDoubleFrom( callPeriod.GetEndDate() )
                           + dShiftTimeValue;
    
    if ( callPeriod.HasYield() )
    {
      m_pdYieldToCall[nIdxCall] = callPeriod.GetGuaranteedYield();
      m_pdStrikeRates[nIdxCall] = 0;
    }
    else
    {
      m_pdStrikeRates[nIdxCall] = callPeriod.GetStrike();
      m_pdYieldToCall[nIdxCall] = 0;
    }
   
    m_pdTriggerRates[nIdxCall] = callPeriod.GetTrigger();


    //by default all the calls are active
    m_pbIsActives[nIdxCall] = true;
  }

  m_nNbCalls = nIdxCall;
}

double CBCalls::GetStrikeWithoutCoupon(double dTime, bool bPlus) const
{

  /*
    Note that although this function has current time as parameter, we still 
    shouldn't use GetIndexContinousCallAt(dTime) to compute the index since 
    the call period was shifted. 

    The means that we should use this function carefully. The index of the call 
    must be setup before. 
  */
  size_t nIdxCall     = m_pParams->GetIndexCall();

  // Check if strike rate or yield to call
  double dStrike = 0.0;
  if ( m_pdStrikeRates[nIdxCall] <= 0.0 )
    dStrike = 
      m_pParams->GetClaimFromYield(dTime, m_pdYieldToCall[nIdxCall], bPlus);
  else
  {
    // Same as GetStrike below, but no coupon
    double dAccruedAmount =
      m_pParams->GetCashFlows()->GetAccruedInterest(dTime, bPlus);
 
    double dClaim  = m_pParams->GetCBLike().GetClaim(dTime, bPlus);

    double dBase = dClaim - dAccruedAmount;
    
    dStrike = m_pdStrikeRates[nIdxCall] * dBase;
  
    if (m_bKeepAccrued)
      dStrike += dAccruedAmount;
  }
  
  return dStrike;
}

void CBCalls::GetCallStrikesWithoutCoupon
     (double dTime, const double* /* pdS */, size_t nNbS, 
      const double* /* pdNewSharePrices */, double* pdValues, 
      bool bPlus) const
{
  double dStrike = GetStrikeWithoutCoupon(dTime, bPlus);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] = dStrike;
}


double CBCalls::GetStrike() const
{
  size_t nIdxCall     = m_pParams->GetIndexCall();

  // Check if strike rate or yield to call
  double dStrike = 0.0;
  if (  m_pdStrikeRates[nIdxCall] <= 0.0 )
  {
    double dClaim = m_pParams->GetClaimFromYield(m_pdYieldToCall[ nIdxCall ]);
    dStrike = dClaim + m_pParams->GetCouponAmount();
  }
  else
  {
    // strike is applied to principal = claim minus accrued
    double dClaim = m_pParams->GetClaim();
    double dAccruedAmount = m_pParams->GetAccruedInterest();
    double dBase = dClaim - dAccruedAmount;

    dStrike = m_pdStrikeRates[nIdxCall] * dBase;
  
    // accrued interest and forfeit coupon only apply for strikes
    if (m_bKeepAccrued)
      dStrike += dAccruedAmount;

    if ( !m_bForfeitCoupon )
      dStrike += m_pParams->GetCouponAmount();
  }

  return dStrike;
}

void CBCalls::GetCallStrikes(const double* /* pdS */, size_t nNbS, 
                             const double* /* pdNewSharePrices */, 
                             double* pdValues) const
{
  ASSERT_MSG(m_pParams->GetIndexCall() < INVALIDINDEX,
             "Invalid index when computing call strikes");

  // Use the call value as the max constraint
  double dCallStrike = GetStrike();

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] = dCallStrike;
}

void CBCalls::Clear()
{
  CBCalls calls;

  *this = calls;
}


} // namespace pricing

} // namespace ito33
