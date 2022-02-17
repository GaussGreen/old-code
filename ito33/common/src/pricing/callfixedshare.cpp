///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/callfixedshare.cpp
// Purpose:     fixed share call for PEPS-like
// Author:      Wang
// Created:     2004/08/23
// RCS-ID:      $Id: callfixedshare.cpp,v 1.22 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/callfixedshare.h"

#include "ito33/pricing/callfixedshare.h"

#include "ito33/pricing/cblikeparams.h"


namespace ito33
{

namespace pricing
{

CallFixedShare::CallFixedShare
      (const finance::GeneralizedPEPSLike& peps)
  : CallProvisions()
{
  ASSERT_MSG
      (
        peps.GetGeneralizedPEPSLikeCall() &&
        peps.GetGeneralizedPEPSLikeCall()->GetType()
          == finance::GeneralizedPEPSLikeCallType_FixedShare,
        "Inconsistent data."
      );

  const finance::GeneralizedPEPSLikeCall*
      pCall = peps.GetGeneralizedPEPSLikeCall().get();

  GetCommonCallData(*pCall);
  
  double dShiftTimeValue = 0.0;

  Date startDate = pCall->GetStartDate();
  Date maturityDate = peps.GetMaturityDate();

  Date endDate = pCall->GetEndDate();
  endDate = (endDate > maturityDate) ? maturityDate : endDate;

  if ( startDate > endDate)
    return;

  m_nNbCalls = 1;

  m_pdStartTimes.resize(m_nNbCalls);
  m_pdEndTimes.resize(m_nNbCalls);
  m_pdRatios.resize(m_nNbCalls);
  m_pdTriggerRates.resize(m_nNbCalls);
  m_pbIsActives.resize(m_nNbCalls);
  m_pdYieldToCall.resize(m_nNbCalls);
  
  if ( m_bHasNoticePeriod )
  {
    dShiftTimeValue = - m_dNoticePeriod;
  }

  m_pdStartTimes[0] = GetDoubleFrom(startDate) + dShiftTimeValue;
  m_pdEndTimes[0] = GetDoubleFrom(endDate) + dShiftTimeValue;
  m_pdRatios[0] = peps.GetMinConversionRatio();
  m_pdTriggerRates[0] = pCall->GetTrigger();
  m_pbIsActives[0] = true;
  m_pdYieldToCall[0] = 0.0;
}


CallFixedShare::CallFixedShare
(const shared_ptr<finance::CallFixedShare>& pCall, Date maturityDate) 
  : CallProvisions()
{
  if ( !pCall )
  {
    return;
  }

  GetCommonCallData(*pCall);
  
  double dShiftTimeValue = 0.0;

  if ( m_bHasNoticePeriod )
  {
    dShiftTimeValue = - m_dNoticePeriod;
  }

  Date startDate = pCall->GetStartDate();

  Date endDate = pCall->GetEndDate();
  endDate = (endDate > maturityDate) ? maturityDate : endDate;

  if ( startDate > endDate)
    return;

  m_nNbCalls = 1;

  m_pdStartTimes.resize(m_nNbCalls);
  m_pdEndTimes.resize(m_nNbCalls);
  m_pdRatios.resize(m_nNbCalls);
  m_pdTriggerRates.resize(m_nNbCalls);
  m_pbIsActives.resize(m_nNbCalls);
  m_pdYieldToCall.resize(m_nNbCalls);

  m_pdStartTimes[0] = GetDoubleFrom(startDate) + dShiftTimeValue;
  m_pdEndTimes[0] = GetDoubleFrom(endDate) + dShiftTimeValue;
  m_pdRatios[0] = pCall->GetRatio();
  m_pdTriggerRates[0] = pCall->GetTrigger();
  m_pbIsActives[0] = true;
  m_pdYieldToCall[0] = 0.0;
}

void CallFixedShare::GetCallStrikes
                     (const double* /* pdS */, size_t nNbS, 
                      const double* pdNewSharePrices, double* pdValues) const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT(nIdxCall != INVALIDINDEX);
 
  double dTmp = 0.;

  if (m_bKeepAccrued)
    dTmp += m_pParams->GetAccruedInterest();

  if (!m_bForfeitCoupon)
    dTmp += m_pParams->GetCouponAmount();

  double dRatio = m_pdRatios[nIdxCall];
  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    dRatio *= m_pParams->GetFXRate();

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] = pdNewSharePrices[nIdxS] * dRatio + dTmp;
}

void CallFixedShare::GetCallStrikesWithoutCoupon
     (double dTime, const double* /* pdS */, size_t nNbS, 
      const double* pdNewSharePrices, double* pdValues, 
      bool bPlus) const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT(nIdxCall != INVALIDINDEX);
 
  double dTmp = (m_bKeepAccrued)
              ? m_pParams->GetCashFlows()->GetAccruedInterest(dTime, bPlus)
              : 0;

  double dRatio = m_pdRatios[nIdxCall];
  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    dRatio *= m_pParams->GetFXRate(dTime);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] = pdNewSharePrices[nIdxS] * dRatio + dTmp;
}

void CallFixedShare::Clear()
{
  *this = CallFixedShare();
}

} // namesapce pricing

} // namespace ito33
