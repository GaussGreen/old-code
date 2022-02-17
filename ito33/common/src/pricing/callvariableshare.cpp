///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/callvariableshare.cpp
// Purpose:     variable share call for PEPS-like
// Author:      Wang
// Created:     2004/08/23
// RCS-ID:      $Id: callvariableshare.cpp,v 1.3 2005/05/24 10:33:33 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/bondlike/call.h"

#include "ito33/finance/bondlike/generalizedpepslike.h"

#include "ito33/pricing/callvariableshare.h"
#include "ito33/pricing/mandatory_payoff_structure.h"

#include "ito33/pricing/cblikeparams.h"


namespace ito33
{

namespace pricing
{

CallVariableShare::CallVariableShare
      (const finance::GeneralizedPEPSLike& peps)
  : CallProvisions()
{
  ASSERT_MSG
      (
        peps.GetGeneralizedPEPSLikeCall() &&
        peps.GetGeneralizedPEPSLikeCall()->GetType()
          == finance::GeneralizedPEPSLikeCallType_VariableShare,
        "Inconsistent data."
      );

  const finance::GeneralizedPEPSLikeCall*
      pCall = peps.GetGeneralizedPEPSLikeCall().get();

  GetCommonCallData(*pCall);
  
  double dShiftTimeValue = 0.0;

  if ( m_bHasNoticePeriod )
  {
    dShiftTimeValue = - m_dNoticePeriod;
  }

  Date startDate = pCall->GetStartDate();
  Date maturityDate = peps.GetMaturityDate();

  Date endDate = pCall->GetEndDate();
  endDate = (endDate > maturityDate) ? maturityDate : endDate;

  if ( startDate > endDate)
    return;

  m_nNbCalls = 1;

  m_pdStartTimes.resize(m_nNbCalls);
  m_pdEndTimes.resize(m_nNbCalls);
  m_pData.resize(m_nNbCalls);
  m_pdTriggerRates.resize(m_nNbCalls);
  m_pbIsActives.resize(m_nNbCalls);
  m_pdYieldToCall.resize(m_nNbCalls);

  m_pdStartTimes[0] = GetDoubleFrom(startDate) + dShiftTimeValue;
  m_pdEndTimes[0] = GetDoubleFrom(endDate) + dShiftTimeValue;
  m_pData[0] = MandatoryPayoffStructure
                    (
                      peps.GetDownsideConversionRatio(),
                      peps.GetUpsideBaseConversionRatio(),
                      peps.GetLowerStrike(),
                      peps.GetHigherStrike()
                    );
  m_pdTriggerRates[0] = pCall->GetTrigger();
  m_pbIsActives[0] = true;
  m_pdYieldToCall[0] = 0.0;
}

void CallVariableShare::GetCallStrikes
     (const double* pdS, size_t nNbS, const double* pdNewSharePrices, 
      double* pdValues) const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT(nIdxCall != INVALIDINDEX);
 
  double dTmp = 0.;

  if (m_bKeepAccrued)
    dTmp += m_pParams->GetAccruedInterest();

  if (!m_bForfeitCoupon)
    dTmp += m_pParams->GetCouponAmount();

  if(m_pParams->GetCBLike().IsCrossCurrency())
    m_pData[nIdxCall].GetValues
      (pdS, pdValues, nNbS, pdNewSharePrices, m_pParams->GetFXRate());
  else
    m_pData[nIdxCall].GetValues(pdS, pdValues, nNbS, pdNewSharePrices);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] += dTmp;
}

void CallVariableShare::GetCallStrikesWithoutCoupon
     (double dTime, const double* pdS, size_t nNbS, 
      const double* pdNewSharePrices, double* pdValues, 
      bool bPlus) const
{
  size_t nIdxCall = m_pParams->GetIndexCall();

  ASSERT(nIdxCall != INVALIDINDEX);
 
  double dTmp = (m_bKeepAccrued)
              ? m_pParams->GetCashFlows()->GetAccruedInterest(dTime, bPlus)
              : 0;

  if(m_pParams->GetCBLike().IsCrossCurrency())
    m_pData[nIdxCall].GetValues
      (pdS, pdValues, nNbS, pdNewSharePrices, m_pParams->GetFXRate());
  else
    m_pData[nIdxCall].GetValues(pdS, pdValues, nNbS, pdNewSharePrices);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdValues[nIdxS] += dTmp;
}

void CallVariableShare::Clear()
{
  *this = CallVariableShare();
}

} // namesapce pricing

} // namespace ito33
