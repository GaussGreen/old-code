///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cbputs.cpp
// Purpose:     convertible bond puts class
// Author:      Laurence
// Created:     2004/03/11
// RCS-ID:      $Id: cbputs.cpp,v 1.28 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/dateutils.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/cbputs.h"
#include "ito33/pricing/cblikeparams.h"

#include "ito33/finance/bondlike/putschedule.h"

namespace ito33
{

//Implementation of the AutoPtrDeleter for CBPuts class.
ITO33_IMPLEMENT_AUTOPTR(pricing::CBPuts);

namespace pricing
{

 
CBPuts::CBPuts(const shared_ptr<finance::PutSchedule>& pPuts)
{
  DefaultInit();

  if ( !pPuts )
    return;

  const finance::PutSchedule::Elements& puts = pPuts->GetAll();

  if ( puts.empty() )
    return;

  m_nNbPuts = puts.size();

  m_pdTimes.resize(m_nNbPuts);
  m_pdStrikeRates.resize(m_nNbPuts);
  m_pdYieldToPut.resize(m_nNbPuts);

  // If any put has a yield, then KeepAccrued must be true
  // and ForfeitCoupon must be false.
  bool bHasYield = false;

  finance::PutSchedule::Elements::const_iterator iter = puts.begin();      
  for(size_t i = 0; iter != puts.end(); ++iter, ++i)
  {
    m_pdTimes[i] = GetDoubleFrom(iter->first);
    if ( iter->second.bYield )
      bHasYield = true;
    m_pdStrikeRates[i] = iter->second.dStrike;
    m_pdYieldToPut[i]  = iter->second.dYield;
    // strike < 0 <==> put has a yield
  }

  m_bKeepAccrued = pPuts->GetKeepAccrued();
  m_bForfeitCoupon = pPuts->GetForfeitCoupon();
}

bool CBPuts::GetPutConstraintValues
             (const double* /* pdS */, size_t nNbS, double* pdValues) const
{
  size_t nIdxPut     = m_pParams->GetIndexPut();
  
  if (nIdxPut == INVALIDINDEX)
    return false;
   
  // Check if strike rate or yield to put
  double dStrike = 0.0;
  
  if ( m_pdStrikeRates[nIdxPut] <= 0) // has yield
  {
    double dClaim = m_pParams->GetClaimFromYield(m_pdYieldToPut[nIdxPut]);
    dStrike = dClaim + m_pParams->GetCouponAmount();
  }
  else
  {
    // strike is applied to principal = claim minus accrued
    double dClaim = m_pParams->GetClaim();
    double dAccruedAmount = m_pParams->GetAccruedInterest();  
    double dBase = dClaim - dAccruedAmount;
    dStrike = m_pdStrikeRates[nIdxPut] * dBase;

    // accrued interest and forfeit coupon only apply for strikes
    if (m_bKeepAccrued)
      dStrike += dAccruedAmount;

    if (!m_bForfeitCoupon)
      dStrike += m_pParams->GetCouponAmount();
  }

  
  for (size_t nIdxSpot = 0; nIdxSpot < nNbS; nIdxSpot++)
    pdValues[nIdxSpot] = dStrike;

  return true;
}


} // namespace pricing

} // namespace ito33
