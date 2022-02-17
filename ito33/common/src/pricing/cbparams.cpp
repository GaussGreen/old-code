/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cbparams.cpp
// Purpose:     cb params class
// Author:      Laurence
// Created:     2004/03/16
// RCS-ID:      $Id: cbparams.cpp,v 1.105 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/payoffconstant.h"

#include "ito33/pricing/dividendevents.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbparams.h"
#include "ito33/pricing/cbevent.h"

#include "ito33/numeric/predicatetime.h"

using namespace ito33;
using namespace ito33::pricing;
using namespace ito33::numeric;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::CBParams);
}

void CBParams::Init()
{
  CBLikeParams::Init();

  // Add accretion start date in the cb event manager
  double dAccretionStartTime = GetCB().GetAccretionStartTime();
  if (   !IsBefore(dAccretionStartTime, m_dValuationTime)
       && IsBefore(dAccretionStartTime, m_dStoppingTime) )
  {
    m_CBEventManager.AddEvent
                     (shared_ptr<CBEvent>(
                        new CBEvent(dAccretionStartTime, CBET_StartAccretion)));
  }
}

CBParams* CBParams::Clone() const
{
  // Copy the underlying contract. The new params class will manage the memory
  AutoPtr<pricing::CB> clonedCB( new CB(m_cb) );

  // Construct and setup the cloned cb params
  CBParams* pClonedParams = new CBParams(clonedCB);

  pClonedParams->SetNumParams(m_pNumParams);
  pClonedParams->SetMeshParams(m_pMeshParams);

  pClonedParams->SetYieldCurve( GetYieldCurve() );
  pClonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  pClonedParams->SetForeignCurve( GetForeignCurve() );
  pClonedParams->SetDividends( GetDividends() );
  pClonedParams->SetValuationTime( GetValuationTime() );
  pClonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return pClonedParams;
}

CBParams* CBParams::GetCallNoticeParams()
{
  CBParams* 
    pCBParams = (CBParams *)(CBLikeParams::GetCallNoticeParams() );

  double dMaturityTime = m_pCallNoticeParams->GetCBLike().GetMaturityTime();

  // m_cbLike corresponds to the main cb problem, not the clone.
  // Hence, the GetClaim call uses the full coupon structure when
  // getting the claim at the new (cloned call notice) maturity time
  pCBParams->GetCB().SetRedemptionValue( m_cbLike.GetClaim(dMaturityTime) );
  
  return pCBParams;
}
