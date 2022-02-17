/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/attachedwarrantcbparams.cpp
// Purpose:     Attached warrant cb params
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: attachedwarrantcbparams.cpp,v 1.7 2005/12/27 10:46:09 nabil Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/pricing/attachedwarrantcb.h"
#include "ito33/pricing/attachedwarrantcbparams.h"

using namespace ito33;
using namespace ito33::pricing;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::AttachedWarrantConvertibleBondParams);
}

bool AttachedWarrantConvertibleBondParams::HasPutCallReset() const
{
  // Need to cast to cbcalls to get access to the call strikes
  CBCalls* pCalls = m_warrant.GetCBCalls();
  CBPuts* pPuts = GetPuts();

  double dResetTime = GetResetTime();

  // check if there is a put on the reset date
  size_t nIdx;
  bool bHasPut = false;
  double dPutStrike = 0.0;
  double dPutYield = 0.0;
  for (nIdx = 0; nIdx < pPuts->GetNbPuts(); nIdx++)
  {
    if ( numeric::IsEqual(pPuts->GetTime(nIdx), dResetTime) )
    {
      bHasPut = true;
      dPutStrike = pPuts->GetStrikeRate(nIdx);
      dPutYield = pPuts->GetYieldToPut(nIdx);
      break;
    }
  } // loop over puts

  // check if there is a call on the reset date
  bool bHasCall = false;
  double dCallStrike = 0.0;
  double dCallYield = 0.0;
  for (nIdx = 0; nIdx < pCalls->GetNbCalls(); nIdx++)
  {
    if ( numeric::IsEqualOrGreater(dResetTime, pCalls->GetStartTime(nIdx) )
      && numeric::IsEqualOrLess(dResetTime, pCalls->GetEndTime(nIdx)) )
    {
      bHasCall = true;
      dCallStrike = pCalls->GetStrikeRate(nIdx);
      dCallYield = pCalls->GetYieldToCall(nIdx);
      break;
    }
  } // loop over calls

  // Check if it had a put and call on the reset date
  if ( !bHasPut || !bHasCall )
    return false;

  // For a put and call to always force conversion, the strikes or yields
  // must be the same.  Check for yield first
  if (dPutStrike <= 0.0 && dCallStrike <= 0.)
    return numeric::IsEqual(dPutYield, dCallYield);

  // Otherwise, check the strikes
  return numeric::IsEqual(dPutStrike, dCallStrike);
}

AttachedWarrantConvertibleBondParams*
AttachedWarrantConvertibleBondParams::Clone() const
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<pricing::AttachedWarrantConvertibleBond> 
    clonedWarrant( new AttachedWarrantConvertibleBond(m_warrant) );

  // Construct and setup the cloned cb params
  AttachedWarrantConvertibleBondParams*
    pClonedParams = new AttachedWarrantConvertibleBondParams(clonedWarrant);

  pClonedParams->SetNumParams(      m_pNumParams );
  pClonedParams->SetMeshParams(     m_pMeshParams );
  pClonedParams->SetYieldCurve(     GetYieldCurve() );
  pClonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  pClonedParams->SetForeignCurve(   GetForeignCurve() );
  pClonedParams->SetDividends(      GetDividends() );
  pClonedParams->SetValuationTime(  GetValuationTime() );
  pClonedParams->SetSpotSharePrice( GetSpotSharePrice() );
  pClonedParams->SetStoppingTime(   m_dStoppingTime );

  return pClonedParams;
} // AttachedWarrantConvertibleBondParams::Clone()
