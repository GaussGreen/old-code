/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/resetparams.cpp
// Purpose:     Reset params class
// Author:      David and Yann
// Created:     2004/11/03
// RCS-ID:      $Id: resetparams.cpp,v 1.13 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/reset.h"
#include "ito33/pricing/resetparams.h"
#include "ito33/pricing/reseteventoned.h"

#include "ito33/numeric/predicatetime.h"

using namespace ito33;
using namespace ito33::pricing;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::ResetParams);
}


void ResetParams::Init()
{
  CBParams::Init();
}


ResetParams* ResetParams::Clone() const
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<Reset> clonedReset( new Reset(m_reset) );

  // Construct and setup the cloned reset params
  ResetParams* pClonedParams = new ResetParams(clonedReset); 

  pClonedParams->SetNumParams(      m_pNumParams );
  pClonedParams->SetMeshParams(     m_pMeshParams );
  pClonedParams->SetYieldCurve(     GetYieldCurve() );
  pClonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  pClonedParams->SetForeignCurve(   GetForeignCurve() );
  pClonedParams->SetDividends(      GetDividends() );
  pClonedParams->SetValuationTime(  GetValuationTime() );
  pClonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return pClonedParams;
} // ResetParams::Clone()

bool ResetParams::HasActiveResetDate() const
{
  // Check if there is an event between the valuation time and maturity time
  size_t nNbDates = m_reset.GetResetTimes().size();
  for (size_t nIdx = 0; nIdx < nNbDates; nIdx++)
  {
    double dTime = (m_reset.GetResetTimes())[nIdx];
    if (   numeric::IsAfter(dTime, m_dValuationTime) 
        && numeric::IsBefore(dTime, m_dStoppingTime ))
      return true;
  } // loop over reset times

  return false;
}

void ResetParams::ConstructResetEvents()
{ 
  const std::vector<double>& pdCapRates    = m_reset.GetCapRates();
  const std::vector<double>& pdFloorRates  = m_reset.GetFloorRates();
  const std::vector<double>& pdMultipliers = m_reset.GetMultipliers();
  const std::vector<double>& pdResetTimes  = m_reset.GetResetTimes();

  double dNominal = GetCBLike().GetNominal();

  size_t nSize = pdResetTimes.size();
  size_t nIdx;

  for (nIdx = 0; nIdx < nSize ; nIdx++)
  {
    double dEventTime = pdResetTimes[nIdx];

    // resets at valuation date and maturity are not taken into account
    if (   numeric::IsAfter(dEventTime, m_dValuationTime) 
        && numeric::IsBefore(dEventTime, m_dStoppingTime) )
    {               

      // add to base event manager (will be applied to price)
      m_eventManager.AddEvent
                     (shared_ptr<Event>(
                       new ResetEvent1D(dEventTime,
                                        pdCapRates[nIdx],
                                        pdFloorRates[nIdx],
                                        pdMultipliers[nIdx] * GetFXRate(dEventTime),
                                        dNominal) 
                       ));
                         
    } //end if
  } //end loop

} //ResetParams::ConstructResetEvents()


double ResetParams::GetFirstActiveResetTime() const
{ 
  // Simply loop through the reset times and return the first one after the
  // valuation time
  double dFirstResetTime = 0.0;
  double dValuationTime = GetValuationTime();
  const std::vector<double>& pdResetTimes = GetReset().GetResetTimes();
  size_t nNbTimes = pdResetTimes.size();
  for (size_t nIdx = 0; nIdx < nNbTimes; nIdx++)
  {
    double dTime = pdResetTimes[nIdx];
    if ( numeric::IsAfter(dTime, dValuationTime) )
    {
      dFirstResetTime = dTime;
      break;
    }
  } // loop over reset times  

  return dFirstResetTime;
}
