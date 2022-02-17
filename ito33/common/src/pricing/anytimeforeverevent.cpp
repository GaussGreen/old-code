/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/anytimeforever.cpp
// Purpose:     anytime forever event
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: anytimeforeverevent.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/binarysearch.h"

#include "ito33/pricing/anytimeforeverevent.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"
#include "ito33/pricing/conversionprovisions.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"

namespace ito33
{

namespace pricing
{

    
void AnyTimeForeverEvent::ApplyAtStartTime(PathDepStructure& /*pathDepStruct*/) const
{

}

void AnyTimeForeverEvent::ApplyAtEndTime(PathDepStructure& /*pathDepStruct*/) const
{

}

void AnyTimeForeverEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{
  
  // Make sure path 0 is active before getting the time
  ASSERT_MSG(pathDepStruct.m_pbIsActive[0] == true,
    "Path 0 is not active.");

  // All paths must be at the same time. The trigger is the same for
  // all paths
  double dTime = pathDepStruct.m_path[0].meshes->GetTime();
  finance::TriggerAsPercentageOf
    trigger = m_pConversionProvisions->GetTriggerAsPercentageOf();

  double dTrigger = m_pConversionProvisions->GetConversionPrice(dTime,trigger);
  dTrigger *= m_dTriggerRate;

  // Loop over the call window (if any)
  size_t nIdxCall;
  for (nIdxCall = 0; nIdxCall < m_nCallPaths; nIdxCall++)
  {
    // Path dependent means at least two problems to solve
    ASSERT_MSG(pathDepStruct.m_pbIsActive[nIdxCall] == true,
      "Base path is not active.");

    //copy anytime/forever path into base path if S is above trigger
    size_t nGridSizeBase     = pathDepStruct.m_path[nIdxCall].meshes->GetNbS();
    const double* pdGridBase = pathDepStruct.m_path[nIdxCall].meshes->GetS();
    double* pdValuesBase     = pathDepStruct.m_path[nIdxCall].instdata->m_pdPrices.Get();
  
    size_t nIdxForever = 2*m_nCallPaths + nIdxCall;

    ASSERT_MSG(pathDepStruct.m_pbIsActive[nIdxForever] == true,
      "Anytime/forever path is not active.");

    size_t nGridSizeForever     = pathDepStruct.m_path[nIdxForever].meshes->GetNbS();
    const double* pdGridForever = pathDepStruct.m_path[nIdxForever].meshes->GetS();
    double* pdValuesForever     = pathDepStruct.m_path[nIdxForever].instdata->m_pdPrices.Get();

    Array<double> pdInterpValues(nGridSizeBase);

    Interpolate(pdGridForever, pdValuesForever, nGridSizeForever,
      pdGridBase, pdInterpValues.Get(), nGridSizeBase,
      numeric::ExtrapolationMode_Linear,
      numeric::ExtrapolationMode_Linear);

    // find the trigger in the grid    
    size_t nIdx = BinSearch(pdGridBase, nGridSizeBase, dTrigger);

    if ( nIdx > 0 
          && !numeric::LessThanTrigger(pdGridBase[nIdx - 1], dTrigger) ) 
        nIdx--;

    for (; nIdx < nGridSizeBase; nIdx++)
    {
      // Copy above the trigger
      pdValuesBase[nIdx] = pdInterpValues[nIdx];
    }
  } // loop over call paths
}


} // namespace pricing

} // namespace ito33
