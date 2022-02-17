/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/quarterlyevent.cpp
// Purpose:     quarterly event
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: quarterlyevent.cpp,v 1.13 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/binarysearch.h"

#include "ito33/pricing/quarterlyevent.h"
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

    
void QuarterlyEvent::ApplyAtStartTime(PathDepStructure& ) const
{
  // Nothing to do for discrete event
}

void QuarterlyEvent::TurnOffPaths(PathDepStructure& pathDepStruct) const
{
  // set the quarterly pde paths to be inactive
  size_t nIdxCall;
  for (nIdxCall = 0; nIdxCall < m_nCallPaths; nIdxCall++)
  {
    pathDepStruct.m_pbIsActive[m_nCallPaths + nIdxCall] = false;
  }

}

void QuarterlyEvent::ApplyAtEndTime(PathDepStructure& ) const
{
  // Nothing to do for discrete event
 
}

void QuarterlyEvent::Initialize(PathDepStructure& pathDepStruct) const
{
  
  // Make sure path 0 is active before copying
  ASSERT_MSG(pathDepStruct.m_pbIsActive[0] == true,
    "Path 0 is not active.");

  // copy the base pde path into the quarterly pde path for 
  // each call (day above trigger) path
  size_t nIdxCall;
  for (nIdxCall = 0; nIdxCall < m_nCallPaths; nIdxCall++)
  {
    // set the quaterly pde path to be active
    size_t nIdxQuarterly = m_nCallPaths + nIdxCall;
    pathDepStruct.m_pbIsActive[nIdxQuarterly] = true;

    size_t nGridSizeBase     = pathDepStruct.m_path[nIdxCall].meshes->GetNbS();
    const double* pdGridBase = pathDepStruct.m_path[nIdxCall].meshes->GetS();
    double* pdValuesBase     = pathDepStruct.m_path[nIdxCall].instdata->m_pdPrices.Get();     

    size_t nGridSizeQuarterly     = pathDepStruct.m_path[nIdxQuarterly].meshes->GetNbS();
    const double* pdGridQuarterly = pathDepStruct.m_path[nIdxQuarterly].meshes->GetS();
    double* pdValuesQuarterly     = pathDepStruct.m_path[nIdxQuarterly].instdata->m_pdPrices.Get();

    Interpolate(pdGridBase, pdValuesBase, nGridSizeBase,
      pdGridQuarterly, pdValuesQuarterly, nGridSizeQuarterly,
      numeric::ExtrapolationMode_Linear,
      numeric::ExtrapolationMode_Linear);
  } // loop over call paths
  
}


void QuarterlyEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{

  if (m_bInitializeOnly == true)
  {
    Initialize(pathDepStruct);
    return;
  }

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
  

  // copy the quaterly pde path into the base pde path for 
  // each call (day above trigger) path for values above the trigger
  size_t nIdxCall;
  for (nIdxCall = 0; nIdxCall < m_nCallPaths; nIdxCall++)
  {

    size_t nGridSizeBase     = pathDepStruct.m_path[nIdxCall].meshes->GetNbS();
    const double* pdGridBase = pathDepStruct.m_path[nIdxCall].meshes->GetS();
    double* pdValuesBase     = pathDepStruct.m_path[nIdxCall].instdata->m_pdPrices.Get();
  
    size_t nIdxQuarterly = m_nCallPaths + nIdxCall;
    size_t nGridSizeQuarterly     = pathDepStruct.m_path[nIdxQuarterly].meshes->GetNbS();
    const double* pdGridQuarterly = pathDepStruct.m_path[nIdxQuarterly].meshes->GetS();
    double* pdValuesQuarterly     = pathDepStruct.m_path[nIdxQuarterly].instdata->m_pdPrices.Get();

    Array<double> pdInterpValuesBase(nGridSizeBase);
 
    Interpolate(pdGridQuarterly, pdValuesQuarterly, nGridSizeQuarterly,
      pdGridBase, pdInterpValuesBase.Get(), nGridSizeBase,
      numeric::ExtrapolationMode_Linear,
      numeric::ExtrapolationMode_Linear);

    Array<double> pdInterpValuesQuarterly(nGridSizeQuarterly);

    Interpolate(pdGridBase, pdValuesBase, nGridSizeBase,
      pdGridQuarterly, pdInterpValuesQuarterly.Get(), nGridSizeQuarterly,
      numeric::ExtrapolationMode_Linear,
      numeric::ExtrapolationMode_Linear);


    // find the trigger in the grid
    // copy the quarterly pde into base pde if S is above trigger    
    size_t nIdx = BinSearch(pdGridBase, nGridSizeBase, dTrigger);

     if ( nIdx > 0 
          && !numeric::LessThanTrigger(pdGridBase[nIdx - 1], dTrigger) ) 
        nIdx--;

    for (; nIdx < nGridSizeBase; nIdx++)
    {
      // Copy above the trigger
      pdValuesBase[nIdx] = pdInterpValuesBase[nIdx];
    }

    // copy base pde into quarterly pde if S is below trigger
    // This is done to initialize the next quaterly period
    // nIdx-1 since binary search returns  A[t]<= x <= A[u]
    // returns u
    nIdx = BinSearch(pdGridQuarterly, nGridSizeQuarterly, dTrigger) - 1;

     if ( nIdx > 0 
          && !numeric::LessThanTrigger(pdGridQuarterly[nIdx - 1], dTrigger) ) 
        nIdx--;
  
    for ( ; nIdx < nGridSizeQuarterly; nIdx-- )
    {
      // Copy above the trigger
      pdValuesQuarterly[nIdx] = pdInterpValuesQuarterly[nIdx];
    }
  } // loop over call paths

  // Turn off the paths
  TurnOffPaths(pathDepStruct);
}


} // namespace pricing

} // namespace ito33
