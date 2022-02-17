/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/forwardmeshmanager.cpp
// Purpose:     base mesh manager for forward PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: forwardmeshmanager.cpp,v 1.12 2004/10/08 15:59:28 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/forwardmeshmanager.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

void ForwardMeshManager::ConstructTimeMesh(SpecialTimes& specialTimes)
{
  // Basic sanity checks
  ASSERT_MSG(specialTimes.size() >= 2, "At least 2 time points are required");

  // Check if a constant mesh is requested and act accordingly
  if (m_params.GetMeshParams()->GetUniformTimeGrid() == true)
    MeshManager::ConstructUniformTimeMesh(specialTimes);
  else
    MeshManager::ConstructNonuniformTimeMesh(specialTimes, 1);
}

void ForwardMeshManager::SetupSchemeTypes(SpecialTimes &pdSpecialTimes)
{
  size_t nIdxTime;

  SchemeType schemeType = m_params.GetNumParams()->GetSchemeType();

  m_pSchemeTypes = Array<SchemeType>(m_nNbTimes);

  for (nIdxTime = 0; nIdxTime < m_nNbTimes; nIdxTime++)
    m_pSchemeTypes[nIdxTime] = schemeType; 

  // The first two steps should be implicit to avoid oscillations, and to
  // allow multistep methods to have enough data.  The first time point is 
  // never called for a scheme
  m_pSchemeTypes[1] = SchemeType_Implicit; 
  m_pSchemeTypes[2] = SchemeType_Implicit;
  
  // Take special points (eg. events) into account, as well as the start time
  SpecialTimes::iterator iterTimes = pdSpecialTimes.begin();

  for (nIdxTime = 0; 
       nIdxTime < m_nNbTimes && iterTimes != pdSpecialTimes.end(); 
       ++nIdxTime)
  {
    if ( AreTimesEqual(iterTimes->GetTime(), m_pdTimes[nIdxTime]) )
    {
      if (nIdxTime < m_nNbTimes - 1)
        m_pSchemeTypes[nIdxTime + 1] = SchemeType_Implicit;
      if (nIdxTime < m_nNbTimes - 2)
        m_pSchemeTypes[nIdxTime + 2] = SchemeType_Implicit;

      // Move to the next special time, note that iterTimes is only advanced 
      // when the times are matched
      ++iterTimes;
    }
  }
}

void ForwardMeshManager::SetupRates()
{
  // Allocate the storage for the final rate arrays
  m_pdRates = Array<double>(m_nNbTimes);
  m_pdForeignRates = Array<double>(m_nNbTimes);   
  
  // Get the exp(-r T) factors. Below, we divide by adjacent entries to
  // get exp(-r delta T) factors
  
  Array<double>
    pdRateTmp = Array<double>(m_nNbTimes),
    pdForeignRateTmp = Array<double>(m_nNbTimes);
  
  m_params.GetYieldCurve()->GetDiscountFactor
    (m_pdTimes.Get(), pdRateTmp.Get(), m_nNbTimes);

  m_params.GetForeignCurve()->GetDiscountFactor
    (m_pdTimes.Get(), pdForeignRateTmp.Get(), m_nNbTimes);
  
  double
    dInvDeltaT, dOldInvDeltaT = 0.,
    dTmp1, dTmp2;
  
  for (size_t nIdxTime = 1; nIdxTime < m_nNbTimes; nIdxTime++)
  {
    dInvDeltaT = 1.0 / ( m_pdTimes[nIdxTime] - m_pdTimes[nIdxTime - 1] );

    if (m_pSchemeTypes[nIdxTime] == SchemeType_ThreeLevel)
    {
      dTmp1 = pdRateTmp[nIdxTime - 1] / pdRateTmp[nIdxTime];   
      dTmp2 = pdRateTmp[nIdxTime - 2] / pdRateTmp[nIdxTime];   
      
      m_pdRates[nIdxTime] = 1.5 * (dTmp1 - 1.0) * dInvDeltaT
                          - 0.5 * (dTmp2 - dTmp1) * dOldInvDeltaT;

      dTmp1 = pdForeignRateTmp[nIdxTime - 1] / pdForeignRateTmp[nIdxTime];   
      dTmp2 = pdForeignRateTmp[nIdxTime - 2] / pdForeignRateTmp[nIdxTime];   

      m_pdForeignRates[nIdxTime] = 1.5 * (dTmp1 - 1.0) * dInvDeltaT
                                 - 0.5 * (dTmp2 - dTmp1) * dOldInvDeltaT;
    }
    else if (m_pSchemeTypes[nIdxTime] == SchemeType_Implicit)
    {
      dTmp1 = pdRateTmp[nIdxTime - 1] / pdRateTmp[nIdxTime];   
      m_pdRates[nIdxTime] = (dTmp1 - 1.) * dInvDeltaT;

      dTmp1 = pdForeignRateTmp[nIdxTime - 1] / pdForeignRateTmp[nIdxTime]; 
      m_pdForeignRates[nIdxTime] = (dTmp1 - 1.) * dInvDeltaT;
    }
    else if (m_pSchemeTypes[nIdxTime] == SchemeType_CrankNicolson)
    {

      dTmp1 = pdRateTmp[nIdxTime - 1] / pdRateTmp[nIdxTime];  
      dTmp2 = m_pdRates[nIdxTime - 1];

      m_pdRates[nIdxTime] = (dTmp1 * (1. - 0.5 * dTmp2 / dInvDeltaT)  - 1.)
                          * 2.0 * dInvDeltaT;

      dTmp1 = pdForeignRateTmp[nIdxTime - 1] / pdForeignRateTmp[nIdxTime]; 
      dTmp2 = m_pdForeignRates[nIdxTime - 1]; 

      m_pdForeignRates[nIdxTime] = (dTmp1 * (1. - 0.5 * dTmp2 / dInvDeltaT)  - 1.)
                                 * 2.0 * dInvDeltaT;
    }

    dOldInvDeltaT = dInvDeltaT;
  }
}


