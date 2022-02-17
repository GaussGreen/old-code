/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/backwardmeshmanager.cpp
// Purpose:     base mesh manager for backward PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: backwardmeshmanager.cpp,v 1.14 2004/10/08 15:59:28 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/backwardmeshmanager.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

void BackwardMeshManager::ConstructTimeMesh(SpecialTimes& specialTimes)
{
  // Basic sanity checks
  ASSERT_MSG(specialTimes.size() >= 2, "At least 2 time points are required");

  // Check if a constant mesh is requested and act accordingly
  if (m_params.GetMeshParams()->GetUniformTimeGrid() == true)
    MeshManager::ConstructUniformTimeMesh(specialTimes);
  else
    MeshManager::ConstructNonuniformTimeMesh(specialTimes, 0);

}

void BackwardMeshManager::SetupSchemeTypes(SpecialTimes& specialTimes)
{
  // Predetermine the default scheme that will be used at each step
  m_pSchemeTypes = Array<SchemeType>(m_nNbTimes);

  size_t nIdxTime;

  SchemeType schemeType = m_params.GetNumParams()->GetSchemeType();
  for (nIdxTime = 0; nIdxTime < m_nNbTimes; nIdxTime++)
    m_pSchemeTypes[nIdxTime] = schemeType; 

  // The first two steps should be implicit to avoid oscillations, and to
  // allow multistep methods to have enough data.  The last time point is 
  // never called for a scheme
  m_pSchemeTypes[m_nNbTimes - 2] = SchemeType_Implicit; 
  m_pSchemeTypes[m_nNbTimes - 3] = SchemeType_Implicit;
  
  // Take special points (eg. events) into account, as well as the start time
  SpecialTimes::iterator iterTimes = specialTimes.begin();

  for (nIdxTime = 0; 
       nIdxTime < m_nNbTimes && iterTimes != specialTimes.end(); 
       ++nIdxTime)
  {
    if ( AreTimesEqual(iterTimes->GetTime(), m_pdTimes[nIdxTime]) )
    {
      if (nIdxTime > 0)
        m_pSchemeTypes[nIdxTime - 1] = SchemeType_Implicit;

      if (nIdxTime > 1)
        m_pSchemeTypes[nIdxTime - 2] = SchemeType_Implicit;

      // Move to the next special time, note that iterTimes is only advanced 
      // when the times are matched
      ++iterTimes;
    }
  }
}

void BackwardMeshManager::SetupRates()
{
  // Allocate the storage for the final rate arrays
  m_pdRates = Array<double>(m_nNbTimes);
  m_pdForeignRates = Array<double>(m_nNbTimes);
 
  Array<double>
    pdRateTmp = Array<double>(m_nNbTimes),
    pdForeignRateTmp = Array<double>(m_nNbTimes); 
  
  /*
    $ e^{ - \int_t^T r(s) ds } $ is a trivial solution of the backward PDE,
    so does $ e^{ \int_{ t_{ref} }^t r(s) ds } $ with $t$ the time variable
    since $ e^{ \int_{ t_{ref} }^t r(s) ds } 
          = e^{ - \int_t^T r(s) ds } * e^{ \int_{ t_{ref} }^T r(s) ds $
  */
  m_params.GetYieldCurve()->GetCompoundFactor
    (m_pdTimes.Get(), pdRateTmp.Get(), m_nNbTimes);

  m_params.GetForeignCurve()->GetCompoundFactor
    (m_pdTimes.Get(), pdForeignRateTmp.Get(), m_nNbTimes);
  
  double
    dInvDeltaT, dOldInvDeltaT = 0.,
    dTmp1, dTmp2;
  
  for (size_t nIdxTime = m_nNbTimes - 2; nIdxTime < m_nNbTimes - 1; nIdxTime--)
  {
    dInvDeltaT = 1.0 / ( m_pdTimes[nIdxTime + 1] - m_pdTimes[nIdxTime] );

    if (m_pSchemeTypes[nIdxTime] == SchemeType_ThreeLevel)
    {
      double dC1 = (dInvDeltaT + dOldInvDeltaT);
      double dC2 = dOldInvDeltaT / dInvDeltaT / ( m_pdTimes[nIdxTime + 2] - 
                   m_pdTimes[nIdxTime] );

      dTmp1 = pdRateTmp[nIdxTime + 1] / pdRateTmp[nIdxTime];   
      dTmp2 = pdRateTmp[nIdxTime + 2] / pdRateTmp[nIdxTime + 1]; 

      m_pdRates[nIdxTime] = dTmp1 * (dC1 - dC2 * dTmp2) - dC1 + dC2;
      
      dTmp1 = pdForeignRateTmp[nIdxTime + 1] / pdForeignRateTmp[nIdxTime];   
      dTmp2 = pdForeignRateTmp[nIdxTime + 2] / pdForeignRateTmp[nIdxTime + 1];   

      m_pdForeignRates[nIdxTime] = dTmp1 * (dC1 - dC2 * dTmp2) - dC1 + dC2;
    }
    else if (m_pSchemeTypes[nIdxTime] == SchemeType_Implicit)
    {
      dTmp1 = pdRateTmp[nIdxTime + 1] / pdRateTmp[nIdxTime];   
      m_pdRates[nIdxTime] = (dTmp1 - 1.) * dInvDeltaT;

      dTmp1 = pdForeignRateTmp[nIdxTime + 1] / pdForeignRateTmp[nIdxTime]; 
      m_pdForeignRates[nIdxTime] = (dTmp1 - 1.) * dInvDeltaT;
    }
    else if (m_pSchemeTypes[nIdxTime] == SchemeType_CrankNicolson)
    {
      dTmp1 = pdRateTmp[nIdxTime + 1] / pdRateTmp[nIdxTime];  
      dTmp2 = m_pdRates[nIdxTime + 1];

      m_pdRates[nIdxTime] = (dTmp1 * (1. - 0.5 * dTmp2 / dInvDeltaT)  - 1.) 
                          * 2.0 * dInvDeltaT;      

      dTmp1 = pdForeignRateTmp[nIdxTime + 1] / pdForeignRateTmp[nIdxTime]; 
      dTmp2 = m_pdForeignRates[nIdxTime + 1]; 

      m_pdForeignRates[nIdxTime] = (dTmp1 * (1. - 0.5 * dTmp2 / dInvDeltaT)  - 1.) 
                                 * 2.0 * dInvDeltaT;
    }

    dOldInvDeltaT = dInvDeltaT;
  }
}


