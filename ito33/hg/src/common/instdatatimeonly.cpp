/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/instdatatimeonly.cpp
// Purpose:     instdata class for time only problem using HG model
// Created:     2005/06/09
// RCS-ID:      $Id: instdatatimeonly.cpp,v 1.1 2005/06/09 14:17:21 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/event.h"
#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"

#include "hg/model.h"
#include "hg/instdatatimeonly.h"

namespace ito33
{

namespace hg
{


InstDataTimeOnly::InstDataTimeOnly
    (pricing::Params& params,Model& model, pricing::MeshManager& meshes)
    : pricing::InstDataTimeOnly(params, meshes), m_model(model)
{
}

void InstDataTimeOnly::Init()
{
  m_nNbRegimes = m_model.GetNbRegimes();

  m_pdJumpsToDefault = m_model.GetJumpsToDefault();

  m_pdPrices = Array<double>(m_nNbRegimes);
  m_pdOldPrices = Array<double>(m_nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    m_ppdA[nIdxR1][nIdxR1] = - m_pdJumpsToDefault[nIdxR1];
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      if (nIdxR2 != nIdxR1)
      {
        const Jumps& jumps = m_model.GetJumps(nIdxR1, nIdxR2);
        Jumps::const_iterator pJump;

        m_ppdA[nIdxR1][nIdxR2] = 0;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          m_ppdA[nIdxR1][nIdxR1] -= pJump->GetIntensity();
          m_ppdA[nIdxR1][nIdxR2] += pJump->GetIntensity();
        }
      }
    }
  }

  // Sensitivity
  // If no sensitivity is required, we just set the number of sensitivity
  // to 0
  m_sensitivityDatas.clear();
  if ( m_pbComputeSensitivities.empty() )
  {
    m_nNbSensitivities = 0;
    return;
  }

  size_t nCounterF = 0;

  // No need to compute vol sensitivity, regardless the flag
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    nCounterF++;

  // Default intensities
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    if ( m_pbComputeSensitivities[nCounterF++] )
    {
      TimeOnlySensitivityData data;
      data.m_bIsJumpToDefault = true;
      data.m_dIntensity = m_pdJumpsToDefault[nIdxR1];
      data.m_nIdxR1 = data.m_nIdxR2 = nIdxR1;

      m_sensitivityDatas.push_back(data);
    }
  }

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      // inter regime jumps
      if (nIdxR2 != nIdxR1)
      {
        const Jumps& jumps = m_model.GetJumps(nIdxR1, nIdxR2);
        Jumps::const_iterator pJump;

        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          if ( m_pbComputeSensitivities[nCounterF++] )
          {
            TimeOnlySensitivityData data;
            data.m_bIsJumpToDefault = false;
            data.m_dIntensity = pJump->GetIntensity();
            data.m_nIdxR1 = nIdxR1;
            data.m_nIdxR2 = nIdxR2;
          
            m_sensitivityDatas.push_back(data);
          }

          nCounterF++;
        }
      }
      else
      {
        nCounterF += 2 * m_model.GetJumps(nIdxR1, nIdxR2).size();
      }
    }

  m_nNbSensitivities = m_sensitivityDatas.size();

  m_ppdSensitivities = Array< Array<double> >(m_nNbSensitivities);
  m_ppdOldSensitivities = Array< Array<double> >(m_nNbSensitivities);
  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  {
    m_ppdSensitivities[nIdxD] = Array<double>(m_nNbRegimes);
    m_ppdOldSensitivities[nIdxD] = Array<double>(m_nNbRegimes);
  }
}

void InstDataTimeOnly::SetInitialValue()
{
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    m_pdPrices[nIdxR1] = m_pdOldPrices[nIdxR1] = 0.;

  UpdateRecoveryValue();

  DoEvents(); 

  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
    for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      m_ppdSensitivities[nIdxD][nIdxR1] = 0.;
      m_ppdOldSensitivities[nIdxD][nIdxR1] = 0.;
    }
}

void InstDataTimeOnly::ApplyEvent(const pricing::Event* pEvent)
{
  double dS = 0.;

  if ( (pEvent->GetType()) == ito33::pricing::ET_Payment )
  {
    for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pEvent->ApplyToPrice(&dS, &m_pdPrices[nIdxR1], 1);

    pEvent->ApplyToPrice(&dS, &m_dRecoveryValue, 1);
  }
  
  // nothing to do for sensitivity
}

void InstDataTimeOnly::UpdateBeforeStep()
{
  m_dOldRate = m_dRate;
  m_dOldRecoveryValue = m_dRecoveryValue;

  pricing::InstDataTimeOnly::UpdateBeforeStep();

  UpdateRecoveryValue();
  
  // Swap the prices 
  swap(m_pdOldPrices, m_pdPrices);

  // Swap the sensitivities
  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
    swap(m_ppdOldSensitivities[nIdxD], m_ppdSensitivities[nIdxD]);
}

//todo: the rate is wrong, need to use GetContinuousRate which require the time
void InstDataTimeOnly::Run()
{   
  size_t nIdxR1, nIdxR2;

  m_dOldRate = m_dRate;
  double dRate12 = m_dRate; 
   
  double dRecoveryValue12 = 0.5 * (m_dRecoveryValue + m_dOldRecoveryValue);

  // temoprary array for price
  Array<double> pdTmp(m_nNbRegimes);
  Array<double> pdTmp1(m_nNbRegimes);

  // Temporary arrays for sensitivity
  Array< Array<double> > ppdDerivTmp(m_nNbSensitivities);
  Array< Array<double> > ppdDerivTmp1(m_nNbSensitivities);
  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  {
    ppdDerivTmp[nIdxD] = Array<double>(m_nNbRegimes);
    ppdDerivTmp1[nIdxD] = Array<double>(m_nNbRegimes);
  }

  // First term: y_n 
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    pdTmp[nIdxR1] = m_pdOldPrices[nIdxR1];

    // Add the contribution
    m_pdPrices[nIdxR1] = pdTmp[nIdxR1];
  }

  // Sensitivity. First term
  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      ppdDerivTmp[nIdxD][nIdxR1] = m_ppdOldSensitivities[nIdxD][nIdxR1];

      // Add the contribution
      m_ppdSensitivities[nIdxD][nIdxR1] = ppdDerivTmp[nIdxD][nIdxR1];
    }

  // k_1 = f(x_n, y_n)
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    pdTmp1[nIdxR1] = m_pdJumpsToDefault[nIdxR1] * m_dOldRecoveryValue 
                   - m_dOldRate * pdTmp[nIdxR1];

    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      pdTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    
    pdTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    m_pdPrices[nIdxR1] += pdTmp1[nIdxR1] / 6.; 
  }

  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  { 
    const TimeOnlySensitivityData& data = m_sensitivityDatas[nIdxD];

    double* pdDerivTmp = ppdDerivTmp[nIdxD].Get();
    double* pdDerivTmp1 = ppdDerivTmp1[nIdxD].Get();
  
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      pdDerivTmp1[nIdxR1] = - m_dOldRate * pdDerivTmp[nIdxR1];

      for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
        pdDerivTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdDerivTmp[nIdxR2];
    }

    nIdxR1 = data.m_nIdxR1;
    if ( m_sensitivityDatas[nIdxD].m_bIsJumpToDefault )
      pdDerivTmp1[nIdxR1] += m_dOldRecoveryValue - pdTmp[nIdxR1];
    else
      pdDerivTmp1[nIdxR1] += pdTmp[data.m_nIdxR2] - pdTmp[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      m_ppdSensitivities[nIdxD][nIdxR1] += pdDerivTmp1[nIdxR1] / 6.;
  }
  
  // Second term: k_2 = f(x_n + 0.5 * h, y_n + 0.5 * k_1)
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    pdTmp[nIdxR1] = m_pdOldPrices[nIdxR1] + 0.5 * pdTmp1[nIdxR1];

  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    pdTmp1[nIdxR1] = m_pdJumpsToDefault[nIdxR1] * dRecoveryValue12
                   - dRate12 * pdTmp[nIdxR1];

    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      pdTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    
    pdTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    m_pdPrices[nIdxR1] += pdTmp1[nIdxR1] / 3.;
  }

  // Sensitivity. Second term: k_2 = f(x_n + 0.5 * h, y_n + 0.5 * k_1)
  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  {
    const TimeOnlySensitivityData& data = m_sensitivityDatas[nIdxD];

    double* pdDerivTmp = ppdDerivTmp[nIdxD].Get();
    double* pdDerivTmp1 = ppdDerivTmp1[nIdxD].Get();
    double* pdOldDeriv = m_ppdOldSensitivities[nIdxD].Get();
    
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp[nIdxR1] = pdOldDeriv[nIdxR1] + 0.5 * pdDerivTmp1[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      pdDerivTmp1[nIdxR1] = - dRate12 * pdDerivTmp[nIdxR1];

      for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
        pdDerivTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdDerivTmp[nIdxR2];
    }

    nIdxR1 = data.m_nIdxR1;
    if ( m_sensitivityDatas[nIdxD].m_bIsJumpToDefault )
      pdDerivTmp1[nIdxR1] += dRecoveryValue12 - pdTmp[nIdxR1];
    else
      pdDerivTmp1[nIdxR1] += pdTmp[data.m_nIdxR2] - pdTmp[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      m_ppdSensitivities[nIdxD][nIdxR1] += pdDerivTmp1[nIdxR1] / 3.;
  }

  // Third term: k3 = f(x_n + 0.5 * h, y_n + 0.5 * k_2)
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    pdTmp[nIdxR1] = m_pdOldPrices[nIdxR1] + 0.5 * pdTmp1[nIdxR1];

  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    pdTmp1[nIdxR1] = m_pdJumpsToDefault[nIdxR1] * dRecoveryValue12 
                   - dRate12 * pdTmp[nIdxR1];

    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      pdTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    
    pdTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    m_pdPrices[nIdxR1] += pdTmp1[nIdxR1] / 3.;
  }

  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  {
    const TimeOnlySensitivityData& data = m_sensitivityDatas[nIdxD];

    double* pdDerivTmp = ppdDerivTmp[nIdxD].Get();
    double* pdDerivTmp1 = ppdDerivTmp1[nIdxD].Get();
    double* pdOldDeriv = m_ppdOldSensitivities[nIdxD].Get();
    
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp[nIdxR1] = pdOldDeriv[nIdxR1] + 0.5 * pdDerivTmp1[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      pdDerivTmp1[nIdxR1] = - dRate12 * pdDerivTmp[nIdxR1];

      for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
        pdDerivTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdDerivTmp[nIdxR2];
    }

    nIdxR1 = data.m_nIdxR1;
    if ( m_sensitivityDatas[nIdxD].m_bIsJumpToDefault )
      pdDerivTmp1[nIdxR1] += dRecoveryValue12 - pdTmp[nIdxR1];
    else
      pdDerivTmp1[nIdxR1] += pdTmp[data.m_nIdxR2] - pdTmp[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      m_ppdSensitivities[nIdxD][nIdxR1] += pdDerivTmp1[nIdxR1] / 3.;
  }

  // fourth term: k_4 = f(x_n + h, y_n + k_3)
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    pdTmp[nIdxR1] = m_pdOldPrices[nIdxR1] + pdTmp1[nIdxR1];

  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    pdTmp1[nIdxR1] = m_pdJumpsToDefault[nIdxR1] * m_dRecoveryValue 
                   - m_dRate * pdTmp[nIdxR1];

    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      pdTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    
    pdTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    m_pdPrices[nIdxR1] += pdTmp1[nIdxR1] / 6.;
  }

  for (size_t nIdxD = 0; nIdxD < m_nNbSensitivities; nIdxD++)
  { 
    const TimeOnlySensitivityData& data = m_sensitivityDatas[nIdxD];

    double* pdDerivTmp = ppdDerivTmp[nIdxD].Get();
    double* pdDerivTmp1 = ppdDerivTmp1[nIdxD].Get();
    double* pdOldDeriv = m_ppdOldSensitivities[nIdxD].Get();

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp[nIdxR1] = pdOldDeriv[nIdxR1] + pdDerivTmp1[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    {
      pdDerivTmp1[nIdxR1] = - m_dRate * pdDerivTmp[nIdxR1];

      for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
        pdDerivTmp1[nIdxR1] += m_ppdA[nIdxR1][nIdxR2] * pdDerivTmp[nIdxR2];
    }

    nIdxR1 = data.m_nIdxR1;
    if ( m_sensitivityDatas[nIdxD].m_bIsJumpToDefault )
      pdDerivTmp1[nIdxR1] += m_dRecoveryValue - pdTmp[nIdxR1];
    else
      pdDerivTmp1[nIdxR1] += pdTmp[data.m_nIdxR2] - pdTmp[nIdxR1];

    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      pdDerivTmp1[nIdxR1] *= m_dTimeStep;

    // Add the contribution
    for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
      m_ppdSensitivities[nIdxD][nIdxR1] += pdDerivTmp1[nIdxR1] / 6.;
  }
}

std::vector<double> InstDataTimeOnly::GetSensitivities() const
{
  std::vector<double> pdSensitivities;

  // Vol sensitivity is just zero
  size_t nCounterF = 0;

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    if ( m_pbComputeSensitivities[nCounterF++] )
      pdSensitivities.push_back(0);

  size_t nIdxReal = 0;

  // Default intensity sensitivity
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    if ( m_pbComputeSensitivities[nCounterF++] )
      pdSensitivities.push_back(m_ppdSensitivities[nIdxReal++][0]);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_model.GetJumps(nIdxR1, nIdxR2);

      // Internal jumps sensitivity is just zero
      if (nIdxR1 == nIdxR2)
      {
        Jumps::const_iterator pJump;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          if ( m_pbComputeSensitivities[nCounterF++] )
            pdSensitivities.push_back(0);

          if ( m_pbComputeSensitivities[nCounterF++] )
            pdSensitivities.push_back(0);
        }
      }
      else
      {
        Jumps::const_iterator pJump;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          // intensity sensitivity
          if ( m_pbComputeSensitivities[nCounterF++] )
            pdSensitivities.push_back(m_ppdSensitivities[nIdxReal++][0]);

          // amplitude sensitivity is just zero
          if ( m_pbComputeSensitivities[nCounterF++] )
            pdSensitivities.push_back(0);  
        }
      }
    }

  return pdSensitivities;
}


} // namespace hg

} // namespace ito33
