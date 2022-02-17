/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/varianceswapinstdata.cpp
// Purpose:     Implementation of HG variance swap instdata class
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswapinstdata.cpp,v 1.3 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/src/option/optioninstdata.cpp
   @brief Implementation of HG OptionInstData class
 */

#include "ito33/sharedptr.h"

#include "ito33/finance/payoff.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"

#include "hg/varianceswapinstdata.h"

namespace ito33
{

namespace hg
{

  class Model;

VarianceSwapInstData::VarianceSwapInstData(pricing::VarianceSwapParams& params,
                                           Model& model, 
                                           pricing::VarianceSwapMeshManager& meshes)
                             : BackwardInstData(params, model, meshes), 
                               m_varianceSwapParams(params), 
                               m_varianceSwapMeshes(meshes)
{
}

void VarianceSwapInstData::Init()
{
  // Get the space mesh
  m_pdS = m_varianceSwapMeshes.GetS();

  m_nNbS = m_varianceSwapMeshes.GetNbS();

  m_pdLogS = m_varianceSwapMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  m_varianceSwapMeshes.ComputeRecoveryValues(); 
  
}

void VarianceSwapInstData::SetInitialValue()
{
  double* pdPrices = m_pdPrices.Get();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    m_varianceSwapParams.GetPayoff()->Get(m_pdS, pdPrices + nIdxR * m_nNbS, m_nNbS);
 
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];
  }

  SetInitialSensitivityValue();

  m_bHasEvent = false;
}

void VarianceSwapInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_varianceSwapMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
