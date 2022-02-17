/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/varianceswap/varianceswapinstdata.cpp
// Purpose:     Implementation of VarianceSwapInstData class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapinstdata.cpp,v 1.3 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/varianceswap/varianceswapinstdata.cpp
    @brief Implementation of VarianceSwapInstData class
*/

#include "ito33/sharedptr.h"

#include "ito33/finance/payoff.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"

#include "ihg/varianceswapinstdata.h"

namespace ito33
{

namespace ihg
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

  Alloc(m_nNbS);

  m_varianceSwapMeshes.ComputeRecoveryValues(); 
  
}

void VarianceSwapInstData::SetInitialValue()
{

  m_varianceSwapParams.GetPayoff()->Get(m_pdS, m_pdPrices.Get(), m_nNbS);
 
  m_bHasEvent = false;
}

void VarianceSwapInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_varianceSwapMeshes.GetRecoveryValue();
}


} // namespace ihg

} // namespace ito33
