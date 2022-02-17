/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forward/forwardoptioninstdata.cpp
// Purpose:     Implementation of HG ForwardOptionInstData class
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptioninstdata.cpp,v 1.7 2006/03/20 14:54:10 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/src/forward/forwardoptioninstdata.cpp
   @brief Implementation of HG ForwardOptionInstData class
 */

#include "ito33/finance/payoff.h"

#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/forwardoptionmeshmanager.h"

#include "hg/forwardoptioninstdata.h"

namespace ito33
{

namespace hg
{

  class Model;

ForwardOptionInstData::ForwardOptionInstData
    (pricing::ForwardOptionParams& params,
     Model& model, 
     pricing::ForwardOptionMeshManager& meshes)
   : BackwardInstData(params, model, meshes), 
     m_forwardOptionParams(params), 
     m_forwardOptionMeshes(meshes)
{
}

void ForwardOptionInstData::Init()
{
  // Get the space mesh
  m_pdS = m_forwardOptionMeshes.GetS();

  m_nNbS = m_forwardOptionMeshes.GetNbS();

  m_pdLogS = m_forwardOptionMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);
}

void ForwardOptionInstData::SetupFlags(const finance::ComputationalFlags& flags)
{
  BackwardInstData::SetupFlags(flags);

  m_bComputeFugit = false;
}

void ForwardOptionInstData::SetInitialValue()
{
  double* pdPrices = m_pdPrices.Get();

  // Only the first regime has a traditional payoff.  All other regimes
  // start with zero value
  m_forwardOptionParams.GetPayoff()->Get(m_pdS, pdPrices, m_nNbS);

  for (size_t nIdx = m_nNbS; nIdx < m_nNbX; nIdx++)
    pdPrices[nIdx] = 0.0;
 
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];

  SetInitialSensitivityValue();

  m_bHasEvent = false;
}

void ForwardOptionInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_forwardOptionMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
