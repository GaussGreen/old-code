/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/transitionprobabilityinstdata.cpp
// Purpose:     instdata class for transition probability
// Created:     2006/03/29
// RCS-ID:      $Id: transitionprobabilityinstdata.cpp,v 1.1 2006/03/31 17:43:59 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/binarysearch.h"
#include "ito33/pricing/transitionprobabilityparams.h"
#include "ito33/pricing/transitionprobabilitymeshmanager.h"

#include "hg/transitionprobabilityinstdata.h"

namespace ito33
{

namespace hg
{


TransitionProbabilityInstData::TransitionProbabilityInstData
    (pricing::TransitionProbabilityParams& params,
     Model& model,
     pricing::TransitionProbabilityMeshManager& meshes)
   : BackwardInstData(params, model, meshes), 
     m_tpParams(params), 
     m_tpMeshes(meshes)
{
  m_bDualSystemRequired = true;
}

void TransitionProbabilityInstData::Init()
{
  // Get the space mesh
  m_pdS = m_tpMeshes.GetS();

  m_nNbS = m_tpMeshes.GetNbS();

  m_pdLogS = m_tpMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);
}

void TransitionProbabilityInstData::SetInitialValue()
{
  m_nIdxSpot = BinSearch(m_pdS, m_nNbS, m_tpParams.GetSpotSharePrice());

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdPrices[nIdx] = 0;  

  double dStrike = m_tpParams.GetTransitionProbability().GetStrike();
  size_t nIdxStrike = BinSearch(m_pdS, m_nNbS, dStrike);

  m_pdPrices[nIdxStrike] = 1.;

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];

  SetInitialSensitivityValue();

  DoEvents();
}

void TransitionProbabilityInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_tpMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
