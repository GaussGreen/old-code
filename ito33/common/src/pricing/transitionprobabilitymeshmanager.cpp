/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/transitionprobabilitymeshmanager.cpp
// Purpose:     option mesh manager for backward PDE problems
// Created:     2004/02/11
// RCS-ID:      $Id: transitionprobabilitymeshmanager.cpp,v 1.1 2006/03/31 17:43:48 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/transitionprobabilityparams.h"
#include "ito33/pricing/transitionprobabilitymeshmanager.h"

namespace ito33
{
  using namespace numeric;
  using namespace numeric::mesh;

namespace pricing
{

TransitionProbabilityMeshManager::TransitionProbabilityMeshManager
    (TransitionProbabilityParams& params, Model& model)
   : BackwardMeshManagerFix(params, model),
     m_tpParams(params)  { }

void TransitionProbabilityMeshManager::ConstructSpaceMesh()
{
  const double dStrike = m_tpParams.GetTransitionProbability().GetStrike();
  
  OptionSpaceMesh osgGrid;

  osgGrid.SetOptionType(Option_Call);

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                            ( m_params.GetStoppingTime(), dStrike);

  osgGrid.SetDiffusionSize
          ( m_params.GetDiffusionSize(dSquaredTotalVol) );

  double dConvection = m_model.GetConvection
                       ( m_params.GetStoppingTime(), dStrike );

  // Compute the convection size for this problem
  osgGrid.SetConvectionSize
          ( m_params.GetBackwardConvectionSize(dConvection) );

  // Generate a space mesh centered at zero
  std::vector<double> vecGrid;
  
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();

  osgGrid.Build(nNbPoints,
                log( m_params.GetSpotSharePrice() / dStrike ),
                vecGrid);

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  // Re-center the space mesh on the strike
  double dLogStrike = log(dStrike);
    
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx] + dLogStrike;
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
}

void TransitionProbabilityMeshManager::ComputeRecoveryValues()
{
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
    m_pdRecoveryValues[nIdxT] = 0.;
}

void TransitionProbabilityMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();
}

} // namespace pricing

} // namespace ito33
