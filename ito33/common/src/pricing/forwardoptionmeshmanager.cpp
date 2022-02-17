/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/forwardoptionmeshmanager.cpp
// Purpose:     option mesh manager for forward PDE problems
// Created:     2004/03/04
// RCS-ID:      $Id: forwardoptionmeshmanager.cpp,v 1.21 2006/02/24 09:26:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/forwardoptionmeshmanager.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

void ForwardOptionMeshManager::ConstructSpaceMesh()
{

  if (m_params.GetMeshParams()->GetUniformSpaceGrid() == true)
  {
    ConstructUniformSpaceMesh();
    return;
  }
  
  OptionSpaceMesh osgGrid;

  osgGrid.SetOptionType( pricing::Option_Put );

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                            ( 
                              m_forwardOptionParams.GetStoppingTime(),
                              m_params.GetSpotSharePrice() 
                            );

  osgGrid.SetDiffusionSize
          ( m_params.GetDiffusionSize(dSquaredTotalVol) );

  // Compute the convection size for this problem
  osgGrid.SetConvectionSize
          ( m_params.GetForwardConvectionSize(dSquaredTotalVol) );

  // Generate a space mesh centered at zero
  std::vector<double> vecGrid;
  
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();
  nNbPoints *= 2;

  osgGrid.Build(nNbPoints,
                0., // log( m_params.GetSpotSharePrice() / m_params.GetSpotSharePrice() ),
                vecGrid);

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  // Re-center the space mesh on the spot
  double dLogSpot = log( m_params.GetSpotSharePrice() );
  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx] + dLogSpot;
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
}

void ForwardOptionMeshManager::SetupMe()
{
  ForwardMeshManagerFix::SetupMe();

  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  size_t nIdxT;
  for (nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
      m_pdRecoveryValues[nIdxT] = 0.;

}

void ForwardOptionMeshManager::ConstructUniformSpaceMesh()
{

  // Force the points 0.0 and the spot.  Continue until approximately
  // five time the strike.  This means about 1/5 of the points will be
  // between 0 and the strike.

  m_nNbS = m_params.GetNumParams()->GetNbSpaceSteps();
  size_t nNbTmp = m_nNbS / 5;
  double dDeltaT = m_params.GetSpotSharePrice() / nNbTmp;
  
  // Construct the actual mesh
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS[0] = dDeltaT/100.0;
  m_pdLogS[0] = log(dDeltaT / 100.0);
  for (size_t nIdx = 1; nIdx < m_nNbS; nIdx++)
  {    
    m_pdS[nIdx] = nIdx * dDeltaT;
    m_pdLogS[nIdx] = log(m_pdS[nIdx]);
  }

}
