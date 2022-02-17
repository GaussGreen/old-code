/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/eds/edsinstdata.cpp
// Purpose:     EDS instdata class 
// Created:     2005/01/26
// RCS-ID:      $Id: edsinstdata.cpp,v 1.2 2005/11/04 17:07:42 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/edsparams.h"
#include "ito33/pricing/edsmeshmanager.h"

#include "ihg/edsinstdata.h"

namespace ito33
{

namespace ihg
{


EDSInstData::EDSInstData(pricing::EDSParams& params,
                         Model& model,
                         pricing::EDSMeshManager& meshes)
                       : BackwardInstData(params, model, meshes), 
                         m_edsParams(params), 
                         m_edsMeshes(meshes)
{
}

void EDSInstData::Init()
{
  // Get the space mesh
  m_pdS = m_edsMeshes.GetS();

  m_nNbS = m_edsMeshes.GetNbS();

  m_pdLogS = m_edsMeshes.GetLogS();

  Alloc(m_nNbS);

  // Boundary condition at left should always be Dirichlet, and needs to be 
  // set before any call to numoutput::UpdateMe.  The value does not matter.
  m_BoundaryCondition.SetLeft(numeric::BCType_Dirichlet, 0.0);
}

void EDSInstData::SetInitialValue()
{
  size_t nIdxS;

  // todo: We might be able to smooth out a bit the initial condition
  m_pdPrices[0] = 1. - m_edsParams.GetEDS().GetRecoveryRate();

  for (nIdxS = 1; nIdxS < m_nNbS; nIdxS++)
    m_pdPrices[nIdxS] = 0;  

  if (m_bComputeVega)
    for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdVegas[nIdxS] = 0.;

  DoEvents();
}

void EDSInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_edsMeshes.GetRecoveryValue();

  // At the barrier(lower), the boundary condition is Dirichlet with value
  // equal to the recovery value
  m_BoundaryCondition.SetLeft(numeric::BCType_Dirichlet, m_dRecoveryValue);
}


} // namespace ihg

} // namespace ito33

