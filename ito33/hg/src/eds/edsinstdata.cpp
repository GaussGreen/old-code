/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/eds/edsinstdata.cpp
// Purpose:     EDS instdata class 
// Created:     2005/01/31
// RCS-ID:      $Id: edsinstdata.cpp,v 1.9 2005/06/10 14:03:28 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/edsparams.h"
#include "ito33/pricing/edsmeshmanager.h"

#include "hg/edsinstdata.h"

namespace ito33
{

namespace hg
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

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  // The sparse matrix for the jumps needs to know the boundary type.
  // The Dirichlet value does not matter. It is set again in UpdateBeforeStep
  m_dRecoveryValue = 0.0;

  m_BoundaryCondition.SetLeft(numeric::BCType_Dirichlet, m_dRecoveryValue);
}

void EDSInstData::SetInitialValue()
{
  size_t nIdxS;

  // todo: We might be able to smooth out a bit the initial condition
  m_pdPrices[0] = 1. - m_edsParams.GetEDS().GetRecoveryRate();

  for (nIdxS = 1; nIdxS < m_nNbS; nIdxS++)
    m_pdPrices[nIdxS] = 0;  

  for (size_t nIdxR = 1; nIdxR < m_nNbRegimes; nIdxR++)
    for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdPrices[nIdxR * m_nNbS + nIdxS] = m_pdPrices[nIdxS];

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];

  SetInitialSensitivityValue();

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


} // namespace hg

} // namespace ito33
