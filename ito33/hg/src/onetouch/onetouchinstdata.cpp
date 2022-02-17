/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/onetouch/onetouchinstdata.cpp
// Purpose:     OneTouch instdata class 
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchinstdata.cpp,v 1.2 2005/07/05 18:32:01 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/onetouchparams.h"
#include "ito33/pricing/onetouchmeshmanager.h"

#include "hg/onetouchinstdata.h"

namespace ito33
{

namespace hg
{
  
  using namespace numeric;
  using namespace finance;

OneTouchInstData::OneTouchInstData
(pricing::OneTouchParams& params, Model& model, 
 pricing::OneTouchMeshManager& meshes)
: BackwardInstData(params, model, meshes), 
  m_oneTouchParams(params), m_oneTouchMeshes(meshes)
{
}

void OneTouchInstData::Init()
{
  // Get the space mesh
  m_pdS = m_oneTouchMeshes.GetS();

  m_nNbS = m_oneTouchMeshes.GetNbS();

  m_pdLogS = m_oneTouchMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  // The sparse matrix for the jumps needs to know the boundary type.
  // The Dirichlet value does not matter. It is set again in UpdateBeforeStep
  m_dRecoveryValue = 0.0;

  if ( m_oneTouchParams.GetOneTouch().GetBarrierType() == Barrier_UpAndOut )
    m_BoundaryCondition.SetRight(BCType_Dirichlet, m_dRecoveryValue);
  else
    m_BoundaryCondition.SetLeft(BCType_Dirichlet, m_dRecoveryValue);
}

void OneTouchInstData::SetInitialValue()
{
  size_t nIdxS;

  // todo: We might be able to smooth out a bit the initial condition
  if ( m_oneTouchParams.GetOneTouch().GetBarrierType() == Barrier_UpAndOut )
  {
    m_pdPrices[m_nNbS - 1] = 1.;

    for (nIdxS = 0; nIdxS < m_nNbS - 1; nIdxS++)
      m_pdPrices[nIdxS] = 0;  
  }
  else
  {
    m_pdPrices[0] = 1.;

    for (nIdxS = 1; nIdxS < m_nNbS; nIdxS++)
      m_pdPrices[nIdxS] = 0;  
  }

  for (size_t nIdxR = 1; nIdxR < m_nNbRegimes; nIdxR++)
    for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdPrices[nIdxR * m_nNbS + nIdxS] = m_pdPrices[nIdxS];

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];

  SetInitialSensitivityValue();

  DoEvents();
}

void OneTouchInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_oneTouchMeshes.GetRecoveryValue();

  double dBoundaryValue = m_oneTouchMeshes.GetBoundaryValue();

  // At the barrier, the boundary condition is Dirichlet 
  if ( m_oneTouchParams.GetOneTouch().GetBarrierType() == Barrier_UpAndOut )
    m_BoundaryCondition.SetRight(BCType_Dirichlet, dBoundaryValue);
  else
    m_BoundaryCondition.SetLeft(BCType_Dirichlet, dBoundaryValue);
}


} // namespace hg

} // namespace ito33
