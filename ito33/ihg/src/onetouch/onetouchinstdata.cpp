/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/onetouch/onetouchinstdata.cpp
// Purpose:     OneTouch instdata class 
// Created:     2006/08/10
// RCS-ID:      $Id: onetouchinstdata.cpp,v 1.1 2006/08/10 23:12:02 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/onetouchparams.h"
#include "ito33/pricing/onetouchmeshmanager.h"

#include "ihg/onetouchinstdata.h"

namespace ito33
{

namespace ihg
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

  Alloc(m_nNbS);

  // Dirichlet Boundary condition needs to be set before any call to 
  // numoutput::UpdateMe.  The value does not matter.
  if ( m_oneTouchParams.GetOneTouch().GetBarrierType() == Barrier_UpAndOut )
    m_BoundaryCondition.SetRight(BCType_Dirichlet, 0);
  else
    m_BoundaryCondition.SetLeft(BCType_Dirichlet, 0);
}

void OneTouchInstData::SetInitialValue()
{
  // todo: We might be able to smooth out a bit the initial condition
  if ( m_oneTouchParams.GetOneTouch().GetBarrierType() == Barrier_UpAndOut )
  {
    m_pdPrices[m_nNbS - 1] = 1.;

    for (size_t nIdxS = 0; nIdxS < m_nNbS - 1; nIdxS++)
      m_pdPrices[nIdxS] = 0;  
  }
  else
  {
    m_pdPrices[0] = 1.;

    for (size_t nIdxS = 1; nIdxS < m_nNbS; nIdxS++)
      m_pdPrices[nIdxS] = 0;  
  }

  if (m_bComputeVega)
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdVegas[nIdxS] = 0.;

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

} // namespace ihg

} // namespace ito33

