/////////////////////////////////////////////////////////////////////////////
// Name:        forwardcds/forwardcdsinstdata.cpp
// Purpose:     implementation of cds instdata class using forward PDE 
// Author:      David
// RCS-ID:      $Id: forwardcdsinstdata.cpp,v 1.10 2006/04/21 09:26:06 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/payoff.h"

#include "ito33/pricing/forwardcdsparams.h"
#include "ito33/pricing/forwardcdsmeshmanager.h"

#include "ihg/forwardcdsinstdata.h"

namespace ito33
{

namespace ihg
{

ForwardCDSInstData::ForwardCDSInstData(pricing::ForwardCDSParams &params, 
	                 Model &model,
                   pricing::ForwardCDSMeshManager &meshes)
                 : InstData(params, model, meshes), 
	                 m_forwardCDSParams(params), m_forwardCDSMeshes(meshes)
{
  m_bIsHRTimeOnly = false;
}

void ForwardCDSInstData::Alloc(size_t nNbS)
{
  InstData::Alloc(nNbS);
}


void ForwardCDSInstData::Init()
{
  // Get the space mesh
  m_pdS = m_forwardCDSMeshes.GetS();

  m_nNbS = m_forwardCDSMeshes.GetNbS();
										                		
  m_pdLogS = m_forwardCDSMeshes.GetLogS();

  Alloc(m_nNbS); 

  m_nCurrentIndex = 0;
}


void ForwardCDSInstData::SetInitialValue()
{
  // Setup initial conditions.    
  m_forwardCDSParams.GetPayoff()->Get(m_pdS, m_pdPrices.Get(), m_nNbS);

  m_dRecoveryTermOld = 0.0;
  m_dSpreadTermOld = 0.0;
  m_dAccruedTermOld = 0.0;
}


void ForwardCDSInstData::UpdateBeforeStep()
{
  InstData::UpdateBeforeStep();

  m_bIsHRTimeOnly = m_model.IsHazardRateTimeOnly();

  m_nCurrentIndex = m_forwardCDSMeshes.GetCurrentIndex();

  m_dRecovery = 1.0 - m_forwardCDSParams.GetCDSes().GetRecoveryRate();

  m_dSpread = m_forwardCDSMeshes.GetSpread();

  m_dAccruedFraction = m_forwardCDSMeshes.GetAccruedFraction();

}

} // namespace ihg

} // namespace ito33
