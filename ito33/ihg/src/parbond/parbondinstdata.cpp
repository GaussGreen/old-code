/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/parbond/parbondinstdata.cpp
// Purpose:     parbond instdata class 
// Author:      Nabil, ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondinstdata.cpp,v 1.1 2005/06/08 16:00:08 zhang Exp $
// Copyright:   (c) 2003-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/parbond/parbondinstdata.cpp
   @brief Implementation of InstData class for parbond
 */

#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/parbondinstdata.h"

namespace ito33
{

namespace ihg
{


ParBondInstData::ParBondInstData(pricing::ParBondParams& params,
                         Model& model,
                         pricing::ParBondMeshManager& meshes)
                       : BackwardInstData(params, model, meshes), 
                         m_parbondParams(params), 
                         m_parbondMeshes(meshes)
{
}

void ParBondInstData::Init()
{
  // Get the space mesh
  m_pdS = m_parbondMeshes.GetS();

  m_nNbS = m_parbondMeshes.GetNbS();

  m_pdLogS = m_parbondMeshes.GetLogS();

  Alloc(m_nNbS);
}

void ParBondInstData::SetInitialValue()
{
  size_t nIdxS;
 
  for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    m_pdPrices[nIdxS] = 1;  

  if (m_bComputeVega)
    for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdVegas[nIdxS] = 0.;

  DoEvents();
}

void ParBondInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_parbondMeshes.GetRecoveryValue();
}


} // namespace ihg

} // namespace ito33

