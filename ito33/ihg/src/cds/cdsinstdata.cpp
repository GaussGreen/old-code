/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cds/cdsinstdata.cpp
// Purpose:     cds instdata class 
// Author:      Nabil, Wang
// Created:     2003/10/29
// RCS-ID:      $Id: cdsinstdata.cpp,v 1.8 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/cds/cdsinstdata.cpp
   @brief Implementation of InstData class for cds
 */

#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

#include "ihg/cdsinstdata.h"

namespace ito33
{

namespace ihg
{


CDSInstData::CDSInstData(pricing::CDSParams& params,
                         Model& model,
                         pricing::CDSMeshManager& meshes)
                       : BackwardInstData(params, model, meshes), 
                         m_cdsParams(params), 
                         m_cdsMeshes(meshes)
{
}

void CDSInstData::Init()
{
  // Get the space mesh
  m_pdS = m_cdsMeshes.GetS();

  m_nNbS = m_cdsMeshes.GetNbS();

  m_pdLogS = m_cdsMeshes.GetLogS();

  Alloc(m_nNbS);
}

void CDSInstData::SetInitialValue()
{
  size_t nIdxS;
 
  for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    m_pdPrices[nIdxS] = 0;  

  if (m_bComputeVega)
    for (nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdVegas[nIdxS] = 0.;

  DoEvents();
}

void CDSInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_cdsMeshes.GetRecoveryValue();
}


} // namespace ihg

} // namespace ito33

