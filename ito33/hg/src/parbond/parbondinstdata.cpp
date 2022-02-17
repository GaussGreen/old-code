/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/parbond/parbondinstdata.cpp
// Purpose:     ParBond instdata class using HG model
// Created:     2005/02/16
// RCS-ID:      $Id: parbondinstdata.cpp,v 1.1 2005/06/09 15:36:19 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "hg/model.h"
#include "hg/parbondinstdata.h"

namespace ito33
{

namespace hg
{


ParBondInstData::ParBondInstData
    (pricing::ParBondParams& params, Model& model, 
     pricing::ParBondMeshManager& meshes)
   : InstDataTimeOnly(params, model, meshes),  
     m_parBondParams(params), m_parBondMeshes(meshes)
{
}

void ParBondInstData::UpdateRecoveryValue()
{
  m_dRecoveryValue = m_parBondMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
