/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cds/cdsinstdata.cpp
// Purpose:     CDS instdata class using HG model
// Created:     2005/02/16
// RCS-ID:      $Id: cdsinstdata.cpp,v 1.6 2005/06/09 14:16:49 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/event.h"
#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

#include "hg/model.h"
#include "hg/cdsinstdata.h"

namespace ito33
{

namespace hg
{


CDSInstData::CDSInstData
    (pricing::CDSParams& params, Model& model, pricing::CDSMeshManager& meshes)
   : InstDataTimeOnly(params, model, meshes),  
     m_cdsParams(params), m_cdsMeshes(meshes)
{
}

void CDSInstData::UpdateRecoveryValue()
{
  m_dRecoveryValue = m_cdsMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
