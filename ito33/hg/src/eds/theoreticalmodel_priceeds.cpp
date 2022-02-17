/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoriticalmodel_priceeds.cpp
// Purpose:     Implementation of EDS pricing using HG model
// Created:     2005/01/31
// RCS-ID:      $Id: theoreticalmodel_priceeds.cpp,v 1.13 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/eds.h"

#include "ito33/pricing/eds.h"
#include "ito33/pricing/edsparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "hg/model.h"
#include "hg/backwardnumoutput.h"
#include "hg/edspricer.h"
#include "hg/license.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_THIS_MODULE(HGPriceEDS);

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceEDS(const finance::EDS& eds) 
{
  double dFact;
  
  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *eds.GetSessionData();

  pricing::EDS edsPricing(eds);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(eds) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::EDSParams params(edsPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  EDSPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<BackwardNumOutput> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::EDS> 
  regEDS(&TheoreticalModel::PriceEDS);


} // namespace hg

} // namespace ito33

