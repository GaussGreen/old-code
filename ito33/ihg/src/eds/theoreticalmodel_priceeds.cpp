/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_priceeds.cpp
// Purpose:     Implementation of EDS pricing using ihg model
// Created:     2005/01/26
// RCS-ID:      $Id: theoreticalmodel_priceeds.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/eds.h"

#include "ito33/pricing/eds.h"
#include "ito33/pricing/edsparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/edsnumoutput.h"
#include "ihg/edspricer.h"
#include "ihg/license.h"

#include "ito33/ihg/modeloutput.h"
#include "ito33/ihg/theoreticalmodel.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceEDS);

namespace ito33
{

namespace ihg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceEDS(const finance::EDS& eds) 
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *eds.GetSessionData();

  pricing::EDS edsPricing(eds);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(eds) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::EDSParams params(edsPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  GetIHGLicense().Check();

  EDSPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<EDSNumOutput> pNumoutput = pricer.Price();

  /*
  if ( pNumoutput->GetPriceSurface() )
    pNumoutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput = pNumoutput->GetOutput();
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;

}

static RegisterPriceFunction<finance::EDS> 
  regEDS(&TheoreticalModel::PriceEDS);


} // namespace ihg

} // namespace ito33

