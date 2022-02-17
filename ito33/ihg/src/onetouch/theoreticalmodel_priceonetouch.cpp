/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/onetouch/theoriticalmodel_priceonetouch.cpp
// Purpose:     Implementation of OneTouch pricing using IHG model
// Created:     2006/08/11
// RCS-ID:      $Id: theoreticalmodel_priceonetouch.cpp,v 1.2 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/pricing/onetouch.h"
#include "ito33/pricing/onetouchparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/onetouchnumoutput.h"
#include "ihg/onetouchpricer.h"
#include "ihg/license.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/modeloutput.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_THIS_MODULE(IHGPriceOneTouch);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceOneTouch(const finance::OneTouch& oneTouch) 
{
  double dFact;
  
  GetIHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *oneTouch.GetSessionData();

  pricing::OneTouch oneTouchPricing(oneTouch);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(oneTouch) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::OneTouchParams 
    params(oneTouchPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  GetIHGLicense().Check();

  OneTouchPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<OneTouchNumOutput> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput( pNumOutput->GetOutput() );

  return pOutput;
}

static RegisterPriceFunction<finance::OneTouch> 
  regOneTouch(&TheoreticalModel::PriceOneTouch);

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceFXOneTouch(const finance::FXOneTouch& oneTouch)
{
  return PriceOneTouch(oneTouch);
}

static RegisterPriceFunction<finance::FXOneTouch> 
  regFXOneTouch(&TheoreticalModel::PriceFXOneTouch);

} // namespace hg

} // namespace ito33
