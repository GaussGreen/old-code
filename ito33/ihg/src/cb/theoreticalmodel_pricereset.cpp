/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricereset.cpp
// Purpose:     Implementation of reset pricing using ihg model
// Author:      Wang
// Created:     2004/03/22
// RCS-ID:      $Id: theoreticalmodel_pricereset.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/pricing/reset.h"
#include "ito33/pricing/resetparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpricer.h"
#include "ihg/resetpricer.h"
#include "ihg/license.h"
#include "ihg/gethrforexchangeable.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/bondlikeoutput.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceReset);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceReset(const finance::Reset& reset)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *reset.GetSessionData();

  pricing::Reset pricingReset(reset);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(reset) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ResetParams 
    params(pricingReset, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  // add second hazard rate if necessary
  if ( reset.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(reset));

  GetIHGCBLicense().Check();
  
  shared_ptr<BondLikeOutput> pOutput;

  ResetPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBNumOutput> pNumOutput = pricer.Price(); 

  // activate the following code for step by step debug
  // Don't do it when show convergence is enabled.
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  pOutput = pNumOutput->GetOutput();

  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::Reset> 
  regReset(&TheoreticalModel::PriceReset);


} // namespace ihg

} // namespace ito33

