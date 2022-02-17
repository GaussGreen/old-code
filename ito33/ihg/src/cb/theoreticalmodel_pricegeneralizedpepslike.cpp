/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricegeneralizedpepslike.cpp
// Purpose:     Implementation of generalized PEPS-like pricing using ihg model
// Created:     2004/08/20
// RCS-ID:      $Id: theoreticalmodel_pricegeneralizedpepslike.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"

#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/mandatoryparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpricer.h"
#include "ihg/license.h"
#include "ihg/gethrforexchangeable.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/bondlikeoutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceGeneralizedPEPSLike);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput> 
TheoreticalModel::PriceGeneralizedPEPSLike
(const finance::GeneralizedPEPSLike& pepsLike)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *pepsLike.GetSessionData();

  pricing::Mandatory mandatory(pepsLike);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(pepsLike) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::MandatoryParams 
    params(mandatory, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());
  
  // add second hazard rate if necessary
  if ( pepsLike.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(pepsLike));

  GetIHGCBLicense().Check();
  
  shared_ptr<BondLikeOutput> pOutput;

  CBPricer pricer(params, modelParams, *m_pComputFlags);

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

static RegisterPriceFunction<finance::GeneralizedPEPSLike> 
  regGeneralizedPEPSLike(&TheoreticalModel::PriceGeneralizedPEPSLike);


} // namespace ihg

} // namespace ito33
