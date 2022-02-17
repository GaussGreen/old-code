/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricepercslike.cpp
// Purpose:     Implementation of PERCS-like pricing using ihg model
// Author:      Wang
// Created:     2004/08/20
// RCS-ID:      $Id: theoreticalmodel_pricepercslike.cpp,v 1.18 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/percslike.h"

#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/mandatoryparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpricer.h"
#include "ihg/license.h"
#include "ihg/gethrforexchangeable.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/bondlikeoutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPricePERCSLike);

namespace ito33
{

namespace ihg
{

// right now, defined in theoreticalmodel_pricecb.cpp
extern shared_ptr<HazardRateWithTimeComponent>
GetSecondHazardRate(const finance::ConvertibleLike& cb);


shared_ptr<finance::ModelOutput>
TheoreticalModel::PricePERCSLike
                  (const finance::PERCSLike& percsLike)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *percsLike.GetSessionData();

  pricing::Mandatory mandatory(percsLike);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(percsLike) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::MandatoryParams 
    params(mandatory, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());
  
  // add second hazard rate if necessary
  if ( percsLike.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(percsLike));

  GetIHGCBLicense().Check();
  
  shared_ptr<BondLikeOutput> pOutput;

  CBPricer pricer(params, modelParams, *m_pComputFlags);

  pOutput = pricer.Price()->GetOutput();
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::PERCSLike> 
  regPERCSLike(&TheoreticalModel::PricePERCSLike);


} // namespace ihg

} // namespace ito33

