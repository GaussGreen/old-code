/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricecb.cpp
// Purpose:     Implementation of cb pricing using ihg model
// Author:      Wang
// Created:     2004/03/22
// RCS-ID:      $Id: theoreticalmodel_pricecb.cpp,v 1.30 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbparams.h"

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
#include "ito33/ihg/bondlikeoutput.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceCB);

extern const ito33::finance::BondError
    ITO33_EXCHANGEABLE_NO_ISSUER_DEFAULT;

namespace ito33
{

namespace ihg
{

shared_ptr<HazardRateTimeOnly>
GetSecondHazardRate(const finance::ConvertibleLike& cb);

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceConvertibleBond
              (const finance::ConvertibleBond& cb)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *cb.GetSessionData();

  pricing::CB pricingCB(cb);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(cb) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CBParams params(pricingCB, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());


  // add second hazard rate if necessary
  if ( cb.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(cb));

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

static RegisterPriceFunction<finance::ConvertibleBond> 
  regConvertibleBond(&TheoreticalModel::PriceConvertibleBond);

} // namespace ihg

} // namespace ito33
