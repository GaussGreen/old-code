/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_priceattachedwarrantcb.cpp
// Purpose:     Implementation of attached warrant pricing using ihg model
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: theoreticalmodel_priceattachedwarrantcb.cpp,v 1.9 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/pricing/attachedwarrantcb.h"
#include "ito33/pricing/attachedwarrantcbparams.h"

#include "ito33/ihg/bondlikeoutput.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpricer.h"
#include "ihg/attachedwarrantcbpricer.h"
#include "ihg/license.h"
#include "ihg/gethrforexchangeable.h"

#include "ito33/ihg/theoreticalmodel.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceAttachedWarrantConvertibleBond);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceAttachedWarrantConvertibleBond(
    const finance::AttachedWarrantConvertibleBond& warrant)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *warrant.GetSessionData();

  pricing::AttachedWarrantConvertibleBond pricingWarrant(warrant);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(warrant) );
             
  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::AttachedWarrantConvertibleBondParams
      params(pricingWarrant, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  // add second hazard rate if necessary
  if ( warrant.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(warrant));

  GetIHGCBLicense().Check();
  
  shared_ptr<BondLikeOutput> pOutput;

  AttachedWarrantConvertibleBondPricer
    pricer(params, modelParams, *m_pComputFlags);

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

static RegisterPriceFunction<finance::AttachedWarrantConvertibleBond> 
  regWarrant(&TheoreticalModel::PriceAttachedWarrantConvertibleBond);


} // namespace ihg

} // namespace ito33
