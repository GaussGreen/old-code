/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricebond.cpp
// Purpose:     Implementation of bond pricing using ihg model
// Author:      Wang
// Created:     2004/10/06
// RCS-ID:      $Id: theoreticalmodel_pricebond.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bond.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpricer.h"
#include "ihg/license.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/bondlikeoutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceBond);

namespace ito33
{

namespace ihg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceBond(const finance::Bond& bond)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *bond.GetSessionData();

  pricing::CB pricingCB(bond);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(bond) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CBParams 
    params(pricingCB, sessionData, pNumParams, pMeshParams);


  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  GetIHGCBLicense().Check();
  
  shared_ptr<BondLikeOutput> pOutput;

 // CBPricer pricer(params, modelParams, *m_pComputFlags);
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

static RegisterPriceFunction<finance::Bond> 
  regBond(&TheoreticalModel::PriceBond);


} // namespace ihg

} // namespace ito33

