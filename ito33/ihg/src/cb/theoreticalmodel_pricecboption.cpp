/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricecboption.cpp
// Purpose:     Implementation of a cb option pricing using ihg model
// Author:      Nabil
// Created:     2005/06/23
// RCS-ID:      $Id: theoreticalmodel_pricecboption.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/bondlike/cboption.h"

#include "ito33/pricing/cboption.h"
#include "ito33/pricing/cboptionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/model.h"
#include "ihg/cboptionnumoutput.h"
#include "ihg/cboptionpricer.h"
#include "ihg/license.h"
#include "ihg/gethrforexchangeable.h"

#include "ito33/ihg/theoreticalmodel.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceCBOption);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceCBOption(const finance::CBOption& cboption)
{
  double dFact;
  GetIHGCBLicense().Go(&dFact);

  CheckAll();

  const finance::ConvertibleBond& cb = *cboption.GetConvertibleBond(); 
  
  const finance::SessionData& sessionData = *cboption.GetSessionData();

  pricing::CBOption pricingCBOption(cboption);

  // Remark: Here, the maturity of the CB is used and not the one of the
  //         CB Option. This is because the maturity of the CB Option could be
  //         before the maturity of the CB and, also in this case, we have to
  //         give the price of the CB at the CB Option pricer as soon as we 
  //         enter in the CB Option window (in BACKWARD). 
  //         So, we have to price the CB from it maturity even if the 
  //         CB option matrurity is before.
  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                    (
                      *m_pQualityControl,
                      GetDoubleFrom( cb.GetMaturityDate() ) 
                    - GetDoubleFrom( sessionData.GetValuationDate() )
                    )
              );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CBOptionParams 
    params(pricingCBOption, sessionData, pNumParams, pMeshParams);

  //___________________________________________________________________________
  // create ihg model
  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  // add second hazard rate if necessary
  if ( cb.IsExchangeable() )
    modelParams.SetHazardRateOfDerivative(GetHRForExchangeable(cb));

  GetIHGCBLicense().Check();

  CBOptionPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBOptionNumOutput> pNumOutput = pricer.Price(); 

  // activate the following code for step by step debug
  // Don't do it when show convergence is enabled.
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ihg::CBOptionOutput> 
    pCBOptionOutput = pNumOutput->GetOutput();
  
  // cb option output
  pCBOptionOutput->SetPrice(pCBOptionOutput->GetPrice() * dFact);
  
  return pCBOptionOutput;
}

static RegisterPriceFunction<finance::CBOption> 
  regCBOption(&TheoreticalModel::PriceCBOption);


} // namespace ihg

} // namespace ito33
