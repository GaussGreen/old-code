/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/theoriticalmodel_pricecboption.cpp
// Purpose:     Implementation of a cb option pricing using HG model
// Created:     2006/01/19
// RCS-ID:      $Id: theoreticalmodel_pricecboption.cpp,v 1.5 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/bondlike/cboption.h"

#include "ito33/pricing/cboption.h"
#include "ito33/pricing/cboptionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/error.h"

#include "hg/model.h"
#include "hg/cboptionnumoutput.h"
#include "hg/cboptionpricer.h"
#include "hg/license.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/cboptionoutput.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_THIS_MODULE(HGPriceCBOption);

extern const ito33::hg::Error ITO33_HG_EXCHANGEABLE;

namespace ito33
{

namespace hg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceCBOption(const finance::CBOption& cboption)
{
  double dFact;
  GetHGCBLicense().Go(&dFact);

  const finance::ConvertibleBond& cb(*cboption.GetConvertibleBond());

  CHECK_COND(!cb.IsExchangeable(), ITO33_HG_EXCHANGEABLE);
  
  const finance::SessionData& sessionData = *cboption.GetSessionData();

  pricing::CBOption pricingCBOption(cboption);

  // Unlike most instruments, we use the maturity date of the underlying
  // convertible bond to make the num params. This maturity is always after
  // the maturity date of the convertible bond option.
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
  // create HG model
  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGCBLicense().Check();

  CBOptionPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBOptionNumOutput> pNumOutput = pricer.Price(); 

  // activate the following code for step by step debug
  // Don't do it when show convergence is enabled.
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<CBOptionOutput> 
    pCBOptionOutput( pNumOutput->GetCBOptionOutput() );
  
  // cb option output
  pCBOptionOutput->SetPrice(pCBOptionOutput->GetPrice() * dFact);
  
  pCBOptionOutput->SetNumOutput( make_ptr( pNumOutput.release() ) );

  return pCBOptionOutput;
}

static RegisterPriceFunction<finance::CBOption> 
  regCBOption(&TheoreticalModel::PriceCBOption);

} // namespace hg

} // namespace ito33
