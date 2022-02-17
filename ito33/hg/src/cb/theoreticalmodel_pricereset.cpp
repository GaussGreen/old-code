/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_pricereset.cpp
// Purpose:     Implementation of reset pricing using HG model
// Created:     2006/04/17
// RCS-ID:      $Id: theoreticalmodel_pricereset.cpp,v 1.3 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/reset.h"

#include "ito33/pricing/reset.h"
#include "ito33/pricing/resetparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/bondlikeoutput.h"
#include "ito33/hg/error.h"

#include "hg/model.h"
#include "hg/cbnumoutput.h"
#include "hg/resetpricer.h"
#include "hg/license.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceReset);

extern const ito33::hg::Error ITO33_HG_EXCHANGEABLE;

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceReset(const finance::Reset& reset)
{
  double dFact;
  GetHGCBLicense().Go(&dFact);

  CHECK_COND(!reset.IsExchangeable(), ITO33_HG_EXCHANGEABLE);

  const finance::SessionData& sessionData = *reset.GetSessionData();

  pricing::Reset pricingReset(reset);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(reset) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ResetParams 
    params(pricingReset, sessionData, pNumParams, pMeshParams);

  // create HG model
  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGCBLicense().Check();

  ResetPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBNumOutput> pNumOutput = pricer.Price();

  // activate the following code for step by step debug
  // Don't do it when show convergence is enabled.
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<BondLikeOutput> pOutput( pNumOutput->GetBondLikeOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::Reset> 
  regReset(&TheoreticalModel::PriceReset);


} // namespace hg

} // namespace ito33
