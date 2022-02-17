/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoreticalmodel_pricevarianceswap.cpp
// Purpose:     Implementation of variance swap pricing using ihg model
// Created:     2006/02/21
// RCS-ID:      $Id: theoreticalmodel_pricevarianceswap.cpp,v 1.9 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswaplike.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/optionvarianceswap.h"

#include "ito33/finance/error.h"

#include "ito33/pricing/varianceswap.h"
#include "ito33/pricing/varianceswapparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/varianceswappricer.h"
#include "ihg/varianceswapnumoutput.h"
#include "ihg/license.h"

#include "ito33/ihg/modeloutput.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceVarianceSwap);

extern const ito33::finance::Error
  ITO33_VARIANCESWAP_FORWARD_STARTING_GAMMA;

namespace ito33
{

namespace ihg
{

  shared_ptr<finance::ModelOutput> TheoreticalModel::PriceVarianceSwap(
  const finance::VarianceSwap& varianceSwap)
{
  pricing::VarianceSwap varSwap(varianceSwap);

  return DoPriceVarianceSwapLike(varianceSwap, varSwap);
}

shared_ptr<finance::ModelOutput> TheoreticalModel::PriceGammaVarianceSwap(
  const finance::GammaVarianceSwap& varianceSwap)
{
  pricing::VarianceSwap varSwap(varianceSwap);

  return DoPriceVarianceSwapLike(varianceSwap, varSwap);
}

shared_ptr<finance::ModelOutput> TheoreticalModel::PriceOptionVarianceSwap(
  const finance::OptionVarianceSwap& varianceSwap)
{
  pricing::VarianceSwap varSwap(varianceSwap);

  return DoPriceVarianceSwapLike(varianceSwap, varSwap);
}


shared_ptr<finance::ModelOutput> TheoreticalModel::DoPriceVarianceSwapLike(
  const finance::VarianceSwapLike& varianceSwap,
  pricing::VarianceSwap& varSwap)
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *varianceSwap.GetSessionData();

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(varianceSwap) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::VarianceSwapParams 
    params(varSwap, sessionData, pNumParams, pMeshParams);
  
  // TODO: add forward start gamma swap support by PDE solver
  if (    params.IsForwardStarting() 
       && varSwap.GetSwapPayoffType() == finance::SwapPayoff_Gamma )
    throw EXCEPTION(ITO33_VARIANCESWAP_FORWARD_STARTING_GAMMA);

  Model model( GetVolatility(), GetVolatilityForMesh(), GetHazardRate() );
  model.SetPostDefaultVolatility( 
    m_pUnderlyingProcess->GetPostDefaultVolatility() );
  
  GetIHGLicense().Check();

  finance::ComputationalFlags flags(*m_pComputFlags);
  if (    flags.GetAnalysisDate().IsValid()
       && flags.GetAnalysisDate() > varianceSwap.GetTerms()
                                    ->GetStartOfSamplingPeriod()
       && flags.GetAnalysisDate() != sessionData.GetValuationDate() )
    flags.SetAnalysisDate(Date());

  VarianceSwapPricer pricer(params, model, flags);

  shared_ptr<ihg::ModelOutput> pOutput = pricer.Price()->GetModelOutput();

  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::VarianceSwap>
  regVarianceSwap(&ihg::TheoreticalModel::PriceVarianceSwap);

static RegisterPriceFunction<finance::GammaVarianceSwap>
  regGammaVarianceSwap(&ihg::TheoreticalModel::PriceGammaVarianceSwap);

static RegisterPriceFunction<finance::OptionVarianceSwap>
  regOptionVarianceSwap(&ihg::TheoreticalModel::PriceOptionVarianceSwap);



} // namespace ihg

} // namespace ito33

