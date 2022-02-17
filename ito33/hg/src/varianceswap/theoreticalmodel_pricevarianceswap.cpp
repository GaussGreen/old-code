/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/theoreticalmodel_pricevarianceswap.cpp
// Purpose:     Implementation of variance swap pricing using HG model
// Created:     2006/03/05
// RCS-ID:      $Id: theoreticalmodel_pricevarianceswap.cpp,v 1.15 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswaplike.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/conditionalvarianceswap.h"
#include "ito33/finance/optionvarianceswap.h"

#include "ito33/finance/error.h"

#include "ito33/pricing/varianceswap.h"
#include "ito33/pricing/varianceswapparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "hg/model.h"
#include "hg/varianceswapnumoutput.h"
#include "hg/varianceswappricer.h"
#include "hg/varianceswap_closedform.h"
#include "hg/license.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceVarianceSwap);

extern const ito33::finance::Error
  ITO33_VARIANCESWAP_LOGFORMULA_NOT_APPLICABLE,
  ITO33_VARIANCESWAP_FORWARD_STARTING_GAMMA;

namespace ito33
{

namespace hg
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

shared_ptr<finance::ModelOutput> TheoreticalModel::PriceConditionalVarianceSwap(
  const finance::ConditionalVarianceSwap& varianceSwap)
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

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *varianceSwap.GetSessionData();

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(varianceSwap) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::VarianceSwapParams 
    params(varSwap, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  finance::ComputationalFlags flags(*m_pComputFlags);
  if (    flags.GetAnalysisDate().IsValid()
       && flags.GetAnalysisDate() > varianceSwap.GetTerms()
                                    ->GetStartOfSamplingPeriod()
       && flags.GetAnalysisDate() != sessionData.GetValuationDate() )
     flags.SetAnalysisDate(Date());


  // Use analytic solver if requested and if appropriate variance swap type
  // Although Return_Actual should be possible, it's not supported currently
  // for gamma swap
  if ( flags.GetUseAnalyticSolvers() && params.IsAnalytical() 
       && (  ( varSwap.GetSwapPayoffType() != finance::SwapPayoff_Gamma )
          || ( varSwap.GetReturnType() != finance::Return_Actual ) ) )
  {
    VarianceSwapClosedForm closedForm(params, modelParams, flags);
    AutoPtr<BackwardNumOutput> pNumOutput( closedForm.Price() );
    shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );

    pOutput->SetNumOutput( make_ptr( pNumOutput.release() ) );

    return pOutput;
  }

  // TODO: add forward start gamma swap support by PDE solver
  if (    params.IsForwardStarting() 
       && varSwap.GetSwapPayoffType() == finance::SwapPayoff_Gamma )
    throw EXCEPTION(ITO33_VARIANCESWAP_FORWARD_STARTING_GAMMA);

  // Use normal PDE solver
  VarianceSwapPricer pricer(params, modelParams, flags);

  AutoPtr<VarianceSwapNumOutput> pNumOutput = pricer.Price();

  shared_ptr<ModelOutput> pOutput = pNumOutput->GetModelOutput();
  
  pOutput->SetNumOutput( shared_ptr<VarianceSwapNumOutput>(pNumOutput.release())  );
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

shared_ptr<finance::ModelOutput> 
TheoreticalModel::PriceVarianceSwapByLog
(const finance::VarianceSwap& varianceSwap)
{
  // Check that we can actually use the log formula
  shared_ptr<finance::VarianceSwapTerms> pTerms( varianceSwap.GetTerms() );

  CHECK_COND(    pTerms->GetCapMultiplier() <= 0
              && pTerms->GetUpCorridorBarrier() <= 0
              && pTerms->GetDownCorridorBarrier() <= 0
              && pTerms->GetSwapType() == finance::Swap_Variance
              && pTerms->GetReturnType() == finance::Return_Log,
              ITO33_VARIANCESWAP_LOGFORMULA_NOT_APPLICABLE );

  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *varianceSwap.GetSessionData();

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(varianceSwap) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::VarianceSwap varSwap(varianceSwap);

  pricing::VarianceSwapParams 
    params(varSwap, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  finance::ComputationalFlags flags(*m_pComputFlags);
  if (    flags.GetAnalysisDate().IsValid()
       && flags.GetAnalysisDate() > varianceSwap.GetTerms()
                                    ->GetStartOfSamplingPeriod()
       && flags.GetAnalysisDate() != sessionData.GetValuationDate() )
     flags.SetAnalysisDate(Date());


  // Just use analytical solver for this special pricing
  VarianceSwapClosedForm closedForm(params, modelParams, flags);
  AutoPtr<BackwardNumOutput> pNumOutput( closedForm.PriceByLog() );
  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );

  pOutput->SetNumOutput( make_ptr( pNumOutput.release() ) );

  return pOutput;
}

static RegisterPriceFunction<finance::VarianceSwap>
  regVarianceSwap(&TheoreticalModel::PriceVarianceSwap);

static RegisterPriceFunction<finance::GammaVarianceSwap>
  regGammaVarianceSwap(&TheoreticalModel::PriceGammaVarianceSwap);

static RegisterPriceFunction<finance::OptionVarianceSwap>
  regOptionVarianceSwap(&TheoreticalModel::PriceOptionVarianceSwap);

} // namespace hg

} // namespace ito33
