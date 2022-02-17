/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/theoreticalmodel_pricevarianceswaption.cpp
// Purpose:     Implementation of variance swap pricing using HG model
// Created:     2006/07/21
// RCS-ID:      $Id: theoreticalmodel_pricevarianceswaption.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/densematrix.h"
#include "ito33/numeric/interpolationmatrix.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/varianceswaption.h"

#include "ito33/pricing/contract.h"
#include "ito33/pricing/params.h"

#include "hg/numoutputtimeonly.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceVarianceSwaption);


namespace ito33
{

namespace hg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceVarianceSwaption
(const finance::VarianceSwaption& varianceSwaption)
{
  // Create a variance swap using the terms and an arbitrary vol strike
  double dVolStrikeTmp = 0.1;

  shared_ptr<finance::VarianceSwapTerms> pTerms( varianceSwaption.GetTerms() );

  finance::VarianceSwap varianceSwapTmp(pTerms, dVolStrikeTmp);

  const finance::SessionData& 
    sessionData( *varianceSwaption.GetSessionData() );

  // Create a new session data since we need to price the variance swap
  // at the maturity date of the swaption 
  
  // Use a temporary equity for previous share price
  shared_ptr<finance::Equity>
    pEquityTmp( new finance::Equity( *sessionData.GetEquity() ) );

  pEquityTmp->SetPreviousSharePrice( sessionData.GetSpotSharePrice() );

  shared_ptr<finance::SessionData>
    pSessionDataTmp( new finance::SessionData(sessionData) );

  // Use the maturity date of the swaption as valuation date
  pSessionDataTmp->SetValuationDate( varianceSwaption.GetMaturityDate() );

  // And the new equity
  pSessionDataTmp->SetEquity(pEquityTmp);

  // Define the session data
  varianceSwapTmp.SetSessionData(pSessionDataTmp);

  // Need analysis data at the maturity date of the swaption
  shared_ptr<finance::ComputationalFlags>
    pFlags( new finance::ComputationalFlags );

  // We should probably respect as well the internal flags, but leave it for
  // now, waiting for a global review on internal flags
  pFlags->SetAnalysisDate( varianceSwaption.GetMaturityDate() );

  varianceSwapTmp.SetComputationalFlags(pFlags);

  // Use a cloned model to avoid possible flags problem
  shared_ptr<finance::TheoreticalModel> pModel( Clone() );

  // Price the variance swap, use default flags
  shared_ptr<finance::ModelOutput> pOutputTmp( pModel->Compute(varianceSwapTmp) );

  // The yield curve, for discounting
  shared_ptr<finance::YieldCurve> pYC( pSessionDataTmp->GetYieldCurve() );

  // Code below is similar to ComputeImpliedVolStrike but we need here
  // the implied vol strike at each regime.
  double dDisCount = pYC->GetForwardDiscountFactor
                         ( GetDoubleFrom( pSessionDataTmp->GetValuationDate() ),
                           GetDoubleFrom( varianceSwapTmp.GetMaturityDate() ) );

  const double dSlope = - dDisCount;
  const double dExpectedPrice = 0.;

  shared_ptr<BackwardNumOutput> 
    pNumOutputTmp( static_pointer_cast<BackwardNumOutput> 
                   ( pOutputTmp->GetNumOutput() ) );

  const std::vector<double>&
    pdSTmp( pNumOutputTmp->GetSpotsAtAnalysisDate() ),
    pdPricesTmp( pNumOutputTmp->GetPricesAtAnalysisDate() );

  const size_t nNbS = pdSTmp.size();

  const double dSpot = sessionData.GetSpotSharePrice();

  AutoPtr<numeric::InterpolationMatrix>
    pMatrix( new numeric::InterpolationMatrix(&dSpot, 1, &pdSTmp[0], nNbS, 1) );
  
  // Interpolate for price at the spot
  size_t nNbRegimes = m_pUnderlyingProcess->GetNbRegimes();

  Array<double> pdTmp(nNbRegimes);
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pMatrix->ProductMatrixVector(&pdPricesTmp[0] + nIdxR1 * nNbS,
                                 &pdTmp[nIdxR1]);

  if ( pTerms->GetSwapType() == finance::Swap_Variance )
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      double dVolStrike = (dExpectedPrice - pdTmp[nIdxR1]) / dSlope;

      dVolStrike += dVolStrikeTmp * dVolStrikeTmp;
    
      pdTmp[nIdxR1] = dVolStrike > 0 ? sqrt(dVolStrike) : 0;
    }
  else
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdTmp[nIdxR1] = (dExpectedPrice - pdTmp[nIdxR1]) / dSlope
                    + dVolStrikeTmp;
    }

  // Payoff
  const double dStrike = varianceSwaption.GetStrike();
  const double dSqrStrike = dStrike * dStrike;
  const double dV0 = m_pUnderlyingProcess->GetPostDefaultVolatility();
  const double dSqrV0 = dV0 * dV0;
  double dDefaultValue;
  if ( varianceSwaption.GetOptionType() == finance::Option_Call )
  {
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdTmp[nIdxR1] > dStrike 
         ? dDisCount * (pdTmp[nIdxR1] * pdTmp[nIdxR1] - dSqrStrike) : 0;

    dDefaultValue = dV0 > dStrike ? dDisCount * (dSqrV0 - dSqrStrike) : 0;
  }
  else
  {
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = dStrike > pdTmp[nIdxR1] 
         ? dDisCount * (dSqrStrike - pdTmp[nIdxR1] * pdTmp[nIdxR1]) : 0;

    dDefaultValue = dStrike > dV0 ? dDisCount * (dSqrStrike - dSqrV0) : 0;
  }

  // Now solve the time only equation from real valuation date to maturity 
  // date of the swaption
  double dValuationTime = GetDoubleFrom( sessionData.GetValuationDate() );
  double dMaturityTime = GetDoubleFrom( varianceSwaption.GetMaturityDate() );
  
  numeric::DenseMatrix
    transitionProba(m_pUnderlyingProcess->ComputeRegimeTransitionProba
                                          ( dMaturityTime - dValuationTime) );

  double dNewDiscount = pYC->GetForwardDiscountFactor
                             ( dValuationTime, dMaturityTime );

  std::vector<double> pdPrices(nNbRegimes);
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    pdPrices[nIdxR1] = transitionProba[nIdxR1][nNbRegimes] * dDefaultValue;
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
      pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];

    pdPrices[nIdxR1] *= dNewDiscount;
  }

  // Use base contract and params class, nothing specific to varianceswaption
  // is actually required
  pricing::Contract varswaption(varianceSwaption);
    
  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(varianceSwaption) );
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::Params params(varswaption, sessionData, pNumParams, pMeshParams);
  if ( m_pComputFlags->GetAnalysisDate().IsValid() )
    params.SetAnalysisTime( GetDoubleFrom(m_pComputFlags->GetAnalysisDate()) );

  shared_ptr<NumOutputTimeOnly> pNumOutput( new NumOutputTimeOnly(params) );

  pNumOutput->SetValueAfterDefault(dDefaultValue);

  pNumOutput->SetPrices(pdPrices);
  
  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );

  pOutput->SetNumOutput( pNumOutput );

  return pOutput;
}

static RegisterPriceFunction<finance::VarianceSwaption>
  regVarianceSwaption(&TheoreticalModel::PriceVarianceSwaption);


} // namespace hg

} // namespace ito33
