/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/logcontract/logcontract_closedform.cpp
// Purpose:     HG variance swap pricer for vanilla vs
// Created:     2006/07/06
// RCS-ID:      $Id: logcontract_closedform.cpp,v 1.1 2006/07/19 17:39:51 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/logcontractparams.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/densematrix.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/logcontract_closedform.h"

namespace ito33
{

namespace hg
{

LogContractClosedForm::LogContractClosedForm
                       (pricing::LogContractParams& params, 
                        Model& model,
                        const finance::ComputationalFlags& flags)
                      : m_params(params), 
                        m_model(model),
                        m_flags(flags)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}

AutoPtr<NumOutputTimeOnly> LogContractClosedForm::Price()
{
  // Get the model parameters
  const size_t nNbRegimes = m_model.GetNbRegimes();
  const std::vector<double>& vols( m_model.GetVolatilities() );

  // Compute the model terms
  Array<double> sumVolsHat(nNbRegimes);
   
  Jumps::const_iterator jump;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    sumVolsHat[nIdxR1] = vols[nIdxR1] * vols[nIdxR1];

    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
      for (jump = jumps.begin(); jump != jumps.end(); ++jump)
      {
        double dIntensity = jump->GetIntensity();
        double dAmplitude = jump->GetAmplitude();
        double dx = log(1 + dAmplitude);

        sumVolsHat[nIdxR1] += 2 * dIntensity * (dAmplitude - dx);
      }
    }
  }

  pricing::LogContract& logContract(m_params.GetLogContract());

  double dValuationTime     = m_params.GetValuationTime();
  double dStartSamplingTime = logContract.GetT0();
  double dStopSamplingTime  = logContract.GetMaturityTime();
  double dS0 = logContract.GetS0();

  if ( numeric::IsBefore(dStartSamplingTime, dValuationTime) )
    dStartSamplingTime = dValuationTime;

  // The time step that can be used to control the accuracy of the formula
  // For now, a value similar to the one in variance swap will be used
  // We may use in the future the parameters in NumParams
  double dDeltaT = 1. / 252;

  size_t nNbSampling = int( (dStopSamplingTime - dStartSamplingTime) / dDeltaT );
  if ( nNbSampling < 1 )
    nNbSampling = 1;

  // Recompute the time step
  dDeltaT = (dStopSamplingTime - dStartSamplingTime) / nNbSampling;

  // Terms that doesn't depend on the sampling period
  Array<double> dailyReturn(nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    double dSum = 0;
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      if ( nIdxR2 != nIdxR1 )
      {
        const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
        
        for (jump = jumps.begin(); jump != jumps.end(); ++jump)
          dSum += jump->GetIntensity() 
                * (sumVolsHat[nIdxR2] - sumVolsHat[nIdxR1]);        
      }
    }
    
    dailyReturn[nIdxR1] = dDeltaT 
                        * (- 0.5 * sumVolsHat[nIdxR1] - 0.25 * dDeltaT * dSum);
  }
    
  // Terms concerning the yield curves
  Array<double> pdTimes(nNbSampling + 1);
  for (size_t nIdxT = 0; nIdxT < nNbSampling + 1; nIdxT++)
    pdTimes[nIdxT] = dStartSamplingTime + nIdxT * dDeltaT;

  // The discount, last value should not be used
  Array<double> pdRates1(nNbSampling + 1);
  m_params.GetYieldCurve()->GetForwardDiscountFactor
           ( pdTimes.Get(), pdRates1.Get(), nNbSampling + 1 );

  Array<double> pdRates2(nNbSampling + 1);
  m_params.GetForeignCurve()->GetForwardDiscountFactor
           ( pdTimes.Get(), pdRates2.Get(), nNbSampling + 1 );
  
  // pdRates1 = e^{\int(r - r_f)}
  // pdRates2 = \int(r - r_f);
  for (size_t nIdxT = 0; nIdxT < nNbSampling; nIdxT++)
  { 
    pdRates1[nIdxT] = pdRates2[nIdxT] / pdRates1[nIdxT];
    pdRates2[nIdxT] = log(pdRates1[nIdxT]);
  }

  // Compute the price at start of sampling period
  numeric::DenseMatrix
    transitionProba(m_model.ComputeRegimeTransitionProba(dDeltaT));

  std::vector<double> pdPrices(nNbRegimes);
  Array<double> pdTmp(nNbRegimes);

  for (size_t nIdxT = nNbSampling - 1; nIdxT > 0; nIdxT--)
  {
    // E_nIdxT
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdPrices[nIdxR1]
                    + dailyReturn[nIdxR1] + pdRates2[nIdxT];

    // Multiply by the transition proba matrix
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = 0.;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    }
  }

  // The first period, will be different if not forward starting
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] += dailyReturn[nIdxR1] + pdRates2[0];

  // if not forward starting
  if ( numeric::IsEqualOrAfter(dValuationTime, dStartSamplingTime) )
  {
    double dLog = log( m_params.GetSpotSharePrice() / dS0 );

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdPrices[nIdxR1] += dLog;
  }

  double dDisCount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                              ( dStartSamplingTime, dStopSamplingTime );
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] = dDisCount * pdPrices[nIdxR1];

  // Now look at if the log contract is forward starting
  if ( numeric::IsBefore(dValuationTime, dStartSamplingTime) )
  {
    numeric::DenseMatrix 
      transitionProba(m_model.ComputeRegimeTransitionProba
                              (dStartSamplingTime - dValuationTime));

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdPrices[nIdxR1];

    double dNewDiscount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                              ( dValuationTime, dStartSamplingTime );

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
      
      pdPrices[nIdxR1] *= dNewDiscount;
    }
  }

  AutoPtr<NumOutputTimeOnly> pNumOutput(new NumOutputTimeOnly(m_params));

  pNumOutput->SetPrices(pdPrices);

  return pNumOutput;
}

} // namespace hg

} // namespace ito33
