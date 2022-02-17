/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/varianceswap_closedform.cpp
// Purpose:     HG variance swap pricer for vanilla vs
// Created:     2006/07/06
// RCS-ID:      $Id: varianceswap_closedform.cpp,v 1.16 2006/08/21 11:14:44 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"

#include "ito33/numeric/densematrix.h"

#include "hg/model.h"
#include "hg/numoutputanalytical.h"
#include "hg/numoutputtimeonly.h"
#include "hg/varianceswap_closedform.h"

namespace ito33
{

namespace hg
{

VarianceSwapClosedForm::VarianceSwapClosedForm
                        (pricing::VarianceSwapParams& params, 
                         Model& model,
                         const finance::ComputationalFlags& flags)
                       : m_params(params), 
                         m_model(model),
                         m_flags(flags)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}

AutoPtr<BackwardNumOutput> VarianceSwapClosedForm::Price()
{
  pricing::VarianceSwap& vs( m_params.GetVarianceSwap() );

  if ( vs.GetSwapPayoffType() == finance::SwapPayoff_Gamma )
    return PriceGammaSwap();

  // Get the model parameters
  const size_t nNbRegimes = m_model.GetNbRegimes();
  const std::vector<double>& vols( m_model.GetVolatilities() );
  const double dVolInDefault = m_model.GetPostDefaultVolatility();

  // Compute the model terms
  Array<double>
    sumVols(nNbRegimes), sumVolsHat(nNbRegimes), sumVolsTilde(nNbRegimes);
   
  Jumps::const_iterator jump;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    sumVols[nIdxR1] = sumVolsHat[nIdxR1] = sumVolsTilde[nIdxR1] 
                    = vols[nIdxR1] * vols[nIdxR1];

    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
      for (jump = jumps.begin(); jump != jumps.end(); ++jump)
      {
        double dIntensity = jump->GetIntensity();
        double dAmplitude = jump->GetAmplitude();
        double dx = log(1 + dAmplitude);

        sumVols[nIdxR1] += dIntensity * dAmplitude * dAmplitude;
        sumVolsTilde[nIdxR1] += dIntensity * dx * dx;
        sumVolsHat[nIdxR1] += 2 * dIntensity * (dAmplitude - dx);
      }
    }
  }
  
  double dStartSamplingTime = vs.GetStartTimeOfSamplingPeriod();
  double dStopSamplingTime = vs.GetMaturityTime();
  size_t nNbSampling = vs.GetNbSamplingReturns();
  size_t nNbSamplesUsed = vs.GetNbSamplesUsed();
  double dA = vs.GetAnnualReturnFrequency();
  finance::ReturnType returnType = vs.GetReturnType();

  double dStrike = vs.GetVolatilityStrike() * vs.GetVolatilityStrike();
  bool bForwardStarting = m_params.IsForwardStarting();
  double dValuationTime = m_params.GetValuationTime();

  if ( !bForwardStarting )
  {
    dStartSamplingTime = dValuationTime;
    nNbSampling -= nNbSamplesUsed;
  }

  // The length of a sampling period
  double dDeltaT = (dStopSamplingTime - dStartSamplingTime) / nNbSampling;

  // Terms that doesn't depend on the sampling period
  
  // Only used by log since 0 for actual 
  Array<double> dailyReturn(nNbRegimes);

  // coefficients, multiplied by log(S/S_i) or (S/S_i)^2
  Array<double> coefficient(nNbRegimes);

  if ( returnType == finance::Return_Log )
  {
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      double dSum1 = 0, dSum2 = 0;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
      {
        if ( nIdxR2 != nIdxR1 )
        {
          const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
          
          for (jump = jumps.begin(); jump != jumps.end(); ++jump)
          {
            double dIntensity = jump->GetIntensity();
            double dAmplitude = jump->GetAmplitude();
            double dx = log(1 + dAmplitude);

            dSum1 += dIntensity
                   * (  sumVolsTilde[nIdxR2] - sumVolsTilde[nIdxR1]
                      - dx * (sumVolsHat[nIdxR2] - sumVolsHat[nIdxR1]) ); 
    
            dSum2 += dIntensity * (sumVolsHat[nIdxR2] - sumVolsHat[nIdxR1]);
          }         
        }
      }
      
      dailyReturn[nIdxR1] = dDeltaT
                          * ( sumVolsTilde[nIdxR1] + 0.5 * dDeltaT * dSum1 );

      coefficient[nIdxR1] = - dDeltaT 
                          * ( sumVolsHat[nIdxR1] + 0.5 * dDeltaT * dSum2 );
    }
  }
  else
  {
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      double dSum = 0.;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
      {     
        if ( nIdxR2 != nIdxR1 )
        {
          const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
          for (jump = jumps.begin(); jump != jumps.end(); ++jump)
          {
            double dIntensity = jump->GetIntensity();
            double dAmplitude = jump->GetAmplitude();
            dSum += dIntensity * (1. + dAmplitude) * (1. + dAmplitude)
                  * (sumVols[nIdxR2] - sumVols[nIdxR1]);
          }
        }
      }
      
      coefficient[nIdxR1] = dDeltaT * (sumVols[nIdxR1] + 0.5 * dDeltaT * dSum);
    }
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
    transitionProba( m_model.ComputeRegimeTransitionProba(dDeltaT) );
  
  double dDefaultTerm = dVolInDefault * dVolInDefault / dA;

  std::vector<double> pdPrices(nNbRegimes);
  Array<double> pdTmp(nNbRegimes);

  for (size_t nIdxT = nNbSampling - 1; nIdxT > 0; nIdxT--)
  {
    // E_nIdxT - dDefaultTerm
    if ( returnType == finance::Return_Log )
    {
      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {
        double dTmp = pdRates2[nIdxT] - 0.5 * dDeltaT * sumVolsHat[nIdxR1];
        pdTmp[nIdxR1] = pdPrices[nIdxR1]
                      + dailyReturn[nIdxR1] + dTmp * dTmp - dDefaultTerm;
      }
    }
    else
    {
      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {
        double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVols[nIdxR1];
        pdTmp[nIdxR1] = pdPrices[nIdxR1]
                      + coefficient[nIdxR1]
                      + 2. * (1 - pdRates1[nIdxT] + pdRates2[nIdxT])
                      + 2. * dTmp * dTmp - dDefaultTerm;
      }
    }

    // Multiply by the transition proba matrix
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = 0.;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    }
  }

  double dX0 = 0.;
  if ( !bForwardStarting )
    dX0 = vs.GetCurrentVolatility() * vs.GetCurrentVolatility()
        * nNbSamplesUsed / dA;

  double dTmp = dX0 + dDefaultTerm * nNbSampling;
  double dRecoveryValue = vs.GetPayoffValue(dTmp / vs.GetNbSamplingReturns());

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] += dX0 + dDefaultTerm * (nNbSampling - 1);

  // Scalar value for non forward starting
  double dDelta = 0, dGamma = 0;

  // Analysis data variables for non forward starting
  bool bRequestAnalysisData = false;
  size_t nNbS = 1;
  std::vector<double> pdLogS, pdS, pdSPrices;

  // The first period, will be different if not forward starting
  if ( bForwardStarting )
  {
    const size_t nIdxT = 0;

    if ( returnType == finance::Return_Log )   
    {
      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {
        double dTmp = pdRates2[nIdxT] - 0.5 * dDeltaT * sumVolsHat[nIdxR1]; 
      
        pdPrices[nIdxR1] += dailyReturn[nIdxR1] + dTmp * dTmp;
      }  
    }
    else
    {
      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {
        double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVols[nIdxR1];
        
        pdPrices[nIdxR1] += coefficient[nIdxR1]
                          + 2. * (1 - pdRates1[nIdxT] + pdRates2[nIdxT])
                          + 2. * dTmp * dTmp;
      }
    }
  }
  else
  {   
    // To be consistent with our convention, we'll just take
    // \delta t' = \delta t

    // Check if analysis data can be supported
    if ( numeric::AreTimesEqual( dValuationTime, m_params.GetAnalysisTime() ) )
    {
      bRequestAnalysisData = true;
      pdLogS = m_params.GenerateSpaceMesh(m_model);
      nNbS = pdLogS.size();

      pdS.resize(nNbS);
      for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
        pdS[nIdxS] = exp(pdLogS[nIdxS]);

      pdSPrices.resize(nNbRegimes * nNbS);
    }

    const size_t nIdxT = 0;

    double dSpot = m_params.GetSpotSharePrice();
    if ( returnType == finance::Return_Log )   
    {  
      double dLogPreSpot = log( vs.GetPreviousSharePrice() );
      double dLog0 = log(dSpot) - dLogPreSpot;

      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {    
        double dTmp = pdRates2[nIdxT] - 0.5 * dDeltaT * sumVolsHat[nIdxR1];
        double dTmp1 = coefficient[nIdxR1] + 2 * pdRates2[nIdxT];

        if ( bRequestAnalysisData )
        {
          double *pdTmp = &pdSPrices[0] + nIdxR1 * nNbS;
          for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
          {
            double dLog = pdLogS[nIdxS] - dLogPreSpot;
            pdTmp[nIdxS] = pdPrices[nIdxR1] + dLog * dLog + dLog * dTmp1
                        + dailyReturn[nIdxR1] + dTmp * dTmp;       
          }
        }

        pdPrices[nIdxR1] += dLog0 * dLog0 + dLog0 * dTmp1
                          + dailyReturn[nIdxR1] + dTmp * dTmp;

        if ( nIdxR1 == 0 )
        {
          double dInvSpot = 1. / dSpot;
          dDelta = dInvSpot * ( 2. * dLog0 + dTmp1 );
          dGamma = - dInvSpot * dDelta + 2. * dInvSpot * dInvSpot;
        }
      }
    }
    else
    {
      double dInvPreSpot = 1. / vs.GetPreviousSharePrice();

      for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      {
        double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVols[nIdxR1];
        double dTmp2 = 1. + coefficient[nIdxR1] 
                     + 2. * pdRates2[nIdxT] + 2. * dTmp * dTmp;

        if ( bRequestAnalysisData )
        {
          double *pdTmp = &pdSPrices[0] + nIdxR1 * nNbS;
          for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
          {
            double dQ = pdS[nIdxS] * dInvPreSpot;
            pdTmp[nIdxS] = pdPrices[nIdxR1]
                         + 1. - 2. * dQ * pdRates1[nIdxT] + dQ * dQ * dTmp2;    
          }
        }
        
        double dQ0 = dSpot * dInvPreSpot;
        pdPrices[nIdxR1] += 1. - 2. * dQ0 * pdRates1[nIdxT] + dQ0 * dQ0 * dTmp2;

        if ( nIdxR1 == 0 )
        {
          dDelta = - 2. * dInvPreSpot * pdRates1[nIdxT]
                 + 2. * dQ0 * dInvPreSpot * dTmp2;

          dGamma = 2. * dInvPreSpot * dInvPreSpot * dTmp2;
        }
      }
    }
  }

  double dDisCount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                              ( dStartSamplingTime, dStopSamplingTime );

  dRecoveryValue *= dDisCount;

  double dScale = dA / vs.GetNbSamplingReturns()* dDisCount;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] *= dScale;

  dDelta *= dScale;
  dGamma *= dScale;
  
  double dDisCountStrike = dStrike * dDisCount;
  
  // Now look at if the vs is forward starting
  if ( bForwardStarting )
  {
    numeric::DenseMatrix 
      transitionProba(m_model.ComputeRegimeTransitionProba
                              (dStartSamplingTime - dValuationTime));

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdPrices[nIdxR1];

    // Post default term discounted
    double dTmp = dVolInDefault * dVolInDefault * dDisCount;
    double dNewDiscount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                              ( dValuationTime, dStartSamplingTime );

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = transitionProba[nIdxR1][nNbRegimes] * dTmp;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
      pdPrices[nIdxR1] *= dNewDiscount;
    }

    dDisCountStrike *= dNewDiscount;
    dRecoveryValue *= dNewDiscount;
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] -= dDisCountStrike;

  if ( !bForwardStarting )
  {
    AutoPtr<NumOutputAnalytical> 
      pNumOutput( new NumOutputAnalytical(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrice(pdPrices[0]);
    pNumOutput->SetDelta(dDelta);
    pNumOutput->SetGamma(dGamma);

    if ( bRequestAnalysisData )
    {
      for (size_t nIdxX = 0; nIdxX < pdSPrices.size(); nIdxX++)
        pdSPrices[nIdxX] = dScale * pdSPrices[nIdxX] - dDisCountStrike;

      pNumOutput->SetFinalValues(pdS, pdSPrices);
    }

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
  else
  {
    AutoPtr<NumOutputTimeOnly> pNumOutput( new NumOutputTimeOnly(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrices(pdPrices);

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
}

AutoPtr<BackwardNumOutput> VarianceSwapClosedForm::PriceGammaSwap()
{
  // Get the model parameters
  const size_t nNbRegimes = m_model.GetNbRegimes();
  const std::vector<double>& vols( m_model.GetVolatilities() );

  // Compute the model terms
  Array<double> sumVolsTilde2(nNbRegimes), sumVolsHat2(nNbRegimes);
   
  Jumps::const_iterator jump;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    sumVolsTilde2[nIdxR1] = sumVolsHat2[nIdxR1] = vols[nIdxR1] * vols[nIdxR1];
    sumVolsHat2[nIdxR1] += 2. * m_model.GetJumpsToDefault()[nIdxR1];

    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
      for (jump = jumps.begin(); jump != jumps.end(); ++jump)
      {
        double dIntensity = jump->GetIntensity();
        double dAmplitude = jump->GetAmplitude();
        double dx = log(1 + dAmplitude);

        sumVolsTilde2[nIdxR1] += dIntensity * dx * dx * (1. + dAmplitude);
        sumVolsHat2[nIdxR1] += 2 * dIntensity
                             * (dx + dx * dAmplitude - dAmplitude);
      }
    }
  }

  pricing::VarianceSwap& vs(m_params.GetVarianceSwap());
    
  // Currently, only log return is supported
  ASSERT( vs.GetReturnType() == finance::Return_Log );
  
  double dStartSamplingTime = vs.GetStartTimeOfSamplingPeriod();
  double dStopSamplingTime = vs.GetMaturityTime();
  size_t nNbSampling = vs.GetNbSamplingReturns();
  size_t nNbSamplesUsed = vs.GetNbSamplesUsed();
  double dA = vs.GetAnnualReturnFrequency();
  double dStrike = vs.GetVolatilityStrike() * vs.GetVolatilityStrike();
  bool bForwardStarting = m_params.IsForwardStarting();
  double dValuationTime = m_params.GetValuationTime();

  if ( !bForwardStarting )
  {
    dStartSamplingTime = dValuationTime;
    nNbSampling -= nNbSamplesUsed;
  }

  // The length of a sampling period
  double dDeltaT = (dStopSamplingTime - dStartSamplingTime) / nNbSampling;

  // Terms that doesn't depend on the sampling period

  // coefficients, multiplied by S/S_i
  Array<double> coeActual(nNbRegimes);
  
  // coefficients, multiplied by S/S_i log(S/S_i)
  Array<double> coeMixed(nNbRegimes);

  // coeffiecients multiplied by S/S_i log^2(S/S_i) is just 0

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    double dSum1 = 0, dSum2 = 0;
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      if ( nIdxR2 != nIdxR1 )
      {
        const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
        
        for (jump = jumps.begin(); jump != jumps.end(); ++jump)
        {
          double dIntensity = jump->GetIntensity();
          double dAmplitude = jump->GetAmplitude();
          double dx = log(1 + dAmplitude);

          dSum1 += dIntensity * (1. + dAmplitude)
                 * (  sumVolsTilde2[nIdxR2] - sumVolsTilde2[nIdxR1]
                    + dx * (sumVolsHat2[nIdxR2] - sumVolsHat2[nIdxR1]) );
  
          dSum2 += dIntensity * (1. + dAmplitude) 
                 * (sumVolsHat2[nIdxR2] - sumVolsHat2[nIdxR1]);
        }         
      }
    }
    
    coeActual[nIdxR1] = dDeltaT
                      * ( sumVolsTilde2[nIdxR1] + 0.5 * dDeltaT * dSum1 );

    coeMixed[nIdxR1] = dDeltaT 
                     * ( sumVolsHat2[nIdxR1] + 0.5 * dDeltaT * dSum2 );
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
    transitionProba( m_model.ComputeSpotRegimeTransitionProba(dDeltaT) );

  std::vector<double> pdPrices(nNbRegimes);
  Array<double> pdTmp(nNbRegimes);

  for (size_t nIdxT = nNbSampling - 1; nIdxT > 0; nIdxT--)
  {
    // E_nIdxT
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVolsHat2[nIdxR1];
      pdTmp[nIdxR1] = pdPrices[nIdxR1]
                    + coeActual[nIdxR1]
                    + dDeltaT * pdRates2[nIdxT] * sumVolsTilde2[nIdxR1]
                    + dTmp * dTmp;
    }

    // Multiply by the transition proba matrix
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = 0.;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    }

    // Time dependent part of transition proba matrix
    // for now, ignore yield dividend
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdPrices[nIdxR1] *= pdRates1[nIdxT];
  }

  // Scalar value for non forward starting
  double dDelta = 0, dGamma = 0;

  // Analysis data variables for non forward starting
  bool bRequestAnalysisData = false;
  size_t nNbS = 1;
  std::vector<double> pdLogS, pdS, pdSPrices;

  // The first period, will be different if not forward starting
  if ( bForwardStarting )
  {
    const size_t nIdxT = 0;

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVolsHat2[nIdxR1]; 
    
      pdPrices[nIdxR1] += coeActual[nIdxR1]
                        + dDeltaT * pdRates2[nIdxT] * sumVolsTilde2[nIdxR1]
                        + dTmp * dTmp;
    }
  }
  else
  {   
    // To be consistent with our convention, we'll just take
    // \delta t' = \delta t
    double dX0 = vs.GetCurrentVolatility() * vs.GetCurrentVolatility()
               * nNbSamplesUsed / dA;

    // Check if analysis data can be supported
    if ( numeric::AreTimesEqual( dValuationTime, m_params.GetAnalysisTime() ) )
    {
      bRequestAnalysisData = true;
      pdLogS = m_params.GenerateSpaceMesh(m_model);
      nNbS = pdLogS.size();

      pdS.resize(nNbS);
      for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
        pdS[nIdxS] = exp(pdLogS[nIdxS]);

      pdSPrices.resize(nNbRegimes * nNbS);
    }

    const size_t nIdxT = 0;

    double dSpot = m_params.GetSpotSharePrice();
    double dPreSpot = vs.GetPreviousSharePrice();
    
    double dLogPreSpot = log(dPreSpot);
    double dInvPreSpot = 1. / dPreSpot;

    double dLog0 = log(dSpot) - dLogPreSpot;
    double dQ0 = dSpot * dInvPreSpot;
    
    double dCoe = dPreSpot / vs.GetStartSharePrice();

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {    
      double dTmp = pdRates2[nIdxT] + 0.5 * dDeltaT * sumVolsHat2[nIdxR1];
      double dTmpAcutal = coeActual[nIdxR1]
                        + dDeltaT * pdRates2[nIdxT] * sumVolsTilde2[nIdxR1]
                        + dTmp * dTmp;
      double dTmpMixed = coeMixed[nIdxR1]
                       + 2. * pdRates2[nIdxT]
                       + 2. * pdRates2[nIdxT] * pdRates2[nIdxT]
                       + pdRates2[nIdxT] * dDeltaT * sumVolsHat2[nIdxR1];
      double dTmpMixed2 = pdRates1[nIdxT];

      if ( bRequestAnalysisData )
      {
        double *pdTmp = &pdSPrices[0] + nIdxR1 * nNbS;
        for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
        {
          double dLog = pdLogS[nIdxS] - dLogPreSpot;
          double dQ = pdS[nIdxS] * dInvPreSpot;

          pdTmp[nIdxS] = pdPrices[nIdxR1] 
                       + dQ * dTmpAcutal
                       + dQ * dLog * dTmpMixed
                       + dQ * dLog * dLog * dTmpMixed2;       
        }
      }

      pdPrices[nIdxR1] += dQ0 * dTmpAcutal
                        + dQ0 * dLog0 * dTmpMixed
                        + dQ0 * dLog0 * dLog0 * dTmpMixed2;  

      pdPrices[nIdxR1] = dCoe * pdPrices[nIdxR1] + dX0;

      if ( nIdxR1 == 0 )
      {
        double dInvSpot = 1. / dSpot;
        dDelta = dCoe * dInvPreSpot 
               * (   dTmpAcutal 
                   + (1 + dLog0) * dTmpMixed
                   + (dLog0 * dLog0 + 2. * dLog0) * dTmpMixed2 );
        dGamma = dCoe * dInvPreSpot 
               * (   dInvSpot * dTmpMixed
                   + 2. * (dLog0 + 1) * dInvSpot * dTmpMixed2 );
      }
    }

    for (size_t nIdxX = 0; nIdxX < pdSPrices.size(); nIdxX++)
      pdSPrices[nIdxX] = dCoe * pdSPrices[nIdxX] + dX0;
  }

  double dDisCount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                              ( dStartSamplingTime, dStopSamplingTime );

  double dScale = dA / vs.GetNbSamplingReturns() * dDisCount;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] *= dScale;

  dDelta *= dScale;
  dGamma *= dScale;
  
  double dDisCountStrike = dStrike * dDisCount;
  
  // Now look at if the vs is forward starting
  if ( bForwardStarting )
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
      pdPrices[nIdxR1] = 0;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
      pdPrices[nIdxR1] *= dNewDiscount;
    }

    dDisCountStrike *= dNewDiscount;
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] -= dDisCountStrike;
  
  // We don't recovery anything when default occurs
  double dRecoveryValue = 0;

  if ( !bForwardStarting )
  {
    AutoPtr<NumOutputAnalytical> 
      pNumOutput( new NumOutputAnalytical(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrice(pdPrices[0]);
    pNumOutput->SetDelta(dDelta);
    pNumOutput->SetGamma(dGamma);

    if ( bRequestAnalysisData )
    {
      for (size_t nIdxX = 0; nIdxX < pdSPrices.size(); nIdxX++)
        pdSPrices[nIdxX] = dScale * pdSPrices[nIdxX] - dDisCountStrike;

      pNumOutput->SetFinalValues(pdS, pdSPrices);
    }

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
  else
  {
    AutoPtr<NumOutputTimeOnly> pNumOutput( new NumOutputTimeOnly(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrices(pdPrices);

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
}

AutoPtr<BackwardNumOutput> VarianceSwapClosedForm::PriceByLog()
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
  
  pricing::VarianceSwap& vs( m_params.GetVarianceSwap() );

  double dStartSamplingTime = vs.GetStartTimeOfSamplingPeriod();
  double dStopSamplingTime = vs.GetMaturityTime();
  size_t nNbSampling = vs.GetNbSamplingReturns();
  size_t nNbSamplesUsed = vs.GetNbSamplesUsed();
  double dA = vs.GetAnnualReturnFrequency();
  double dStrike = vs.GetVolatilityStrike() * vs.GetVolatilityStrike();
  bool bForwardStarting = m_params.IsForwardStarting();
  double dValuationTime = m_params.GetValuationTime();

  if ( !bForwardStarting )
  {
    dStartSamplingTime = dValuationTime;
    nNbSampling -= nNbSamplesUsed;
  }

  // The length of a sampling period
  double dDeltaT = (dStopSamplingTime - dStartSamplingTime) / nNbSampling;

  // Terms that doesn't depend on the sampling period

  // values will be added directly to return of each period
  Array<double> daylyReturn(nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    double dSum = 0;
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      if ( nIdxR2 != nIdxR1 )
      {
        const Jumps& jumps( m_model.GetJumps(nIdxR1, nIdxR2) );
        
        for (jump = jumps.begin(); jump != jumps.end(); ++jump)
        {
          double dIntensity = jump->GetIntensity();

          dSum += dIntensity * (sumVolsHat[nIdxR2] - sumVolsHat[nIdxR1]);
        }         
      }
    }
    
    daylyReturn[nIdxR1] = - 0.5 * dDeltaT
                        * ( sumVolsHat[nIdxR1] + 0.5 * dDeltaT * dSum );
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
    transitionProba( m_model.ComputeRegimeTransitionProba(dDeltaT) );

  std::vector<double> pdPrices(nNbRegimes);
  Array<double> pdTmp(nNbRegimes);

  for (size_t nIdxT = nNbSampling - 1; nIdxT > 0; nIdxT--)
  {
    // E_nIdxT
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdPrices[nIdxR1] + pdRates2[nIdxT] + daylyReturn[nIdxR1];

    // Multiply by the transition proba matrix
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = 0.;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
    }
  }

  // Scalar value for non forward starting
  double dDelta = 0, dGamma = 0;

  // Analysis data variables for non forward starting
  bool bRequestAnalysisData = false;
  size_t nNbS = 1;
  std::vector<double> pdLogS, pdS, pdSPrices;
  
  double dSumRate = 0.;
  for ( size_t nIdxT1 = 0; nIdxT1 < nNbSampling; nIdxT1++)
    dSumRate += pdRates2[nIdxT1];

  double dX0 = 0.;

  // The first period, will be different if not forward starting
  if ( bForwardStarting )
  {
    const size_t nIdxT = 0;

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdPrices[nIdxR1] += pdRates2[nIdxT] + daylyReturn[nIdxR1];
  }
  else
  {   
    // To be consistent with our convention, we'll just take
    // \delta t' = \delta t
    dX0 = vs.GetCurrentVolatility() * vs.GetCurrentVolatility()
        * nNbSamplesUsed / dA;

    // Check if analysis data can be supported
    if ( numeric::AreTimesEqual( dValuationTime, m_params.GetAnalysisTime() ) )
    {
      bRequestAnalysisData = true;
      pdLogS = m_params.GenerateSpaceMesh(m_model);
      nNbS = pdLogS.size();

      pdS.resize(nNbS);
      for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
        pdS[nIdxS] = exp(pdLogS[nIdxS]);

      pdSPrices.resize(nNbRegimes * nNbS);
    }

    const size_t nIdxT = 0;

    double dSpot = m_params.GetSpotSharePrice();
    double dPreSpot = vs.GetPreviousSharePrice();
    
    double dLogPreSpot = log(dPreSpot);

    double dLog0 = log(dSpot) - dLogPreSpot;

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      if ( bRequestAnalysisData )
      {
        double *pdTmp = &pdSPrices[0] + nIdxR1 * nNbS;
        for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
        {
          double dLog = pdLogS[nIdxS] - dLogPreSpot;

          pdTmp[nIdxS] = pdPrices[nIdxR1]
                       + dLog + pdRates2[nIdxT] + daylyReturn[nIdxR1];  
        }
      }

      pdPrices[nIdxR1] += dLog0 + pdRates2[nIdxT] + daylyReturn[nIdxR1];

      if ( nIdxR1 == 0 )
      {
        double dInvSpot = 1. / dSpot;
        dDelta = - 2. * dInvSpot;
        dGamma = 2. * dInvSpot * dInvSpot;
      }
    }

    for (size_t nIdxX = 0; nIdxX < pdSPrices.size(); nIdxX++)
      pdSPrices[nIdxX] = dX0 - 2. * pdSPrices[nIdxX] + 2. * dSumRate;
  }

  double dDisCount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                            ( dStartSamplingTime, dStopSamplingTime );

  double dScale = dA / vs.GetNbSamplingReturns() * dDisCount;
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    pdPrices[nIdxR1] = dX0 - 2. * pdPrices[nIdxR1] + 2. * dSumRate;
    pdPrices[nIdxR1] *= dScale;
  }

  dDelta *= dScale;
  dGamma *= dScale;
  
  double dDisCountStrike = dStrike * dDisCount;
  
  // Now look at if the vs is forward starting
  if ( bForwardStarting )
  {
    numeric::DenseMatrix 
      transitionProba( m_model.ComputeRegimeTransitionProba
                               ( dStartSamplingTime - dValuationTime ) );

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
      pdTmp[nIdxR1] = pdPrices[nIdxR1];

    double dNewDiscount = m_params.GetYieldCurve()->GetForwardDiscountFactor
                                   ( dValuationTime, dStartSamplingTime );

    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      pdPrices[nIdxR1] = 0;
      for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
        pdPrices[nIdxR1] += transitionProba[nIdxR1][nIdxR2] * pdTmp[nIdxR2];
      pdPrices[nIdxR1] *= dNewDiscount;
    }

    dDisCountStrike *= dNewDiscount;
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    pdPrices[nIdxR1] -= dDisCountStrike;
  
  AutoPtr<BackwardNumOutput> pNumOutputTmp;

  // We don't recovery anything when default occurs
  double dRecoveryValue = 0;

  if ( !bForwardStarting )
  {
    AutoPtr<NumOutputAnalytical> 
      pNumOutput( new NumOutputAnalytical(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrice(pdPrices[0]);
    pNumOutput->SetDelta(dDelta);
    pNumOutput->SetGamma(dGamma);

    if ( bRequestAnalysisData )
    {
      for (size_t nIdxX = 0; nIdxX < pdSPrices.size(); nIdxX++)
        pdSPrices[nIdxX] = dScale * pdSPrices[nIdxX] - dDisCountStrike;

      pNumOutput->SetFinalValues(pdS, pdSPrices);
    }

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
  else
  {
    AutoPtr<NumOutputTimeOnly> pNumOutput( new NumOutputTimeOnly(m_params) );

    pNumOutput->SetValueAfterDefault(dRecoveryValue);

    pNumOutput->SetPrices(pdPrices);

    return AutoPtr<BackwardNumOutput>( pNumOutput.release() );
  }
}

} // namespace hg

} // namespace ito33
