/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/varianceswapmeshmanager.cpp
// Purpose:     variance swap mesh manager for backward PDE problems
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapmeshmanager.cpp,v 1.12 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"
#include "ito33/pricing/optionliketype.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

VarianceSwapMeshManager::VarianceSwapMeshManager(
  VarianceSwapParams& params, 
  Model& model, 
  const std::vector<double>& pdSpaceMesh)
  : BackwardMeshManagerFix(params, model),
    m_varianceSwapParams(params),
    m_bIsConditionalFixed(false)
{
  // Copy the space mesh, and construct the log mesh
  m_nNbS = pdSpaceMesh.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdS[nIdx] = pdSpaceMesh[nIdx];
    m_pdLogS[nIdx] = log( m_pdS[nIdx] );    
  }
}


void VarianceSwapMeshManager::ConstructSpaceMesh()
{
  // Check if the mesh was created in the constructor
  if (m_nNbS > 0)
    return;

  if (m_params.GetMeshParams()->GetUniformSpaceGrid() == true)
  {
    ConstructUniformSpaceMesh();
    return;
  }

  std::vector<double> vecGrid( m_params.GenerateSpaceMesh(m_model) );

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);
  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx];
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
}

void VarianceSwapMeshManager::ConstructUniformSpaceMesh()
{
  // Force the points 0.0 and the spot.  Continue until approximately
  // three times the spot.  This means about 1/3 of the points will be
  // between 0 and the spot.

  m_nNbS = m_params.GetNumParams()->GetNbSpaceSteps();
  size_t nNbTmp = m_nNbS / 3;
  double dDeltaT = m_varianceSwapParams.GetSpotSharePrice() / nNbTmp;
  
  // Construct the actual mesh
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS[0] = 0.0;
  m_pdLogS[0] = log(dDeltaT / 100.0);
  for (size_t nIdx = 1; nIdx < m_nNbS; nIdx++)
  {    
    m_pdS[nIdx] = nIdx * dDeltaT;
    m_pdLogS[nIdx] = log(m_pdS[nIdx]);
  }
}

void VarianceSwapMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  VarianceSwap &varianceSwap = m_varianceSwapParams.GetVarianceSwap();
    
  //Compute discount factor
  Array<double> pdCF(m_nNbTimes);
      
  m_varianceSwapParams.GetYieldCurve()->GetCompoundFactor(m_pdTimes.Get(), 
                                                          pdCF.Get(), 
                                                          m_nNbTimes);
 
  double dTmp = 1. / pdCF[m_nNbTimes - 1];

  // Get the average squared return of the path currently considered
  // Related to X in Philippe's documentation
  double dAvgSqrReturns = m_varianceSwapParams.GetAvgSqrReturn();
  double dTini          = varianceSwap.GetStartTimeOfSamplingPeriod();
  size_t N              = varianceSwap.GetNbSamplingReturns();
  double dDefaultVar    = m_model.GetPostDefaultVolatility() *
                          m_model.GetPostDefaultVolatility() ;

  double dInvAnnualReturnFreq 
    = 1./ double( varianceSwap.GetAnnualReturnFrequency() );

  for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
  {
    double dValue = 0.0;

    if ( m_bIsConditionalFixed )
      dValue = m_varianceSwapParams.ComputeConditionalPayoff();
    else if ( numeric::IsEqualOrBefore( m_pdTimes[nIdxT] , dTini ) )
      dValue = varianceSwap.GetPayoffValue(dInvAnnualReturnFreq * dDefaultVar);
    else
    {
      bool bIsTimeOnEvent = false;
      size_t nObs = GetObservationNumber(m_pdTimes[nIdxT], bIsTimeOnEvent);

      double dX = 0.0;
      if ( bIsTimeOnEvent )
        dX = (nObs - 1) * dAvgSqrReturns 
           + ( N - nObs ) * dInvAnnualReturnFreq * dDefaultVar;
      else
        dX = (nObs - 1) * dAvgSqrReturns 
           + ( N - (nObs - 1) ) * dInvAnnualReturnFreq  * dDefaultVar;

      dValue = varianceSwap.GetPayoffValue( dX / N ); 
    }

  
    m_pdRecoveryValues[nIdxT] = dValue * pdCF[nIdxT] * dTmp;
  }

}

size_t VarianceSwapMeshManager::GetObservationNumber(double dTime,
                bool &bIsTimeOnEvent)
{
  VarianceSwap &varianceSwap = m_varianceSwapParams.GetVarianceSwap();

  double dValuationTime  = m_varianceSwapParams.GetValuationTime();
  double dStartTime      = varianceSwap.GetStartTimeOfSamplingPeriod();
  double dMaturityTime   = varianceSwap.GetMaturityTime();  

  size_t nSamplingNumber     =  varianceSwap.GetNbSamplesUsed();
  size_t nNbSamplingReturns  =  varianceSwap.GetNbSamplingReturns();
  size_t nNbSamplesRemaining = nNbSamplingReturns;

  // Adjust if already past the sampling start date
  if ( numeric::IsAfter(dValuationTime, dStartTime) )
  {
    dStartTime = dValuationTime;
    nNbSamplesRemaining = nNbSamplingReturns 
                          - varianceSwap.GetNbSamplesUsed();
  }

  // Calculate what will be a "trading day"
  double dStep = (dMaturityTime - dStartTime) / nNbSamplesRemaining;

  // starting time
  double dT = dStartTime;  
  
  while ( numeric::IsBefore(dT, dTime) )
  {
    // Update for next observation
    nSamplingNumber++;
    dT += dStep;
  } // end while loop creating observation
  
  // current time is on an observation time
  if ( numeric::AreTimesEqual(dT, dTime) )
    bIsTimeOnEvent = true;
  else
    bIsTimeOnEvent = false;

  return nSamplingNumber;
}

void VarianceSwapMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();
}

size_t VarianceSwapMeshManager::GetNbRequestedTimes() const
{
  // Variance swaps need a minimum number of timesteps between
  // observations for short lived contracts in order to accurately
  // model the returns.  For long-dated contracts, the returns average
  // out, and accuracy is decent. Further, for long-dated contracts
  // (contracts with many samples), the time mesh will automatically
  // be generated with 2 (or more) timesteps per observation, so this
  // check only affects short-dated contracts.  The magic number 80 was 
  // determined by experiment, in particular, for vol swaps (actual return) 
  // with 5 or 6 samples (5 or 6 days to maturity).  One percent pricing
  // accuracy at the default time/space mesh settings was desired.
  size_t nNbTimes = m_params.GetNumParams()->GetNbTimeSteps();
  
  if (nNbTimes < 80)
    nNbTimes = 80;

  return nNbTimes;
}
