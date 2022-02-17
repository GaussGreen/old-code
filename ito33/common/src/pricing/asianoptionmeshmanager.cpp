/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/asianoptionmeshmanager.cpp
// Purpose:     Asian option mesh manager for backward PDE problems
// Author:      Ito 33 Canada
// Created:     2005/06/02
// RCS-ID:      $Id: asianoptionmeshmanager.cpp,v 1.5 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"
#include "ito33/constants.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/asianoptionparams.h"
#include "ito33/pricing/asianoptionmeshmanager.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

AsianOptionMeshManager::AsianOptionMeshManager
  (AsianOptionParams& params, Model& model)
 : OptionMeshManager(params, model), m_asianOptionParams(params)
{ 
}

void AsianOptionMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  AsianOption &asianOpt = m_asianOptionParams.GetAsianOption();

  // Average of the path not current average
  double dAverage   = m_asianOptionParams.GetAverage(); 
  double dStartTime = asianOpt.GetAverageStartTime();

  // Compute discount factor
  Array<double> pdCF(m_nNbTimes);
      
  m_asianOptionParams.GetYieldCurve()->GetCompoundFactor(m_pdTimes.Get(), 
                                                         pdCF.Get(), 
                                                         m_nNbTimes);
  double dTmp = 1. / pdCF[m_nNbTimes - 1];
  size_t nObs    = 0;
  double dAdjust = 0.0;

  for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
  {
    // We are before the start of the observation including it
    // Average is simply 0 then
    if ( numeric::IsEqualOrBefore(m_pdTimes[nIdxT], dStartTime ) )
    {
      dAverage = 0.0;
    }
    // Anytime after T_ini
    else
    {
      nObs    = GetObservationNumber( m_pdTimes[nIdxT] );
      dAdjust = double(nObs) / double( asianOpt.GetNbSamplingAverages() );
    }
    
    switch ( asianOpt.GetOptionType() )     
    {
    case pricing::AsianOption_FloatingStrikeCall: //max( S - A, 0.0 )
      m_pdRecoveryValues[nIdxT] = 0.;
      break;
        
    case pricing::AsianOption_FloatingStrikePut: //max( A - S, 0.0)
      {
        if ( asianOpt.GetExerciseType() ==  ExerciseType_American )
          m_pdRecoveryValues[nIdxT] = dAverage;
        else
          m_pdRecoveryValues[nIdxT] = dAverage * dAdjust * pdCF[nIdxT] * dTmp; 

        break;
      }
    case pricing::AsianOption_FixedStrikeCall: //max(A - K, 0.0)
      {
        double dStrike = asianOpt.GetStrike();

        if ( asianOpt.GetExerciseType() ==  ExerciseType_American )
          m_pdRecoveryValues[nIdxT] = std::max(dAverage - dStrike, 0.0);
        else 
          m_pdRecoveryValues[nIdxT] = std::max(dAverage * dAdjust - dStrike, 0.0) 
                                      * pdCF[nIdxT] * dTmp;

        break;
      }
    case pricing::AsianOption_FixedStrikePut://max( K - A, 0.0)
      {
        double dStrike = asianOpt.GetStrike();

        if ( asianOpt.GetExerciseType() ==  ExerciseType_American )
          m_pdRecoveryValues[nIdxT] = std::max(dStrike - dAverage, 0.0); 
        else 
          m_pdRecoveryValues[nIdxT] = 
            std::max( dStrike - dAverage * dAdjust, 0.0 ) * pdCF[nIdxT] * dTmp;

        break;
      }
    default:  
      FAIL("Asian Option pricer doesn't handle current option type.");    
    }
  }

}

size_t AsianOptionMeshManager::GetObservationNumber(double dTime)
{
  AsianOption &asianOpt = m_asianOptionParams.GetAsianOption();

  double dValuationTime  = m_asianOptionParams.GetValuationTime();
  double dStartTime      = asianOpt.GetAverageStartTime();
  double dEndTime        = asianOpt.GetAverageEndTime();  

  size_t nSamplingNumber     =  asianOpt.GetNbSamplesUsed();
  size_t nNbSamplingAverages =  asianOpt.GetNbSamplingAverages();
  size_t nNbSamplesRemaining = nNbSamplingAverages;
  size_t nNbSamplesUsed      = asianOpt.GetNbSamplesUsed();

  // When passed the end time simply return
  // the total number of observation
  if ( numeric::IsAfter(dTime, dEndTime) )
    return nNbSamplingAverages;

  // Adjust if already past the sampling start date
  if ( numeric::IsAfter(dValuationTime, dStartTime) )
  {
    dStartTime = dValuationTime;
    nNbSamplesRemaining = nNbSamplingAverages - nNbSamplesUsed;
  }

  // Calculate what will be a "trading day"
  double dStep = (dEndTime - dStartTime) / (nNbSamplesRemaining );

  // starting time
  double dT = dStartTime;  
  
  while ( numeric::IsEqualOrBefore(dT, dTime) )
  {
    // Update for next observation
    nSamplingNumber++;
    dT += dStep;
  } // end while loop creating observation

  return nSamplingNumber - 1;
}
