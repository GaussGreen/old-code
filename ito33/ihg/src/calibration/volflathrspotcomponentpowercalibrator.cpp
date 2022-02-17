/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volflathrspotcomponentpowercalibrator.cpp
// Purpose:     calibrator on a vol flat and a hr with spot component power
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: volflathrspotcomponentpowercalibrator.cpp,v 1.9 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"

#include "ito33/finance/option.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratecombo.h"

#include "ihg/impliedvol.h"
#include "ihg/hrspotcomponentpowercalibrator.h"
#include "ihg/volflathrspotcomponentpowercalibrator.h"
#include "ihg/volflathrtimecomponentcalibrator.h"


extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

using ito33::numeric::Exception;


void VolFlatHRSpotComponentPowerCalibrator::Calibrate
                (const finance::Derivative& derivative, 
                 const finance::TermStructureDerivative& tsDerivative)
{ 
  // Over-all strategy
  // 1) Find good starting guess for vol
  // 2) Fit option and cds term structure by iterating:
  //     - find hazard rate (alpha, beta, time components) to fit the cds 
  //       term structure using using the fixed vol
  //     - find new fixed vol to fit the option using calibrated hazard rate
  //     - stop when the prices have converged
  //
  // It is also possible to stop iterating when the vol stops changing.
  // However, convergence of prices may occur before a vol tolerance is
  // reached, and the vols may never converge due to pricing and
  // calibration errors during the hazard rate fitting.

  // Init output objects to default values in case calibration fails
  m_pVolatility = make_ptr( new VolatilityFlat(0.0) );
  double dS0 = derivative.GetSessionData()->GetSpotSharePrice();
  m_pHRSpotComponentPower = make_ptr( new HRSpotComponentPower(0.0, dS0) );
  m_pHazardRate = make_ptr( new HazardRateCombo(m_pHRSpotComponentPower) );

  // Try to get a good first guess for the vol.  At least two approaches:
  // 1) Fit vol and hazard rate with time component. Should give a very
  //    good approx. to vol, but can be slow, and may not work
  // 2) Find implied vol for random guess at hazard rate.  Should be
  //    very quick, but may not be as accurate
  // Try the first approach
  double dVolatility = 0.0;
  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivative.GetAll();
  finance::TermStructureDerivative tsDerivTmp;
  tsDerivTmp.Add(pElements.front());
  tsDerivTmp.Add(pElements.back());

  try
  {
    VolFlatHRTimeComponentCalibrator calibratorVolFlatHRTime;
    calibratorVolFlatHRTime.Calibrate(derivative, tsDerivTmp,
                                    shared_ptr<HazardRateWithTimeComponent>());
    dVolatility = calibratorVolFlatHRTime.GetVolatility()->GetValue();
    
    // The vol is often estimated too high for small vols, and too low
    // for large vols. Adjust down if the hr shows structure. Do not
    // adjust up, since we would rather be too low than too high
    shared_ptr<HazardRateWithTimeComponent> pHRTmp
      = calibratorVolFlatHRTime.GetHazardRate();
    std::vector<double> pdVals = pHRTmp->GetTimeComponentValues();
    double dTmp1 = pdVals[0];
    double dTmp2 = pdVals[1];
    double dRatio = dTmp2/dTmp1;
    if (dRatio < 1.0 && dVolatility < 0.61)
      dVolatility *= dRatio;

  }
  catch (const ito33::numeric::Exception& /* e */)
  {
    // Estimating vol should only fail for weird cases.  The code below is 
    // more than likely to also fail, but might as well try a boundary case.
    dVolatility = 0.0;
  }

  // Cache the market prices for the IsConverged function in case
  // the market prices are computed (eg. option with implied vol set)
  // instead of simply stored.
  m_dDerivMarketPrice = derivative.GetMarketPrice();
  m_pdTSMarketPrices.resize( pElements.size() );
  size_t nCounter = 0;
  finance::TermStructureDerivative::Elements::const_iterator ppDeriv;
  for ( ppDeriv = pElements.begin(); ppDeriv != pElements.end(); ++ppDeriv)
  {
    m_pdTSMarketPrices[nCounter] = (*ppDeriv)->GetMarketPrice();
    nCounter++;
  }


  // Iteration section. During the iterations, save the current 
  // guesses of the vol and hazard rate in case something fails.  
  // The user may want to retrieve these "last best guesses".
  const size_t nNbMaxIterations = 40;
  size_t nIdxIteration = 0;
  bool bRestarted = false;
  bool bConverged = false;
  std::vector<double> pdVolatilities(nNbMaxIterations+1);
  pdVolatilities[0] = dVolatility;

  HazardRateSpotComponentPowerCalibrator calibratorSP;

  m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );

  while ( nIdxIteration < nNbMaxIterations && !bConverged)
  { 

    nIdxIteration++;

    // Fit the spot component (alpha and beta) to the term structure
    try
    {      
      m_pHazardRate = calibratorSP.Calibrate(tsDerivative, m_pVolatility);
      m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower> 
                                ( m_pHazardRate->GetSpotComponent() );
    }
    catch(const ito33::numeric::Exception& /* e */)
    {
      // Even if it failed, it may have been close.  Keep trying these
      // outer iterations with the last best guess
      m_pHazardRate = calibratorSP.GetHazardRate();
      m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower>
                                ( m_pHazardRate->GetSpotComponent() );

      // If we are still having problems after 10 iterations, something
      // is probably wrong
      if (nIdxIteration > 10 && bRestarted == false)
      {
        dVolatility = 0.0;
        m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );
        bRestarted = true;
        continue;
      }
    }

    // Check for convergence
    if (nIdxIteration > 1)
    {
      bConverged = IsConverged(derivative, tsDerivative);
      if (bConverged)
        continue;
    }


    // Fit the single derivative
    try
    {      
      ImpliedVol impliedVol(derivative, m_pHazardRate);
      dVolatility = impliedVol.Compute();
    }
    catch(const ito33::numeric::Exception& e)
    {
      // The first iteration can have problems due to a bad first guess for 
      // the volatility. Can also have problems if the vol is near zero, 
      // since the hazard rate may be set too high while iterating. Try again
      // with guess of zero
      if ( bRestarted )
          throw e;
      else
      {
        dVolatility = 0.0;
        bRestarted = true;
      }
    }

    // Use the new volatility
    m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );

    // Check for convergence
    bConverged = IsConverged(derivative, tsDerivative);
    
    // Sometimes the solution oscillates around the exact answer. This
    // might help, and costs very little
    pdVolatilities[nIdxIteration] = dVolatility;
    if (nIdxIteration == 20 && bConverged == false)
    {
      double dVolatilityTmp = 1.0/6.0 *
          (pdVolatilities[nIdxIteration] + pdVolatilities[nIdxIteration-1]
         + pdVolatilities[nIdxIteration-2] + pdVolatilities[nIdxIteration-3]
         + pdVolatilities[nIdxIteration-4] + pdVolatilities[nIdxIteration-5]);
      if ( fabs(dVolatilityTmp - dVolatility) < 0.05)
        m_pVolatility = make_ptr( new VolatilityFlat(dVolatilityTmp) );
    }

    // Sometimes the convergence is just really slow from above. Restart
    // from below
    if (nIdxIteration == 30 && bConverged == false && bRestarted == false)
    {
      dVolatility = 0.0;
      m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );
      bRestarted = true;
    }


  } // while not converged

  if ( !bConverged ) 
    throw EXCEPTION_MSG
          (
            ITO33_MAX_ITER,
            TRANS("Can't calibrate the actual instruments with current model.")
          );

} //VolFlatHRSpotComponentPowerCalibrator::Calibrate


bool VolFlatHRSpotComponentPowerCalibrator::IsConverged
                (const finance::Derivative& derivative, 
                 const finance::TermStructureDerivative& tsDerivative)
{
  // Setup a model for pricing
  TheoreticalModel model;
  model.SetVolatility(m_pVolatility);
  model.SetHazardRate(m_pHazardRate);
  model.SetExternalFlagsToDefaults();

  // Make sure the relative error of each contract is less than a tolerance.
  // Since we are using normal pricing params, the tolerance must be slightly
  // larger than the usual pricing accuracy.

  // Check the single derivative first
  double dTol = 1.e-3;
  double dPrice = model.Compute(derivative)->GetPrice();
  double dScale = fabs( m_dDerivMarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  double dError = fabs( dPrice - m_dDerivMarketPrice );
  if ( dError/dScale > dTol)
    return false;

  // Check the term structure prices
  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivative.GetAll();
  finance::TermStructureDerivative::Elements::const_iterator ppDeriv;
  size_t nCounter = 0;
  for ( ppDeriv = pElements.begin(); ppDeriv != pElements.end(); ++ppDeriv)
  {
    const finance::Derivative& deriv = *(*ppDeriv);

    double dDerivPrice = model.Compute(deriv)->GetPrice();
    
    dScale = fabs( m_pdTSMarketPrices[nCounter] );
    if (dScale < 1.e-6)
      dScale = 1.0;

    dError = fabs( dDerivPrice - m_pdTSMarketPrices[nCounter] );
    if ( dError/dScale > dTol)
      return false;

    nCounter++;
  }

  // If we make it here, all contracts have converged
  return true;
  
}


} // namespace ihg

} // namespace ito33
