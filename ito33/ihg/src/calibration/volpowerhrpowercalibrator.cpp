/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volflathrspotcomponentpowercalibrator.cpp
// Purpose:     calibrator on a vol flat and a hr with spot component power
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: volpowerhrpowercalibrator.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/volatilitypower.h"

#include "ihg/impliedvol.h"
#include "ihg/hrspotcomponentpowercalibrator.h"
#include "ihg/volpowercalibrator.h"
#include "ihg/volflathrspotcomponentpowercalibrator.h"
#include "ihg/volflathrtimecomponentcalibrator.h"
#include "ihg/hrtimecomponentcalibrator.h"

#include "ihg/volpowerhrpowercalibrator.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

using ito33::numeric::Exception;


void VolPowerHRPowerCalibrator::Calibrate
                (const finance::Derivative& deriv1,
                 const finance::Derivative& deriv2,
                 const finance::TermStructureDerivative& tsDerivs)
{ 
  ASSERT_MSG(!tsDerivs.GetAll().empty(),
    "Can't calibrate with an empty derivative term structure.");

  // Over-all strategy
  // 1) Find good starting guess for hazard rate
  // 2) Fit the two derivatives and the term structure by iterating:
  //     - find new power vol to fit the two derivatives using the
  //       hazard rate
  //     - find new hazard rate (alpha, beta, time components) to fit the 
  //       term structure using using the current power vol guess
  //     - stop when the prices have converged
  //
  // See also the volflathrwithspotcomponentpower calibrator

  // Init output objects to default values in case calibration fails
  double dS0 = deriv1.GetSessionData()->GetSpotSharePrice();
  m_pVolatility = make_ptr( new VolatilityPower(0.0, 0.0, dS0) );
  m_pHRSpotComponentPower = make_ptr( new HRSpotComponentPower(0.0, dS0) );
  m_pHazardRate = make_ptr( new HazardRateCombo(m_pHRSpotComponentPower) );

  // Try to get a good first guess for the hazard rate. Assume that the 
  // first derivative gives a better estimate for a flat vol, which is
  // a reasonable proxy for alpha in the power volatility
  double dVolatility = 0.0;
  try
  {

    //VolFlatHRTimeComponentCalibrator calibratorVolFlatHRTime;
    //calibratorVolFlatHRTime.Calibrate(deriv1, tsDerivs, NULL);
    //dVolatility = calibratorVolFlatHRTime.GetVolatility()->GetValue();

    VolFlatHRSpotComponentPowerCalibrator calibratorVolFlatHRPower;
    calibratorVolFlatHRPower.Calibrate(deriv1, tsDerivs);
    m_pHazardRate = calibratorVolFlatHRPower.GetHazardRate();
    m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower> 
                              ( m_pHazardRate->GetSpotComponent() );
    dVolatility = calibratorVolFlatHRPower.GetVolatility()->GetValue();
    m_pVolatility = make_ptr( new VolatilityPower(dVolatility, 0.0, dS0) );
    
  }
  catch (const ito33::numeric::Exception& /* e */)
  {
    // Estimating vol should only fail for weird cases.  The code below is 
    // more than likely to also fail, but might as well try a boundary case.
    //dVolatility = 0.0;

    // Estimating the hr should only fail for weird cases.  The code below is 
    // more than likely to also fail, but might as well try a boundary case.
    // Have already set hazard rate and volatility to default zero values
  }

  // Cache the market prices for the IsConverged function in case
  // the market prices are computed (eg. option with implied vol set)
  // instead of simply stored.
  m_dDeriv1MarketPrice = deriv1.GetMarketPrice();
  m_dDeriv2MarketPrice = deriv2.GetMarketPrice();

  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivs.GetAll();
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
  bool bSpotComponentConverged = false;
  bool bConverged = false;

  HazardRateSpotComponentPowerCalibrator calibratorHR;
  VolPowerCalibrator calibratorVol;

  while ( nIdxIteration < nNbMaxIterations && !bConverged)
  { 

    nIdxIteration++;

    // Fit the volatility using the two derivatives
    try
    {
      if (nIdxIteration >= 8 && nIdxIteration % 4 == 0)
        m_pVolatility = calibratorVol.Calibrate(deriv2, 
            *tsDerivs.GetAll().front(), m_pHazardRate);
      else if (nIdxIteration >= 8 && nIdxIteration % 2 == 0)
        m_pVolatility = calibratorVol.Calibrate(deriv1, 
            *tsDerivs.GetAll().back(), m_pHazardRate);
      else
        m_pVolatility = calibratorVol.Calibrate(deriv1, deriv2, m_pHazardRate);
    }
    catch(const ito33::numeric::Exception& e)
    {
      // If the vol was already reset to zero, and we are still having problems,
      // just throw the exception
      if ( bRestarted )
        throw e;
      else
      {
        // If it failed, use the last best guess.  Only do this so many times
        // before resetting the vol and essentially starting over
        if (nIdxIteration > 16)
        {
          m_pVolatility = make_ptr( new VolatilityPower(0.0, 0.0, dS0) );
          bRestarted = true;
        }
        else
        {
          m_pVolatility = calibratorVol.GetVolatility();
        }
      } // if restarted
    }

    // Sometimes the derivative prices are more sensitive than the
    // term structure prices.  In this case, the prices are converged
    // at this point. It the code were to continue, a 'loose' hazard
    // rate is found that fits the term structure, but convergence
    // to the other two derivatives is lost.  Give the code
    // some time to compute a hazard rate though
    if (nIdxIteration > 8)
    {
      bConverged = bSpotComponentConverged;
      if (bSpotComponentConverged)
        bConverged = IsConverged(deriv1, deriv2, tsDerivs);
  
      if (bConverged)
        continue;
    }

    // Fit the hazard rate using the term structure
    bSpotComponentConverged = false;
    try
    {      
      m_pHazardRate = calibratorHR.Calibrate(tsDerivs, m_pVolatility);
      m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower>
                                ( m_pHazardRate->GetSpotComponent() );   

      bSpotComponentConverged = true;
    }
    catch(const ito33::numeric::Exception& /* e */)
    {
      // Even if it failed, it may have been close.  Keep trying these
      // outer iterations with the last best guess
      m_pHazardRate = calibratorHR.GetHazardRate();
      m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower> 
                                ( m_pHazardRate->GetSpotComponent() );
    }

    // Check for convergence
    bConverged = bSpotComponentConverged;
    if (bSpotComponentConverged)
      bConverged = IsConverged(deriv1, deriv2, tsDerivs);

    // If it hasn't converged by 30 iterations, try resetting the hazard rate
    // to see what happens
    if (nIdxIteration == 30 && bConverged == false)
    {
      m_pHRSpotComponentPower = make_ptr( new HRSpotComponentPower(0.0, dS0) );
      m_pHazardRate = make_ptr( new HazardRateCombo(m_pHRSpotComponentPower) );
    }      

  } // while not converged

  if ( !bConverged ) 
    throw EXCEPTION_MSG
          (
            ITO33_MAX_ITER,
            TRANS("Cannot calibrate a power volatility and spot component"
                  " power hazard rate to the specified derivatives and"
                  " termstructure.")
          );

} //VolPowerHRPowerCalibrator::Calibrate


bool VolPowerHRPowerCalibrator::IsConverged
                (const finance::Derivative& deriv1, 
                 const finance::Derivative& deriv2, 
                 const finance::TermStructureDerivative& tsDerivs)
{
  ASSERT_MSG(!tsDerivs.GetAll().empty(),
    "Can't calibrate with an empty derivative term structure.");

  // Setup a model for pricing
  TheoreticalModel model;
  model.SetVolatility(m_pVolatility);
  model.SetHazardRate(m_pHazardRate);
  model.SetExternalFlagsToDefaults();

  // Make sure the relative error of each contract is less than a tolerance.
  // Since we are using normal pricing params, the tolerance must be slightly
  // larger than the usual pricing accuracy.

  double dTol = 1.e-3;

  // Check the first derivative
  double dDeriv1Price = model.Compute(deriv1)->GetPrice();
  double dScale = fabs( m_dDeriv1MarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  double dError = fabs( dDeriv1Price - m_dDeriv1MarketPrice );
  if ( dError/dScale > dTol)
    return false;

  // Check the second derivative
  double dDeriv2Price = model.Compute(deriv2)->GetPrice();
  dScale = fabs( m_dDeriv2MarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  dError = fabs( dDeriv2Price - m_dDeriv2MarketPrice );
  if ( dError/dScale > dTol)
    return false;

  // Check the term structure prices
  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivs.GetAll();
  finance::TermStructureDerivative::Elements::const_iterator ppDeriv;
  size_t nCounter = 0;
  for ( ppDeriv = pElements.begin(); ppDeriv != pElements.end(); ++ppDeriv)
  {
    const finance::Derivative& deriv = *(*ppDeriv);

    double dPrice = model.Compute(deriv)->GetPrice();
    
    dScale = fabs( m_pdTSMarketPrices[nCounter] );
    if (dScale < 1.e-6)
      dScale = 1.0;

    dError = fabs( dPrice - m_pdTSMarketPrices[nCounter] );
    if ( dError/dScale > dTol)
      return false;

    nCounter++;
  }

  // If we make it here, all contracts have converged
  return true;
  
}


} // namespace ihg

} // namespace ito33
