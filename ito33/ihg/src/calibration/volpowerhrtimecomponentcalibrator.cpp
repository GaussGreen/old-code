/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volpowerhrtimecomponentcalibrator.cpp
// Purpose:     calibrator on a power vol and a hr with time component
// Author:      ITO 33
// Created:     2005/01/05
// RCS-ID:      $Id: volpowerhrtimecomponentcalibrator.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004-  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructure.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilityflat.h"

#include "ihg/hrtimecomponentcalibrator.h"
#include "ihg/volflathrtimecomponentcalibrator.h"
#include "ihg/volpowerhrtimecomponentcalibrator.h"
#include "ihg/volpowercalibrator.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{


using ito33::numeric::Exception;

void VolPowerHRTimeComponentCalibrator::Calibrate
     (
       const finance::Derivative& derivative1, 
       const finance::Derivative& derivative2, 
       const finance::TermStructureDerivative& tsDerivative,
       shared_ptr<HazardRateWithTimeComponent> pHazardRate
     )
{ 

  // Init output objects to default values in case calibration fails
  double dS0 = derivative1.GetSessionData()->GetSpotSharePrice();
  m_pVolatility = make_ptr( new VolatilityPower(0.0, 0.0, dS0) );

  if (pHazardRate)
    m_pHazardRate = pHazardRate;
  else
  {
    Date valuationDate = derivative1.GetSessionData()->GetValuationDate();
    double dValue = 0.0;
    m_pHazardRate = make_ptr( new HazardRateTimeOnly
                                  (&valuationDate, &dValue, 1) );
  }

  // Try to get a good first guess for the hazard rate.  Assume that the
  // first derivative can give a good approximation to the volatility
  // alpha.  If this fails, use the default settings from above.
  try 
  {
    VolFlatHRTimeComponentCalibrator calibratorHR;
    calibratorHR.Calibrate(derivative1, tsDerivative, pHazardRate);
    m_pHazardRate = calibratorHR.GetHazardRate();
    double dVolatility = calibratorHR.GetVolatility()->GetValue();
    m_pVolatility = make_ptr( new VolatilityPower(dVolatility, 0.0, dS0) );
  }
  catch(const ito33::numeric::Exception&)
  {
    // do nothing, which implies using the default hazard rate from above
  }

  // Cache the market prices for the IsConverged function in case
  // the market prices are computed (eg. option with implied vol set)
  // instead of simply stored.
  m_dDeriv1MarketPrice = derivative1.GetMarketPrice();
  m_dDeriv2MarketPrice = derivative2.GetMarketPrice();

  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivative.GetAll();
  m_pdTSMarketPrices.resize( pElements.size() );
  size_t nCounter = 0;
  finance::TermStructureDerivative::Elements::const_iterator ppDeriv;
  for ( ppDeriv = pElements.begin(); ppDeriv != pElements.end(); ++ppDeriv)
  {
    m_pdTSMarketPrices[nCounter] = (*ppDeriv)->GetMarketPrice();
    nCounter++;
  }

  // prepare for iterating
  const size_t nNbMaxIterations = 40;
  size_t nIdxIteration = 0;
  bool bRestarted = false;
  bool bConverged = false;

  HazardRateTimeComponentCalibrator calibratorHR;
  VolPowerCalibrator calibratorVol;

  while ( nIdxIteration < nNbMaxIterations && !bConverged )
  { 

    // Fit the volatility using the two derivatives
    try
    {
      m_pVolatility = calibratorVol.Calibrate(derivative1, 
                                              derivative2, 
                                              m_pHazardRate);
    }
    catch(const ito33::numeric::Exception& e)
    {
      // Can have problems if the vol is near zero, since the hazard rate 
      // may be set too high. Try again (once only) with guess of zero.
      if ( bRestarted )
        throw e;
      else
      {
        // If it failed on the first iteration, use the original guess 
        // for vol, and don't restart
        if (nIdxIteration > 1)
        {
          m_pVolatility = make_ptr( new VolatilityPower(0.0, 0.0, dS0) );
          bRestarted = true;
        }
      } // if restarted
    }

    // Fit the hazard rate using the derivative term structure
    try
    {      
      m_pHazardRate = calibratorHR.Calibrate(tsDerivative, m_pVolatility, pHazardRate);
    }
    catch(const ito33::numeric::Exception& e)
    {
      // Calibrating time only hazard rate should be easy.  If 
      // it fails, something is really wrong
      throw e;
    }

    // Check errors
    bConverged = IsConverged(derivative1, derivative2, tsDerivative);

    nIdxIteration++;
  }

  CHECK_COND( nIdxIteration < nNbMaxIterations, ITO33_MAX_ITER);
}


bool VolPowerHRTimeComponentCalibrator::IsConverged
                (const finance::Derivative& derivative1, 
                 const finance::Derivative& derivative2, 
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
  double dTol = 1.e-3;
  
  // Check the first derivative 
  double dDeriv1Price = model.Compute(derivative1)->GetPrice();
  double dScale = fabs( m_dDeriv1MarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  double dError = fabs( dDeriv1Price - m_dDeriv1MarketPrice );
  if ( dError/dScale > dTol)
    return false;

  // Check the second derivative
  double dDeriv2Price = model.Compute(derivative2)->GetPrice();
  dScale = fabs( m_dDeriv2MarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  dError = fabs( dDeriv2Price - m_dDeriv2MarketPrice );
  if ( dError/dScale > dTol)
    return false;

  // Check the term structure prices
  const finance::TermStructureDerivative::Elements& 
    pElements = tsDerivative.GetAll();
  finance::TermStructureDerivative::Elements::const_iterator ppDerivs;
  size_t nCounter = 0;
  for ( ppDerivs = pElements.begin(); ppDerivs != pElements.end(); ++ppDerivs)
  {
    const finance::Derivative& deriv = *(*ppDerivs);

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
