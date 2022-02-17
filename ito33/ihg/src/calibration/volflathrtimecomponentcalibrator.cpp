/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volflathrtimecomponentcalibrator.cpp
// Purpose:     calibrator on a vol flat and a hr with time component
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: volflathrtimecomponentcalibrator.cpp,v 1.15 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructure.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ihg/impliedvol.h"
#include "ihg/hrtimecomponentcalibrator.h"
#include "ihg/volflathrtimecomponentcalibrator.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{


using ito33::numeric::Exception;

void VolFlatHRTimeComponentCalibrator::Calibrate
     (
       const finance::Derivative& derivative, 
       const finance::TermStructureDerivative& tsDerivative,
       shared_ptr<HazardRateWithTimeComponent> pHazardRate
     )
{ 

  // Init output objects to default values in case calibration fails
  m_pVolatility = make_ptr( new VolatilityFlat(0.0) );
  if (pHazardRate)
    m_pHazardRate = pHazardRate;
  else
  {
    Date valuationDate = derivative.GetSessionData()->GetValuationDate();
    double dValue = 0.0;
    m_pHazardRate = make_ptr( new HazardRateTimeOnly
                                  (&valuationDate, &dValue, 1) );
  }

  // Creates a calibrator
  HazardRateTimeComponentCalibrator hrCalibrator;

  // If there is no spot component (passed in hazard rate), the cds prices
  // are independent of the volatility.  Hence, just need to fit the cds
  // prices using an arbitrary vol value, then fit the option using the
  // new hazard rate
  if ( !pHazardRate ) // hazard rate is time only
  {
    // Calibrate the hazard rate

    // Creates a temporary Vol
    shared_ptr<Volatility> pVolatility( new VolatilityFlat(0.2) );

    m_pHazardRate = hrCalibrator.Calibrate(tsDerivative, pVolatility);
  
    // Calibrate the vol flat
    ImpliedVol impliedVol(derivative, m_pHazardRate);

    m_pVolatility = make_ptr( new VolatilityFlat( impliedVol.Compute() ) );

    return;
  }
  
  // Iterate over flat vol

  const size_t nNbMaxIterations = 40;

  /*
    We are using the implied vol computed with a zero hazard rate. This gives 
    the upper bound of the to be calibrated vol. If the upper bound is limited,
    this will gurantee the convergence of the following algorithm (suppose that 
    there is at most one cds which has a smaller maturity than the option). 
    
    Otherwise, we can begin with the lower bound of the vol. For the moment, 
    we use just zero, but the lower bound can be computed by using impliedvol 
    with the hazard rate implied by a zero volatility. This should also work.
  */
  double dVolatility;

  try 
  { 
    // upper bound
    shared_ptr<HazardRate> pHRFlatTmp(new HazardRateFlat(0.));
    ImpliedVol impliedVol(derivative, pHRFlatTmp);
    
    dVolatility = impliedVol.Compute();
  }
  catch(const ito33::numeric::Exception&)
  {
    dVolatility = 0;
  }

  // Cache the market prices for the IsConverged function in case
  // the market prices are computed (eg. option with implied vol set)
  // instead of simply stored.
  m_dDerivMarketPrice = derivative.GetMarketPrice();

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

  // initialize
  m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );
 
  size_t nIdxIteration = 0;
  bool bConverged = false;

  while ( nIdxIteration < nNbMaxIterations && !bConverged )
  { 

    // Fit first the term structure
    m_pHazardRate 
      = hrCalibrator.Calibrate(tsDerivative, m_pVolatility, pHazardRate);
    
    // Fit then the derivative    
    ImpliedVol impliedVol(derivative, m_pHazardRate);
    dVolatility = impliedVol.Compute();
    m_pVolatility = make_ptr( new VolatilityFlat(dVolatility) );
      
    // Check errors
    bConverged = IsConverged(derivative, tsDerivative);

    nIdxIteration++;
  }

  if ( nIdxIteration == nNbMaxIterations) 
    throw EXCEPTION_MSG
          (
            ITO33_MAX_ITER,
            TRANS("Can't calibrate the actual instruments with current model.")
          );

}


bool VolFlatHRTimeComponentCalibrator::IsConverged
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
  double dTol = 1.e-3;
  
  // Check the derivative first
  double dDerivPrice = model.Compute(derivative)->GetPrice();
  double dScale = fabs( m_dDerivMarketPrice );
  if (dScale < 1.e-6)
      dScale = 1.0;
 
  double dError = fabs( dDerivPrice - m_dDerivMarketPrice );
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
