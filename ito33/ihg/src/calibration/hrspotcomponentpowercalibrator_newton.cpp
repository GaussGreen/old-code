/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrspotcomponentpowercalibrator.cpp
// Purpose:     calibrate with a spot component power Hazard rate
// Author:      Ito33
// Created:     2004/22/11
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_newton.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrspotcomponentpowercalibrator.cpp
   @brief calibration of time component of HRSpotComponentPower

 */
#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/useexception.h"
#include "ito33/error.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/newton2d.h"
#include "ito33/numeric/numericerror.h"


#include "ito33/finance/error.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/hrtimecomponentcalibrator.h"
#include "ihg/hrspotcomponentpowercalibrator_newton.h"

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;
extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

  using numeric::Exception;

void HazardRateSpotComponentPowerCalibratorNewton::operator () 
   (double dAlpha, double dBeta, double &dPrice1, double &dPrice2,
    double &dF)
 {

   shared_ptr<HazardRate> pHazardRate( new HazardRatePower(dAlpha,dBeta,m_dS0) ) ;

   m_theoreticalModel.SetHazardRate(pHazardRate);
   m_theoreticalModel.SetVolatility(m_pVol);

   dPrice1 = m_theoreticalModel.Compute(*m_pDeriv1)->GetPrice();
   dPrice2 = m_theoreticalModel.Compute(*m_pDeriv2)->GetPrice();

   dPrice1 -=  m_dMarketPrice1;
   dPrice2 -=  m_dMarketPrice2;

   dF = .5*( dPrice1*dPrice1*m_dScale1*m_dScale1 
           + dPrice2*dPrice2*m_dScale2*m_dScale2 );

  // Save best values in case calibration fails but the user still wants 
  // the best guess
  if (dF < m_dObjectif)
  {
    m_dObjectif = dF;
    m_dAlpha = dAlpha;
    m_dBeta = dBeta;
  }

 } //HazardRateSpotComponentPowerCalibratorNewton::operator ()

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorNewton::Calibrate
    (
      const finance::TermStructureDerivative& tsDerivative,
      const shared_ptr<Volatility>& pVolatility,
      size_t nIterMax,
      bool bUsePreviousSolution
    )
{ 
  ASSERT_MSG(!tsDerivative.GetAll().empty(),
    "Can't calibrate with an empty derivative term structure.");

  // Init parameters used when the calibration fails
  m_bConverged = false;
  m_tsDerivative  = tsDerivative;
  m_dObjectif = 1.e99;

  // Init everything else
  m_pVol          = pVolatility;
  m_pDeriv1       = tsDerivative.GetAll().front().get();
  m_pDeriv2       = tsDerivative.GetAll().back().get();
  m_dS0           = m_pDeriv1->GetSessionData()->GetSpotSharePrice();
  m_dMarketPrice1 = m_pDeriv1->GetMarketPrice();;
  m_dMarketPrice2 = m_pDeriv2->GetMarketPrice();

  m_dScale1 = fabs(m_dMarketPrice1);
  if (m_dScale1 > 1.e-6)
    m_dScale1 = 1.0/m_dScale1;
  else
    m_dScale1 = 1.0;

  m_dScale2 = fabs(m_dMarketPrice2);
  if (m_dScale2 > 1.e-6)
    m_dScale2 = 1.0/m_dScale2;
  else
    m_dScale2 = 1.0;


  double dAlphaLowerBound = 0.0;
  double dAlphaUpperBound = 10;
  double dBetaLowerBound  = 0.0;
  double dBetaUpperBound  = 2.0;

  // Check if the caller wants to use the previous solution as initial
  // guess.  Typically, this should be set to false the first time
  // calibrate is called.  If it is not, the default values set in
  // the constructor will still be used. Do not pass in m_dAlpha
  // and m_dBeta in case the calibration fails.
  double dAlpha;
  double dBeta;
  if ( !bUsePreviousSolution )
  {
    dAlpha = .1;
    dBeta  = .2;
  } 
  else
  {
    dAlpha = m_dAlpha;
    dBeta = m_dBeta;
  }

  //newton solver
  numeric::Newton2D solverN(dAlphaLowerBound, dAlphaUpperBound,
                           dBetaLowerBound, dBetaUpperBound,
                           1.e-8, nIterMax);

  numeric::NumericError err = solverN(*this, dAlpha, dBeta);


  // Don't know what to do when there is too many iterations, just report 
  // failure here
  if (   err == numeric::ITO33_NOT_CONVERGED 
      || err == numeric::ITO33_TOO_MANY_ITERATION )
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  // If the calibration failed, m_dAlpha and m_dBeta have already been set
  // in the () operator. Otherwise, we need to set them here.
  m_dAlpha = dAlpha;
  m_dBeta = dBeta;

  shared_ptr<HRSpotComponentPower>
    pSC(new HRSpotComponentPower(m_dBeta, m_dS0));

  shared_ptr<HazardRateCombo> 
    pHRCombo(new HazardRateCombo(pSC));

  HazardRateTimeComponentCalibrator calibratorHRTC;

  // The calibrator should update the hazard rate passed in (if any)
  calibratorHRTC.Calibrate(tsDerivative, pVolatility, pHRCombo);

  m_pHazardRate = pHRCombo;
  m_bConverged = true;

  return pHRCombo;  
}

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorNewton::GetHazardRate()
{
  // if the calibration worked, return the result
  if (m_bConverged)
    return m_pHazardRate;

  // if the calibration failed, make the best guess
  shared_ptr<HRSpotComponentPower>
    pSC(new HRSpotComponentPower(m_dBeta, m_dS0));

  shared_ptr<HazardRateCombo> 
    pHRCombo(new HazardRateCombo(pSC));

  HazardRateTimeComponentCalibrator calibratorHRTC;

  try
  {
    // The calibrator should update the hazard rate passed in (if any)
    calibratorHRTC.Calibrate(m_tsDerivative, m_pVol, pHRCombo);
  }
  catch(const ito33::numeric::Exception& /* e */)
  {
    // do nothing.  Just return whatever was calibrated so far.
  }

  // converged or not, if this function is called twice in a row, we can't
  // do any better
  m_bConverged = true;

  m_pHazardRate = pHRCombo;

  return pHRCombo;
    
}

} // namespace ihg

} // namespace ito33
