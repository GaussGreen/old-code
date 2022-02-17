/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrspotcomponentpowercalibrator_nag.cpp
// Purpose:     calibrate with a spot component power Hazard rate using NAG
// Created:     2004/22/11
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_beta.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrspotcomponentpowercalibrator_nag.cpp
   @brief calibration of time component of HRSpotComponentPower
 */
#include <cmath>

#include "ito33/array.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/hrspotcomponentpowercalibrator_beta.h"
#include "ihg/hrtimecomponentcalibrator.h"

#include "ito33/numeric/nonlinearsolver.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

  using numeric::Exception;

/// Calibrate an alpha with a given beta
class AlphaCalibrator
{
public:

  AlphaCalibrator(const finance::Derivative& deriv,
                  shared_ptr<Volatility> pVolatility, double dBeta)
    : m_deriv(deriv), m_dBeta(dBeta) 
  {
    m_model.SetVolatility(pVolatility);
    m_model.SetExternalFlagsToDefaults();

    m_dS0 = m_deriv.GetSessionData()->GetSpotSharePrice();

    m_dMarketPrice = m_deriv.GetMarketPrice(); 
    m_dScale = fabs( m_dMarketPrice );
    if (m_dScale > 1.e-6)
      m_dScale = 1. / m_dScale;
    else
      m_dScale = 1.;
  }

  double operator()(double dAlpha)
  {
    m_model.SetHazardRate(shared_ptr<HazardRate>(
                new HazardRatePower(dAlpha, m_dBeta, m_dS0)));

    return (m_model.Compute(m_deriv)->GetPrice() - m_dMarketPrice) * m_dScale;
  }     

  double Calibrate() 
  {
    numeric::RegulaFalsi calibrator(1.e-4, 100);
    
    return calibrator(*this, 0, 10.0);
  }


private:

  const finance::Derivative& m_deriv;
  TheoreticalModel m_model;
  double m_dBeta;
  double m_dS0;

  double m_dScale;
  double m_dMarketPrice;

  NO_COPY_CLASS(AlphaCalibrator);
};

double HazardRateSpotComponentPowerCalibratorBeta::operator()(double dBeta)
{
  AlphaCalibrator alphaCalibrator(*m_pDeriv1, m_pVolatility, dBeta);

  double dAlpha = alphaCalibrator.Calibrate();

  m_model.SetHazardRate(shared_ptr<HazardRate>(
            new HazardRatePower(dAlpha, dBeta, m_dS0)));

  double dError = m_model.Compute(*m_pDeriv2)->GetPrice()
                - m_dDeriv2MarketPrice;

  double dRelError = m_dScale * dError;

  // Save best values in case calibration fails but the user still wants 
  // the best guess
  if (dRelError < m_dObjectif)
  {
    m_dObjectif = dRelError;
    m_dAlpha = dAlpha;
    m_dBeta = dBeta;
  }

  return dRelError;
}

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorBeta::Calibrate
    (
      const finance::TermStructureDerivative& tsDerivative,
      const shared_ptr<Volatility>& pVolatility,
      size_t nIterMax
    )
{
  ASSERT_MSG(!tsDerivative.GetAll().empty(),
    "Can't calibrate with an empty derivative term structure.");

  // Init params that are saved even if calibration fails
  m_dAlpha = 0.0;
  m_dBeta = 0.0;
  m_dObjectif = 1.e99;
  m_bConverged = false;
  m_tsDerivative  = tsDerivative;

  // Save data for the objective function  
  m_pDeriv1 = tsDerivative.GetAll().front().get();
  m_pDeriv2 = tsDerivative.GetAll().back().get();

  m_dDeriv1MarketPrice = m_pDeriv1->GetMarketPrice();
  m_dDeriv2MarketPrice = m_pDeriv2->GetMarketPrice();
  
  m_pVolatility = pVolatility;
  m_model.SetVolatility(m_pVolatility);

  m_dS0 = m_pDeriv1->GetSessionData()->GetSpotSharePrice();
  
  m_dScale = fabs( m_dDeriv2MarketPrice );
  if (m_dScale > 1.e-6)
    m_dScale = 1. / m_dScale;
  else
    m_dScale = 1.;

  // the tolerance on the resulted function. This can't be too high, otherwise
  // it doesn't make much sense since it should somehow be limited by
  // the tolerance of the cds pricing precision
  const double dTolerance = 5 * 1.e-4;

  numeric::RegulaFalsi calibrator(dTolerance, nIterMax);

  double dOldBeta = 0.;
  double dOldFuncValue = operator()(dOldBeta);
  
  double dStep = 0.1;
  
  double dBeta = dOldBeta + dStep;
  double dFuncValue = operator()(dBeta);

  while ( !numeric::IsEqual(dBeta, 2) )
  { 
    // Check if the guess value is good enough
    if (dOldFuncValue * dFuncValue < 0 || fabs(dOldFuncValue) < dTolerance
                                       || fabs(dFuncValue) < dTolerance)
      break;

    // save the old values before update
    dOldBeta = dBeta;
    dOldFuncValue = dFuncValue;

    // update 
    dBeta += dStep;   
    
    // due to numerical round-off error, reset dBeta to 2 
    if (dBeta > 2)
      dBeta = 2;
    dFuncValue = operator()(dBeta);
  }

  m_dBeta = dOldBeta;
  dBeta = calibrator(*this, dOldBeta, dBeta);

  shared_ptr<HRSpotComponentPower>
    pSC(new HRSpotComponentPower(dBeta, m_dS0));

  shared_ptr<HazardRateCombo> 
    pHRCombo(new HazardRateCombo(pSC));

  HazardRateTimeComponentCalibrator calibratorHRTC;

  // The calibrator should update the hazard rate passed in (if any)
  calibratorHRTC.Calibrate(tsDerivative, m_pVolatility, pHRCombo);

  m_pHazardRate = pHRCombo;
  m_bConverged = true;

  return pHRCombo;
}

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorBeta::GetHazardRate()
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
    calibratorHRTC.Calibrate(m_tsDerivative, m_pVolatility, pHRCombo);
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
