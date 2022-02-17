/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrspotcomponentpowercalibrator_nag.cpp
// Purpose:     calibrate with a spot component power Hazard rate using NAG
// Created:     2004/22/11
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_nag.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrspotcomponentpowercalibrator_nag.cpp
   @brief calibration of HR with SpotComponentPower by using nag minimization
 */
#include <cmath>

#include "ito33/array.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"
#include "ito33/numeric/calibrator_nag.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"

#include "ihg/hrspotcomponentpowercalibrator_nag.h"
#include "ihg/hrtimecomponentcalibrator.h"

extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

  using numeric::Exception;

void HazardRateSpotComponentPowerCalibratorNAG::ComputeObjectif
     (const double *pdX, double& dObjectif, double*, bool)
{
  // extract the data needed to construct the hr spot component
  double dAlpha = pdX[0];
  double dBeta = pdX[1];

  // For some unknown reason, NAG can set beta to be < 0, even though
  // 0 is the lower bound. Might be round-off error
  if (dBeta < 0.0)
    dBeta = 0.0;
  if (dBeta > 2.0)
    dBeta = 2.0;

  shared_ptr<HazardRate> 
    pHazardRate( new HazardRatePower(dAlpha, dBeta, m_dS0) );

  m_theoreticalModel.SetHazardRate(pHazardRate);
  m_theoreticalModel.SetVolatility(m_pVolatility);

  double dPrice1 = m_theoreticalModel.Compute(*m_pDeriv1)->GetPrice();
  double dPrice2 = m_theoreticalModel.Compute(*m_pDeriv2)->GetPrice();

  // Calculate relative error
  double dError1 = (m_dDeriv1MarketPrice - dPrice1) * m_dScale1;

  double dError2 = (m_dDeriv2MarketPrice - dPrice2) * m_dScale2;
  
  // return least squares error
  dObjectif = dError1 * dError1 + dError2 * dError2;

  // Save best values in case calibration fails but the user still wants 
  // the best guess
  if (dObjectif < m_dObjectif)
  {
    m_dObjectif = dObjectif;
    m_dAlpha = dAlpha;
    m_dBeta = dBeta;
  }

}

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorNAG::Calibrate
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

  // Save data for the objective function , use only the first and the last cds
  m_pDeriv1 = tsDerivative.GetAll().front().get();
  m_pDeriv2 = tsDerivative.GetAll().back().get();

  m_dDeriv1MarketPrice = m_pDeriv1->GetMarketPrice();
  m_dDeriv2MarketPrice = m_pDeriv2->GetMarketPrice();

  // The scales
  m_dScale1 = fabs( m_dDeriv1MarketPrice );
  if (m_dScale1 > 1.e-5)
    m_dScale1 = 1.0/m_dScale1;
  else
    m_dScale1 = 1.0;

  m_dScale2 = fabs( m_dDeriv2MarketPrice );
  if (m_dScale2 > 1.e-5)
    m_dScale2 = 1.0/m_dScale2;
  else
    m_dScale2 = 1.0/m_dScale2;

  m_dS0 = m_pDeriv1->GetSessionData()->GetSpotSharePrice();

  m_pVolatility = pVolatility;

  // Setup the NAG parameters
  size_t nNbX = 2;
  Array<double> pdG(nNbX);
  Array<double> pdX(nNbX);
  
  // Setup the bounds of the parameters
  Array<double> pdLowerBounds(nNbX);
  Array<double> pdUpperBounds(nNbX);
  pdLowerBounds[0] = 0.0;
  pdLowerBounds[1] = 0.0;
  pdUpperBounds[0] = 1;
  pdUpperBounds[1] = 1.99;
  
  // calibrate a flat hazard rate and use it as initial guess
  finance::TermStructureDerivative tsDerivativeTmp;

  tsDerivativeTmp.Add(tsDerivative.GetAll().front());

  HazardRateTimeComponentCalibrator calibratorHRTC;

  double dAlpha = calibratorHRTC.Calibrate(tsDerivativeTmp, pVolatility)
                ->GetTimeComponentValues()[0]; 
  
  pdX[0] = dAlpha;
  pdX[1] = 0;

  
  // Tolerance used to check if calibrated result is good enough, regardless
  // the nag error code. This can't be too small and should be limited by
  // the square of the accuracy of the cds pricing.
  const double TOLERANCEOBJECTIF = 1.e-8;

  // Construct the nag calibrator
  numeric::CalibratorNAG calibrator(nNbX, pdLowerBounds.Get(), 
                       pdUpperBounds.Get(), 1.e-6, nIterMax);
    
  // calibrate with the first initial guess
  calibrator.Calibrate(*this, pdX.Get());
  
  // if failed with the initial guess, try with another by setting beta to 
  // upper bound
  double dFinalObjectif; 

  // check if the end result is good enough
  dFinalObjectif = calibrator.GetObjectif();
  if (calibrator.GetErrorCode() && dFinalObjectif > TOLERANCEOBJECTIF) 
  {
    // reset the intial guess, and start again
    pdX[0] = dAlpha;
    pdX[1] = 1.99999;
    calibrator.Calibrate(*this, pdX.Get());
  
    // check if the end result is good enough
    dFinalObjectif = calibrator.GetObjectif();

    if (calibrator.GetErrorCode() && dFinalObjectif > TOLERANCEOBJECTIF)
      throw EXCEPTION_MSG(ITO33_MAX_ITER, 
                          TRANS("Nag calibration failed!")); 
  }

  // Get the calibrated result
  dAlpha = pdX[0];
  
  double dBeta = pdX[1];

  shared_ptr<HRSpotComponentPower>
    pSC(new HRSpotComponentPower(dBeta, m_dS0));

  shared_ptr<HazardRateCombo> 
    pHRCombo(new HazardRateCombo(pSC));

  // The calibrator should update the hazard rate passed in (if any)
  calibratorHRTC.Calibrate(tsDerivative, m_pVolatility, pHRCombo);

  m_pHazardRate = pHRCombo;
  m_bConverged = true;

  return pHRCombo;
}

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorNAG::GetHazardRate()
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
