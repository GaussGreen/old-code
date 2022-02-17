/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrspotcomponentpowercalibrator_mixed.cpp
// Purpose:     calibrate with a spot component power hazard rate using
//              Newton. If it fails, use beta iteration
// Created:     2004/12/21
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_mixed.cpp,v 1.5 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrspotcomponentpowercalibrator_mixed.cpp
   @brief calibration of time component of HRSpotComponentPower
 */
#include "ito33/numeric/exception.h"

#include "ihg/hrspotcomponentpowercalibrator_mixed.h"
#include "ihg/hrspotcomponentpowercalibrator_beta.h"
#include "ihg/hrspotcomponentpowercalibrator_newton.h"


namespace ito33
{

namespace ihg
{

  using numeric::Exception;

shared_ptr<HazardRateCombo> 
HazardRateSpotComponentPowerCalibratorMixed::Calibrate
    (
      const finance::TermStructureDerivative& tsDerivative,
      const shared_ptr<Volatility>& pVolatility,
      size_t nIterMax
    )
{

  // Try calibrating with Newton (the quickest, but arguably least stable
  // algorithm currently implemented).  If it fails, fall back to the
  // beta iteration algorithm.

  try
  {
    m_pHazardRate = calibratorNewton.Calibrate(tsDerivative, pVolatility, 15);
    return m_pHazardRate;
  }
  catch (const ito33::numeric::Exception& /* e */)
  {
    // try next method. 
  }


  HazardRateSpotComponentPowerCalibratorBeta calibratorBeta;
  try
  {
    m_pHazardRate = calibratorBeta.Calibrate(tsDerivative, pVolatility, nIterMax);
  }
  catch (const ito33::numeric::Exception& e)
  {
    // Both methods failed. Throw the exception this time. Save the hazard rate
    // from the Newton method though
    m_pHazardRate = calibratorNewton.GetHazardRate();
    throw e;
  }

  return m_pHazardRate;
}


} // namespace ihg

} // namespace ito33
