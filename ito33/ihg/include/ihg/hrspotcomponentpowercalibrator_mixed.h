/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hrspotcomponentpowercalibrator_mixed.h
// Purpose:     calibrate spot component power of HRSpotComponentPower
// Author:      Ito33
// Created:     2004/12/21
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_mixed.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/hrspotcomponentpowercalibrator_mixed.h
   @brief calibrate spot component power of HRSpotComponentPower
 */

#ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_MIXED_H_
#define _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_MIXED_H_

#include "ito33/sharedptr.h"
#include "ito33/finance/termstructurederivative.h"
#include "ihg/hrspotcomponentpowercalibrator_newton.h"

namespace ito33
{

namespace ihg
{

  class Volatility;
  class HazardRateCombo;

/// class for spot power hazard rate calibration using a term structure
class HazardRateSpotComponentPowerCalibratorMixed
{

public:
  
  /**
     Calibrates a term structure with known volatility

     @param tsDeriv the derivative term structure to calibrate
     @param pVolatility the known volatility
     @param nIterMax maximum number of iterations

     @return the calibrated spot component power hazard rate combo pointer
   */
  shared_ptr<HazardRateCombo> 
  Calibrate(const finance::TermStructureDerivative& tsDeriv,
            const shared_ptr<Volatility>& pVolatility,
            size_t nIterMax = 50);

  /**
     Return the hazard rate calibrated in the previous call to Calibrate.
     Even if Calibrate failed, this function should return the best guess
     to the calibrated hazard rate.

     @return the calibrated spot component power hazard rate combo pointer
   */
  shared_ptr<HazardRateCombo> GetHazardRate()
  {
    return m_pHazardRate; 
  } 

private:

  /// Keep a copy of the Newton calibrator so it can store previous guesses
  HazardRateSpotComponentPowerCalibratorNewton calibratorNewton;

  /// The final calibrated hazard rate
  shared_ptr<HazardRateCombo> m_pHazardRate;

}; 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_MIXED_H_

