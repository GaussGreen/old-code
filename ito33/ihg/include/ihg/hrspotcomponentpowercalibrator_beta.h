/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hrspotcomponentpowercalibrator_beta.h
// Purpose:     calibrate spot component power of HRSpotComponentPower
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_beta.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/hrspotcomponentpowercalibrator_beta.h
   @brief calibrate spot component power of HRSpotComponentPower
 */

#ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_BETA_H_
#define _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_BETA_H_

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

  class Volatility;
  class HazardRateCombo;

/// class for spot power hazard rate calibration using a term structure
class HazardRateSpotComponentPowerCalibratorBeta
{

public:

  HazardRateSpotComponentPowerCalibratorBeta()
  {
    m_model.SetExternalFlagsToDefaults();
  }

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
     Objective function, for a minimization routine
   */
  double operator()(double dBeta);

  /**
     Return the hazard rate calibrated in the previous call to Calibrate.
     Even if Calibrate failed, this function should return the best guess
     to the calibrated hazard rate.

     @return the calibrated spot component power hazard rate combo pointer
   */
  shared_ptr<HazardRateCombo> GetHazardRate(); 

private:

  /// Helper Model used during calibration. This class updates the hazard rate
  TheoreticalModel m_model;

  /// First derivative contract used for spot component power calibration
  const finance::Derivative* m_pDeriv1;

  /// Second derivative contract used for spot component power calibration
  const finance::Derivative* m_pDeriv2;

  /// Cached market price of first derivative
  double m_dDeriv1MarketPrice;

  /// Cached market price of 2nd derivative
  double m_dDeriv2MarketPrice;

  /// Current spot price
  double m_dS0;

  /// Scale used to compute relative error
  double m_dScale;

  /// The volatility used when calibrating the hr spot component
  shared_ptr<Volatility> m_pVolatility;

  /// The final (or last guess before failure) alpha value
  double m_dAlpha;

  /// The final (or last guess before failure) beta value
  double m_dBeta;

  /// The final (or last guess before failure) objective function value
  double m_dObjectif;

  /// Flag indicating if the calibration worked
  bool m_bConverged;

  /// The final calibrated hazard rate
  shared_ptr<HazardRateCombo> m_pHazardRate;

  /// The term structure to calibrate
  finance::TermStructureDerivative m_tsDerivative;

}; 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_BETA_H_

