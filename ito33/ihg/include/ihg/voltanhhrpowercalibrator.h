/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltanhhrpowercalibrator.h
// Purpose:     calibrator on a tanh vol and a hr with spot component tanh
// Author:      ITO33
// Created:     2005/01/03
// RCS-ID:      $Id: voltanhhrpowercalibrator.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/voltanhhrpowercalibrator.h
    @brief calibrator on a tanh flat and a hr with spot component tanh
 */

#ifndef _ITO33_IHG_VOLTANHHRPOWERCALIBRATOR_H_
#define _ITO33_IHG_VOLTANHHRPOWERCALIBRATOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"

namespace ito33
{

 namespace finance
{
  class Derivative;
}


namespace ihg
{
  class SpotComponent;
  class VolatilityTanh;
  class HazardRateCombo;
  class HRSpotComponentTanh;

/**
   Calibrator on a vol tanh and a hr with spot component tanh. When
   calibration fails, the last guesses to the vol and hr are saved and
   can still be retrieved.
 */
class VolTanhHRPowerCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with two derivatives and a derivative term structure

     @param deriv1 the first derivative to help calibrate a tanh vol
     @param deriv1 the second derivative to help calibrate a tanh vol
     @param dScale the scale to the tanh function of volatility
     @param dS0 shift to the tanh function of volatility
     @param tsDerivs a term structure to help calibrate the tanh hazard rate
   */
  void Calibrate(const finance::Derivative& deriv1,
                 const finance::Derivative& deriv2, 
                 double dScale,
                 double dS0,
                 const finance::TermStructureDerivative& tsDerivs);

  /**
     Gets the calibrated tanh volatility.

     @return a shared pointer to the calibrated tanh volatility
   */
  shared_ptr<VolatilityTanh> GetVolatility() const { return m_pVolatility; }

  
  /**
     Gets the calibrated hazard rate with spot and time component

     @return shared pointer to calibrated hazard rate with spot/time component
   */
  shared_ptr<HazardRateCombo> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /*
    Gets the hazard rate spot component power calibrated to the first 
    and last derivative of the term structure

    @return a shared pointer to the calibrated hazard rate component tanh
  */
  shared_ptr<HRSpotComponentPower> GetHRSpotComponentPower() const
  {
    return m_pHRSpotComponentPower;
  }

private:

  /*
    Helper function to determine if the iterative scheme has converged

    @return true if converged, false otherwise
  */
  bool IsConverged(const finance::Derivative& deriv1,
                   const finance::Derivative& deriv2, 
                   const finance::TermStructureDerivative& tsDeriv);

  /// The calibrated tanh volatility
  shared_ptr<VolatilityTanh> m_pVolatility;

  /// The calibrated hazard rate with spot and time component
  shared_ptr<HazardRateCombo> m_pHazardRate;

  /// The calibrated hazard rate spot component power
  shared_ptr<HRSpotComponentPower> m_pHRSpotComponentPower;

  /// Cached market price of the first derivative to calibrate
  double m_dDeriv1MarketPrice;

  /// Cached market price of the 2nd derivative to calibrate
  double m_dDeriv2MarketPrice;

  /// Cached market prices of the term structure to calibrate
  std::vector<double> m_pdTSMarketPrices;

}; // class VolTanhHRPowerCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLTANHHRPOWERCALIBRATOR_H_
