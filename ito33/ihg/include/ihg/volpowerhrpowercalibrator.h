/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volpowerhrpowercalibrator.h
// Purpose:     calibrator on a power vol and a hr with spot component power
// Author:      ITO33
// Created:     2005/01/03
// RCS-ID:      $Id: volpowerhrpowercalibrator.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/volpowerhrpowercalibrator.h
    @brief calibrator on a power flat and a hr with spot component power
 */

#ifndef _ITO33_IHG_VOLPOWERHRPOWERCALIBRATOR_H_
#define _ITO33_IHG_VOLPOWERHRPOWERCALIBRATOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"

namespace ito33
{

 namespace finance
{
  class ITO33_DLLDECL Derivative;
}


namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityPower;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HRSpotComponentPower;

/**
   Calibrator on a vol power and a hr with spot component power. When
   calibration fails, the last guesses to the vol and hr are saved and
   can still be retrieved.
 */
class VolPowerHRPowerCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with two derivatives and a derivative term structure

     @param deriv1 the first derivative to help calibrate a power vol
     @param deriv1 the second derivative to help calibrate a power vol
     @param tsDerivs a term structure to help calibrate the power hazard rate
   */
  void Calibrate(const finance::Derivative& deriv1,
                 const finance::Derivative& deriv2, 
                 const finance::TermStructureDerivative& tsDerivs);

  /**
     Gets the calibrated power volatility.

     @return a shared pointer to the calibrated power volatility
   */
  shared_ptr<VolatilityPower> GetVolatility() const { return m_pVolatility; }

  
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

    @return a shared pointer to the calibrated hazard rate component power
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

  /// The calibrated power volatility
  shared_ptr<VolatilityPower> m_pVolatility;

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

}; // class VolPowerHRPowerCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLPOWERHRPOWERCALIBRATOR_H_
