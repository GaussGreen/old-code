/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volflathrspotcomponentpowercalibrator.h
// Purpose:     calibrator on a vol flat and a hr with spot component power
// Author:      ITO33
// Created:     2004/11/23
// RCS-ID:      $Id: volflathrspotcomponentpowercalibrator.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/volflathrspotcomponentpowercalibrator.h
    @brief calibrator on a vol flat and a hr with spot component power
 */

#ifndef _ITO33_IHG_VOLFLATHRSPOTCOMPONENTPOWERCALIBRATOR_H_
#define _ITO33_IHG_VOLFLATHRSPOTCOMPONENTPOWERCALIBRATOR_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/common.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}


namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityFlat;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HRSpotComponentPower;

/**
   Calibrator on a vol flat and a hr with spot component power. When
   calibration fails, the last guesses to the vol and hr are saved and
   can still be retrieved.
 */
class VolFlatHRSpotComponentPowerCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with a derivative and a term structure

     @param derivative the derivative to help calibrate a flat vol
     @param tsDerivative a term structure to help calibrate the hazard rate
   */
  void Calibrate(const finance::Derivative& derivative, 
                 const finance::TermStructureDerivative& tsDerivative);

  /**
     Gets the calibrated flat volatility.

     @return a shared pointer to the calibrated flat volatility
   */
  shared_ptr<VolatilityFlat> GetVolatility() const { return m_pVolatility; }

  
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
  bool IsConverged(const finance::Derivative& derivative, 
                   const finance::TermStructureDerivative& tsDerivative);

  /// The calibrated flat volatility
  shared_ptr<VolatilityFlat> m_pVolatility;

  /// The calibrated hazard rate with spot and time component
  shared_ptr<HazardRateCombo> m_pHazardRate;

  /// The calibrated hazard rate spot component power
  shared_ptr<HRSpotComponentPower> m_pHRSpotComponentPower;

  /// Cached market price of the derivative to calibrate
  double m_dDerivMarketPrice;

  /// Cached market prices of the term structure to calibrate
  std::vector<double> m_pdTSMarketPrices;

}; // class VolFlatHRSpotComponentPowerCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLFLATHRTIMECOMPONENTCALIBRATOR_H_
