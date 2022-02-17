/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltanhhrtimecomponentcalibrator.h
// Purpose:     calibrator on a tanh vol and a hr with time component
// Author:      ITO 33
// Created:     2005/01/05
// RCS-ID:      $Id: voltanhhrtimecomponentcalibrator.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/voltanhhrtimecomponentcalibrator.h
    @brief calibrator on a vol tanh and a hr with time component
 */

#ifndef _IHG_VOLTANHHRTIMECOMPONENTCALIBRATOR_H_
#define _IHG_VOLTANHHRTIMECOMPONENTCALIBRATOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"

namespace ito33
{

namespace ihg
{
  class SpotComponent;
  class VolatilityTanh;
  class HazardRateWithTimeComponent;

/**
   calibrator on a vol tanh and a hr with time component
 */
class VolTanhHRTimeComponentCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with two derivatives and a term structure of derivative

     @param derivative1 first derivative to help calibrate a tanh vol
     @param derivative2 second derivative to help calibrate a tanh vol
     @param dScale The scale to the tanh function
     @param dS0 The shift to the tanh function
     @param tsDerivative a term structure of derivative
     @param pHazardRate the known spot component of the hazard rate
   */
  void Calibrate(const finance::Derivative& derivative1, 
                 const finance::Derivative& derivative2,
                 double dScale,
                 double dS0,
                 const finance::TermStructureDerivative& tsDerivative,
                 shared_ptr<HazardRateWithTimeComponent> pHazardRate);

  /**
     Gets the calibrated tanh volatility.

     @return a shared pointer to the calibrated tanh volatility
   */
  shared_ptr<VolatilityTanh> GetVolatility() const { return m_pVolatility; }

  
  /**
     Gets the calibrated hazard rate with time component

     @return a shared pointer to the calibrated hazard rate with time component
   */
  shared_ptr<HazardRateWithTimeComponent> GetHazardRate() const
  {
    return m_pHazardRate;
  }


private:

  /*
    Helper function to determine if the iterative scheme has converged

    @return true if converged, false otherwise
  */
  bool IsConverged(const finance::Derivative& derivative1, 
                   const finance::Derivative& derivative2, 
                   const finance::TermStructureDerivative& tsDerivative);

  /// The calibrated tanh volatility
  shared_ptr<VolatilityTanh> m_pVolatility;

  /// The calibrated hazard rate with time component
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

  /// Cached market price of the first derivative to calibrate
  double m_dDeriv1MarketPrice;

  /// Cached market price of the 2nd derivative to calibrate
  double m_dDeriv2MarketPrice;

  /// Cached market prices of the term structure to calibrate
  std::vector<double> m_pdTSMarketPrices;

}; // class VolTanhHRTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_VOLTANHHRTIMECOMPONENTCALIBRATOR_H_
