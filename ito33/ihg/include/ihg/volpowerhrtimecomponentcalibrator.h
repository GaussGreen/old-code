/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volpowerhrtimecomponentcalibrator.h
// Purpose:     calibrator on a power vol and a hr with time component
// Author:      ITO 33
// Created:     2005/01/05
// RCS-ID:      $Id: volpowerhrtimecomponentcalibrator.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/volpowerhrtimecomponentcalibrator.h
    @brief calibrator on a vol power and a hr with time component
 */

#ifndef _ITO33_IHG_VOLPOWERHRTIMECOMPONENTCALIBRATOR_H_
#define _ITO33_IHG_VOLPOWERHRTIMECOMPONENTCALIBRATOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/common.h"

namespace ito33
{

namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityPower;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;

/**
   calibrator on a vol power and a hr with time component
 */
class VolPowerHRTimeComponentCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with two derivatives and a term structure of derivative

     @param derivative1 first derivative to help calibrate a power vol
     @param derivative2 second derivative to help calibrate a power vol
     @param tsDerivative a term structure of derivative
     @param pHazardRate the known spot component of the hazard rate
   */
  void Calibrate(const finance::Derivative& derivative1, 
                 const finance::Derivative& derivative2, 
                 const finance::TermStructureDerivative& tsDerivative,
                 shared_ptr<HazardRateWithTimeComponent> pHazardRate);

  /**
     Gets the calibrated power volatility.

     @return a shared pointer to the calibrated power volatility
   */
  shared_ptr<VolatilityPower> GetVolatility() const { return m_pVolatility; }

  
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

  /// The calibrated power volatility
  shared_ptr<VolatilityPower> m_pVolatility;

  /// The calibrated hazard rate with time component
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

  /// Cached market price of the first derivative to calibrate
  double m_dDeriv1MarketPrice;

  /// Cached market price of the 2nd derivative to calibrate
  double m_dDeriv2MarketPrice;

  /// Cached market prices of the term structure to calibrate
  std::vector<double> m_pdTSMarketPrices;

}; // class VolPowerHRTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLPOWERHRTIMECOMPONENTCALIBRATOR_H_
