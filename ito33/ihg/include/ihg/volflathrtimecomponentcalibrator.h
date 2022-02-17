/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volflathrtimecomponentcalibrator.h
// Purpose:     calibrator on a vol flat and a hr with time component
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: volflathrtimecomponentcalibrator.h,v 1.11 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/volflathrtimecomponentcalibrator.h
    @brief calibrator on a vol flat and a hr with time component
 */

#ifndef _ITO33_IHG_VOLFLATHRTIMECOMPONENTCALIBRATOR_H_
#define _ITO33_IHG_VOLFLATHRTIMECOMPONENTCALIBRATOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/common.h"

namespace ito33
{

namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityFlat;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;

/**
   calibrator on a vol flat and a hr with time component
 */
class VolFlatHRTimeComponentCalibrator
{
public:

  /// Default ctor and dtor is ok

  /**
     Calibrates with a derivative and a term structure of derivative

     @param derivative a derivative to help calibrate a vol flat
     @param tsDerivative a term structure of derivative
     @param pHazardRate the known hazard rate of the model
   */
  void Calibrate(const finance::Derivative& derivative, 
                 const finance::TermStructureDerivative& tsDerivative,
                 shared_ptr<HazardRateWithTimeComponent> pHazardRate);

  /**
     Gets the calibrated flat volatility.

     @return a shared pointer to the calibrated flat volatility
   */
  shared_ptr<VolatilityFlat> GetVolatility() const { return m_pVolatility; }

  
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
  bool IsConverged(const finance::Derivative& derivative, 
                   const finance::TermStructureDerivative& tsDerivative);

  /// The calibrated flat volatility
  shared_ptr<VolatilityFlat> m_pVolatility;

  /// The calibrated hazard rate with time component
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

  /// Cached market price of the derivative to calibrate
  double m_dDerivMarketPrice;

  /// Cached market prices of the term structure to calibrate
  std::vector<double> m_pdTSMarketPrices;

}; // class VolFlatHRTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLFLATHRTIMECOMPONENTCALIBRATOR_H_
