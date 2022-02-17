/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volpowercalibrator_newton.h
// Purpose:     calibrate power volatility
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: volpowercalibrator_newton.h,v 1.4 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004-  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/volpowercalibrator_newton.h
   @brief calibrate power volatility by using Newton iteration
 */

#ifndef _ITO33_IHG_VOLPOWERCALIBRATOR_NEWTON_H_
#define _ITO33_IHG_VOLPOWERCALIBRATOR_NEWTON_H_

#include "ito33/finance/derivative.h"
#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

  class VolatilityPower;
  class HazardRate;

/// class for spot power hazard rate calibration using a term structure
class VolPowerCalibratorNewton
{

public:

  // default constructor
  VolPowerCalibratorNewton()
    : m_dAlpha(0.0),
      m_dBeta(0.0)
  {
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates power volatility to two derivatives

     @param deriv1 first derivative used for calibration
     @param deriv2 second derivative used for calibration
     @param pHazardRate the known hazard rate
     @param nIterMax maximum number of iterations
     @param bUsePreviousSolution use the previously stored solution or not

     @return the calibrated power volatility
   */
  shared_ptr<VolatilityPower> 
  Calibrate(const finance::Derivative& deriv1,
            const finance::Derivative& deriv2,
            const shared_ptr<HazardRate>& pHazardRate,
            size_t nIterMax = 15,
            bool bUsePreviousSolution = true);


  /** 
      Used by the non-linear newton algorithm
      to find the root for the hazard rate.

      @param dAlpha The current guess for alpha 
      @param dBeta The current guess for beta
      @param dPrice1 (output) computed price of the first instrument
      @param dPrice2 (output) computed price of the second instrument
      @param dF (output) objective function value
      @param dPrice1Error (output) error of first computed price
      @param dPrice2Error (output) error of second computed price
   */
  void operator () (double dAlpha, double dBeta, 
                    double &dPrice1,double &dPrice2,
                    double &dF);

  /**
     Return the volatility calibrated in the previous call to Calibrate.
     Even if Calibrate failed, this function should return the best guess
     to the calibrated volatility.

     @return the calibrated power volatility
   */
  shared_ptr<VolatilityPower> GetVolatility();

private:

  /// Helper Model used during calibration. This class updates the hazard rate
  ihg::TheoreticalModel m_theoreticalModel;

   /// first derivative used for spot component power calibration
  const finance::Derivative* m_pDeriv1;

  /// second derivative used for spot component power calibration
  const finance::Derivative* m_pDeriv2;

  /// market price of first instrument
  double m_dMarketPrice1;

  /// market price of second instrument
  double m_dMarketPrice2;

  /// scale value to get relative error of first instrument
  double m_dScale1;

  /// scale value to get relative error of second instrument
  double m_dScale2;

  /// Spot value
  double m_dS0;

  /// hazard rate used for pricing
  shared_ptr<HazardRate> m_pHazardRate;

  // store value of alpha for future use
  double m_dAlpha;

  // store value of beta for future use
  double m_dBeta;

  /// The smallest objective function obtained
  double m_dObjectif;
}; 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLPOWERCALIBRATOR_NEWTON_H_

