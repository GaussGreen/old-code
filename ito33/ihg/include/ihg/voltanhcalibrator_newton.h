/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltanhcalibrator_newton.h
// Purpose:     calibrate Tanh volatility
// Created:     2005/02/04
// RCS-ID:      $Id: voltanhcalibrator_newton.h,v 1.7 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/voltanhcalibrator_newton.h
   @brief calibrate tanh volatility by using Newton iteration
 */

#ifndef _ITO33_IHG_VOLTANHCALIBRATOR_NEWTON_H_
#define _ITO33_IHG_VOLTANHCALIBRATOR_NEWTON_H_

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}

namespace ihg
{

  class ITO33_IHG_DLLDECL VolatilityTanh;
  class ITO33_IHG_DLLDECL HazardRate;


/**
   Class for volatility calibration where volatility is a tanh function of
   spot.
 */
class VolTanhCalibratorNewton
{

public:

  // default constructor
  VolTanhCalibratorNewton(double dScale, double dS0) 
                        : m_dLeft(0.0), m_dRight(0.0),
                          m_dScale(dScale), m_dS0(dS0)
  {
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates tanh volatility to two derivatives.

     @param deriv1 first derivative used for calibration
     @param deriv2 second derivative used for calibration
     @param pHazardRate the known hazard rate
     @param nIterMax maximum number of iterations
     @param bUsePreviousSolution use the previously stored solution or not

     @return the calibrated volatility
   */
  shared_ptr<VolatilityTanh> 
  Calibrate(const finance::Derivative& deriv1,
            const finance::Derivative& deriv2,
            const shared_ptr<HazardRate>& pHazardRate,
            size_t nIterMax = 15,
            bool bUsePreviousSolution = true);


  /** 
      Used by the non-linear newton algorithm, evaluates the function values.

      @param dLeft The current guess for left limit 
      @param dRight The current guess for right limit
      @param dPrice1 (output) computed price of the first instrument
      @param dPrice2 (output) computed price of the second instrument
      @param dF (output) objective function value
   */
  void operator () (double dLeft, double dLower, 
                    double& dPrice1,double& dPrice2,
                    double& dF);

  /**
     Return the volatility calibrated in the previous call to Calibrate.
     Even if Calibrate failed, this function should return the best guess
     to the calibrated volatility.

     @return the calibrated tanh volatility
   */
  shared_ptr<VolatilityTanh> GetVolatility();


private:

  /// Helper function to find a good initial guess
  //void FindInitialGuess();

  /// Helper Model used during calibration. This class updates the hazard rate
  TheoreticalModel m_theoreticalModel;

   /// first derivative used for calibration
  const finance::Derivative* m_pDeriv1;

  /// second derivative used for calibration
  const finance::Derivative* m_pDeriv2;

  /// market price of first instrument
  double m_dMarketPrice1;

  /// market price of second instrument
  double m_dMarketPrice2;

  /// scale value to get relative error of first instrument
  double m_dScale1;

  /// scale value to get relative error of second instrument
  double m_dScale2;

  /// Scale to the tanh function
  double m_dScale;

  /// Shift to the tanh function
  double m_dS0;

  /// hazard rate used for pricing
  shared_ptr<HazardRate> m_pHazardRate;

  // store value of left limit for future use
  double m_dLeft;

  // store value of right limit for future use
  double m_dRight;

  /// The smallest objective function obtained
  double m_dObjectif;
}; 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLTANHCALIBRATOR_NEWTON_H_
