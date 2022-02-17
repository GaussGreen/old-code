/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hrspotcomponentpowercalibrator_newton.h
// Purpose:     calibrate spot component power of HRSpotComponentPower
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: hrspotcomponentpowercalibrator_newton.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/hrspotcomponentpowercalibrator_newton.h
   @brief calibrate spot component power of HRSpotComponentPower by using 
          newton iteration
 */

#ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_NEWTON_H_
#define _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_NEWTON_H_

#include "ito33/finance/termstructurederivative.h"
#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

  class Volatility;
  class HazardRateCombo;


/// class for spot power hazard rate calibration using a term structure
class HazardRateSpotComponentPowerCalibratorNewton
{

public:

  //default constructor
  HazardRateSpotComponentPowerCalibratorNewton()
    : m_dAlpha(0.1),
      m_dBeta(0.2)
  {
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates a term structure with known volatility

     @param tsDeriv the derivative term structure to calibrate
     @param pVolatility the known volatility
     @param nIterMax maximum number of iterations
     @param bUsePreviousSolution use the previously stored solution or not

     @return the calibrated spot component power hazard rate pointer
   */
  shared_ptr<HazardRateCombo> 
  Calibrate(const finance::TermStructureDerivative& tsDeriv,
            const shared_ptr<Volatility>& pVolatility,
            size_t nIterMax = 15,
            bool bUsePreviousSolution = true);


  /** 
      Used by the non-linear newton algorithm
      to find the root for the hazard rate.

      @param dAlpha The current guess for alpha 
      @param dBeta The current guess for beta
      @param dPrice1 price of the first instrument
      @price dPrice2 price of the second instrument
      @param dF (output) norm2 of the solution
   */
  void operator () (double dAlpha, double dBeta, 
                    double &dPrice1,double &dPrice2,
                    double &dF);

  /**
     Return the hazard rate calibrated in the previous call to Calibrate.
     Even if Calibrate failed, this function should return the best guess
     to the calibrated hazard rate.

     @return the calibrated spot component power hazard rate combo pointer
   */
  shared_ptr<HazardRateCombo> GetHazardRate();

private:

  /// Helper Model used during calibration. This class updates the hazard rate
  ihg::TheoreticalModel m_theoreticalModel;

   /// first derivative used for spot component power calibration
  const finance::Derivative* m_pDeriv1;

  /// second derivative used for spot component power calibration
  const finance::Derivative* m_pDeriv2;

  ///market price of first instrument
  double m_dMarketPrice1;

  ///market price of second instrument
  double m_dMarketPrice2;

  /// scale value to get relative error of first instrument
  double m_dScale1;

  /// scale value to get relative error of second instrument
  double m_dScale2;

  ///Spot price of all instrument
  double m_dS0;

  ///volatility
  shared_ptr<Volatility> m_pVol;

  //store value of alpha for future use
  double m_dAlpha;

  //store value of beta for future use
  double m_dBeta;

  /// Flag indicating if the calibration worked
  bool m_bConverged;

  /// The smallest objective function obtained
  double m_dObjectif;

  /// The final calibrated hazard rate
  shared_ptr<HazardRateCombo> m_pHazardRate;

  /// The term structure to calibrate
  finance::TermStructureDerivative m_tsDerivative;

}; 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_NEWTON_H_

