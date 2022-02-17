/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltimecomponentcalibrator.h
// Purpose:     calibrate time component of VolatilityWithtimeComponent
// Created:     2005/03/04
// RCS-ID:      $Id: voltimecomponentcalibrator.h,v 1.7 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/voltimecomponentcalibrator.h
   @brief calibration of time component of VolatilityWithtimeComponent.
 */

#ifndef _ITO33_IHG_VOLTIMECOMPONENTCALIBRATOR_H_
#define _ITO33_IHG_VOLTIMECOMPONENTCALIBRATOR_H_

#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL VolatilityWithTimeComponent;


/// class for time only hazard rate calibration using a term structure
class VolTimeComponentCalibrator 
{

public:

  /**
     @param pHazardRate the known hazard rate of the model
   */
  VolTimeComponentCalibrator(const shared_ptr<HazardRate>& pHazardRate)
  {
    m_theoreticalModel.SetHazardRate(pHazardRate);
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates a term structure with known volatility

     
     @param tsDerivative a term structure of derivative
     @param pVolatility the (partially) given volatility

     @return the calibrated time only hazard rate pointer
   */
  shared_ptr<VolatilityWithTimeComponent> 
  Calibrate(const finance::TermStructureDerivative& tsDerivative,
            shared_ptr<VolatilityWithTimeComponent> pVolatility =
                shared_ptr<VolatilityWithTimeComponent>());

  /**
     Helper operator to be used for numerical calibrator

     @return the difference of the theoritical price and the market price
             for the current derivative
   */
  double operator()(double dVolatility);


private:

  /// Helper Model used during calibration. This class updates the hazard rate
  ihg::TheoreticalModel m_theoreticalModel;

  /// the temporary volatility pointer
  shared_ptr<VolatilityWithTimeComponent> m_pVolatility;

  /// the market price of the current derivative to be calibrated  
  double m_dMarketPrice;

  /// the scale on the function
  double m_dScale;

  /// a pointer to the current derivative
  shared_ptr<finance::Derivative> m_pDerivative;

  /// The index
  size_t m_nIdx;

  /// Helper for intermediary results
  std::vector<Date> m_pMaturityDates;

  std::vector<double> m_pdValues;

}; // class VolTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLTIMECOMPONENTCALIBRATOR_H_

