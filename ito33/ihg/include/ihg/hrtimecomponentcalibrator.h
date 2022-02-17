/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hrtimecomponentcalibrator.h
// Purpose:     calibrate time component of HazardRateWithtimeComponent
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: hrtimecomponentcalibrator.h,v 1.9 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/hrtimecomponentcalibrator.h
   @brief calibration of time component of HazardRateWithtimeComponent.
 */


#ifndef _ITO33_IHG_HRTIMECOMPONENTCALIBRATOR_H_
#define _ITO33_IHG_HRTIMECOMPONENTCALIBRATOR_H_

#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/finance/termstructurederivative.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;


/// class for time only hazard rate calibration using a term structure
class HazardRateTimeComponentCalibrator 
{

public:

  HazardRateTimeComponentCalibrator()
  {
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates a term structure with known volatility

     @param pVolatility the known volatility of the model
     @param tsDerivative a term structure of derivative
     @param pHazardRate the (partially) given hazardrate

     @return the calibrated time only hazard rate pointer
   */
  shared_ptr<HazardRateWithTimeComponent> 
    Calibrate(const finance::TermStructureDerivative& tsDerivative,
            const shared_ptr<Volatility>& pVolatility, 
            shared_ptr<HazardRateWithTimeComponent> pHazardRate =
                shared_ptr<HazardRateWithTimeComponent>()
           );

  /**
     Helper operator to be used for numerical calibrator

     @return the difference of the theoritical price and the market price
             for the current derivative
   */
  double operator()(double dHazardRate);


private:

  /// Helper Model used during calibration. This class updates the hazard rate
  TheoreticalModel m_theoreticalModel;

  /// the temporary hazard rate pointer
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

  /// the market price of the current derivative to be calibrated  
  double m_dMarketPrice;

  /// a pointer to the current derivative
  shared_ptr<finance::Derivative> m_pDerivative;

  /// The index
  size_t m_nIdx;

  /// Helper for intermediary results
  std::vector<Date> m_pMaturityDates;

  std::vector<double> m_pdValues;

}; // class HazardRateTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HRTIMECOMPONENTCALIBRATOR_H_

