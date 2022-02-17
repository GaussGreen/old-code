/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltimeonlyhrtimecomponentcalibrator.h
// Purpose:     calibrator on a timeonly vol and a hr with time component
// Created:     2005/07/15
// RCS-ID:      $Id: voltimeonlyhrtimecomponentcalibrator.h,v 1.4 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/voltimeonlyhrtimecomponentcalibrator.h
    @brief calibrator on a time only vol and a hr with time component
 */

#ifndef _ITO33_IHG_VOLTIMEONLYHRTIMECOMPONENTCALIBRATOR_H_
#define _ITO33_IHG_VOLTIMEONLYHRTIMECOMPONENTCALIBRATOR_H_

#include "ito33/list.h"
#include "ito33/sharedptr.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}

namespace ihg
{

  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;
  class ITO33_IHG_DLLDECL VolatilityTimeOnly;

/**
   calibrator on a time only vol and a hr with time component
 */
class VolTimeOnlyHRTimeComponentCalibrator
{
public:

  VolTimeOnlyHRTimeComponentCalibrator()
  {
    m_theoreticalModel.SetExternalFlagsToDefaults();
  }

  /**
     Calibrates with a general derivative list.

     The derivative list is assumed to have an even number of entries,
     to be sorted by maturity, and that 3 (or more) derivatives do not
     have the same maturity date.

     @param derivatives a general derivative list
     @param pHazardRate the known hazard rate of the model (can be NULL)
   */
  void Calibrate(const std::list< shared_ptr<finance::Derivative> >& derivatives, 
                 shared_ptr<HazardRateWithTimeComponent> pHazardRate);

  /**
     Gets the calibrated time only volatility.

     @return a shared pointer to the calibrated time only volatility
   */
  shared_ptr<VolatilityTimeOnly> GetVolatility() const 
  { 
    return m_pVolatility; 
  }
  
  /**
     Gets the calibrated hazard rate with time component

     @return a shared pointer to the calibrated hazard rate with time component
   */
  shared_ptr<HazardRateWithTimeComponent> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /**
     Operator called by Newton's method.

     @param dVol The current time component volatility value
     @param dHR The current time component hazard rate value
     @param dPrice1 The price (minus market price) of first derivative (output)
     @param dPrice2 The price (minus market price) of 2nd derivative (output)
     @param dF The objective function value (output)
   */
  void operator()(double dVol, double dHR, 
                  double &dPrice1, double &dPrice2, 
                  double &dF);


private:

  /// The calibrated flat volatility
  shared_ptr<VolatilityTimeOnly> m_pVolatility;

  /// The calibrated hazard rate with time component
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

  /// The time component (vol/hr) maturity dates
  std::vector<Date> m_pMaturityDates;

  /// The time component vol values
  std::vector<double> m_pdVolValues;

  /// The time component hr values
  std::vector<double> m_pdHRValues;

  /// The current index into the time component vectors
  size_t m_nIdx;

  /// The first derivative in the current calibration period
  shared_ptr<finance::Derivative> m_pDeriv1;

  /// The second derivative in the current calibration period
  shared_ptr<finance::Derivative> m_pDeriv2;

  /// Helper model used during calibration. 
  TheoreticalModel m_theoreticalModel;

  /// The market price of deriv1
  double m_dMarketPrice1;

  /// The market price of deriv2
  double m_dMarketPrice2;

  /// The relative error scale factor for deriv1
  double m_dScale1;

  /// The relative error scale factor for deriv1
  double m_dScale2;

}; // class VolTimeOnlyHRTimeComponentCalibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLTIMEONLYHRTIMECOMPONENTCALIBRATOR_H_
