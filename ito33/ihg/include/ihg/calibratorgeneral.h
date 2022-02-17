/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/calibratorgeneral.h
// Purpose:     General calibration class for IHG
// Created:     2005/07/05
// RCS-ID:      $Id: calibratorgeneral.h,v 1.7 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/calibratorgeneral.h
    @brief General calibration class for IHG.
 */

#ifndef _IHG_CALIBRATORGENERAL_H_
#define _IHG_CALIBRATORGENERAL_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/calibratorgeneral.h"
#include "ito33/ihg/common.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
}

namespace ihg
{

class Translator;
class ITO33_IHG_DLLDECL TheoreticalModel;

/**
    General calibration class.

    Fits volatility and hazard rate to an arbitrary list of contracts. It is
    assumed that all contracts in the list can be priced by the ihg model.
 */
class CalibratorGeneral : public pricing::CalibratorGeneral
{

public:

  CalibratorGeneral(finance::CalibrationProcess& calibrationProcess)
    : pricing::CalibratorGeneral(calibrationProcess)
  {
  }

  /**
      Calibrates the specified derivatives using the specified process
      structure.

      @param pTranslator Translator from underlying process to parameter arrays
      @param derivatives The collection of derivatives for price matching

      @return the calibrated model
   */
  shared_ptr<ihg::TheoreticalModel>
  Calibrate(Translator* pTranslator, const finance::Derivatives& derivatives);

  /**
      Returns the last calibrated underlying process (or the last best guess).

      Useful if the calibration fails, but the best calibrated data
      is still desired.

      @return the last calibrated model
   */
  shared_ptr<TheoreticalModel> GetLastCalibratedProcess() const;
  
  virtual void DoComputeObjectif(const double *pdX, double& dObjectif,
                                 double* pdGradients, bool bComputeGradient);
  
  virtual void operator () (const std::vector<double>& pdParam, double& dF);

protected:

  virtual void RunCalibrator(bool bUserGradients);

  /// The translator 
  Translator* m_pTranslator;

}; // class CalibratorGeneral


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_CALIBRATORGENERAL_H_
