/////////////////////////////////////////////////////////////////////////////
// Name:        hg/calibratorgeneral.h
// Purpose:     General calibration class for HG
// Created:     2005/07/04
// RCS-ID:      $Id: calibratorgeneral.h,v 1.23 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/calibratorgeneralhg.h
    @brief General calibration class for HG
 */

#ifndef _HG_CALIBRATORGENERAL_H_
#define _HG_CALIBRATORGENERAL_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/calibratorgeneral.h"


namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
}


namespace hg
{

class ITO33_HG_DLLDECL UnderlyingProcess;
class Translator;

/**
    General calibration class.

    Fits a regime/jump structure to an arbitrary list of contracts. It is 
    assumed that all contracts in the list can be priced by the hg model.
 */
class CalibratorGeneral : public pricing::CalibratorGeneral
{

public:

  /**
      Ctor uses by default finite element for space discretization since it
      gives better sensitivities.
   */
  CalibratorGeneral(finance::CalibrationProcess& calibrationProcess)
    : pricing::CalibratorGeneral(calibrationProcess),
      m_iDiscretizationMethod(1)
  {
  }

  /**
      Calibrates the specified derivatives using the specified underlying
      process structure.

      @param pTranslator Translator from underlying process to parameter arrays
      @param derivatives The collection of derivatives for price matching

      @return the calibrated underlying process
   */
  shared_ptr<UnderlyingProcess>
  Calibrate(Translator* pTranslator, const finance::Derivatives& derivatives);

  /**
      Returns the last calibrated underlying process (or the last best guess).

      Useful if the calibration fails, but the best calibrated data
      is still desired.

      @return the last calibrated process
   */
  shared_ptr<UnderlyingProcess> GetLastCalibratedProcess() const;
  
  virtual void DoComputeObjectif(const double *pdX, double& dObjectif,
                                 double* pdGradients, bool bComputeGradient);

  virtual void operator () (const std::vector<double>& pdParam, double& dF);


protected:

  virtual void RunCalibrator(bool bUserGradients);

  /// The translator from parameter array into underlying process 
  Translator* m_pTranslator;

  /// The (constant) underlying process used for mesh construction
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcessForMesh;

  /// discretization method: 0 = FD, 1 = FE
  int m_iDiscretizationMethod;

}; // class CalibratorGeneral


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CALIBRATORGENERAL_H_
