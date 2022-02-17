/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/calibratorgeneral.h
// Purpose:     General calibration base class for all models
// Created:     2005/07/04
// RCS-ID:      $Id: calibratorgeneral.h,v 1.13 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/calibratorgeneral.h
    @brief General calibration base class for all models.
 */

#ifndef _ITO33_PRICING_CALIBRATORGENERAL_H_
#define _ITO33_PRICING_CALIBRATORGENERAL_H_

#include "ito33/common.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/calibrationprocess.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
  class ITO33_DLLDECL Derivative;
  class ITO33_DLLDECL Option;

  class ForwardOption;
}


namespace pricing
{

class Translator;

/**
    General calibration class.

    Base class for fitting model parameters to an arbitrary list of derivatives.
 */
class CalibratorGeneral
{

public:

  /**
      Ctor takes a calibration process which can be used to cancel the
      calibration.

      @param calibrationProcess calibration process which can be used to
                                cancel the calibration.
   */
  CalibratorGeneral(finance::CalibrationProcess& calibrationProcess)
  : m_dObjectif(1.e99), m_nNbUnknowns(0),
    m_calibrationProcess(calibrationProcess)
  {
  }

  /// Virtual dtor for base class
  virtual ~CalibratorGeneral() { }

  /**
      Objective function, for a minimization routine. The function
      also checks whether the calibration process has been cancelled
      by the user. In this case, it will throw and stop the calibration.

      @param pdX the current guess for the minimum point
      @param dObjectif (output) the value of the objective function
      @param pdGradients (output) the gradient of the objective function
      @param bComputeGradient whether or not to compute the gradient
   */
  void ComputeObjectif(const double *pdX, double& dObjectif,
                       double* pdGradients, bool bComputeGradient)
  {
    DoComputeObjectif(pdX, dObjectif, pdGradients, bComputeGradient);

    if ( m_calibrationProcess.IsCancelled() )
      ThrowCalibrationCancelledException();
  }

  /**
      The real objective function, for a minimization routine that
      should be implemented in sub-classes.

      @param pdX the current guess for the minimum point
      @param dObjectif (output) the value of the objective function
      @param pdGradients (output) the gradient of the objective function
      @param bComputeGradient whether or not to compute the gradient
   */
  virtual void DoComputeObjectif(const double *pdX, double& dObjectif,
                               double* pdGradients, bool bComputeGradient) = 0;

  /**
      Bracket operator needed by simulated annealing.

      @param pdParam parameter list
      @dF Computed objective function value (output)
   */
  virtual void operator () (const std::vector<double>& pdParam,double& dF) = 0;


protected:

  /**
      Initializes the calibration arrays (initial guess, bounds, etc).

      @param translator the translator to get model params
   */
  virtual void InitializeArrays(const Translator& translator);

  /**
      Constructs the derivative lists for pricing.

      @param derivatives the derivative container
      @param bUseForward whether or not to use a forward equation for options
   */
  void ConstructDerivativeLists(const finance::Derivatives& derivatives, 
                                bool bUseForward);

  /**
      Runs the underlying calibration code (NAG, simulated annealing, etc).

      Virtual in base class so that different models can use different
      methods (eg. hg uses nag, ihg uses simulated annealing)

      @param bUserGradient does the objective function compute gradients
   */
  virtual void RunCalibrator(bool bUserGradients) = 0;

  /// Throw exception when the calibration is cancelled by user
  void ThrowCalibrationCancelledException();

protected:

  /// Forward Option
  shared_ptr<finance::ForwardOption> m_pForwardOption;

  /// General list of derivatives to calibrate
  shared_ptr<finance::Derivatives> m_pDerivatives;

  /// The best objective function value so far
  double m_dObjectif;

  /// The initial guess, and vector passed to the calibrator
  std::vector<double> m_pdX;

  /// The best guess so far
  std::vector<double> m_pdFinalX;

  /// The number of variables used in the minimization
  size_t m_nNbUnknowns;

  /// The lower bounds
  std::vector<double> m_pdLowerBounds;

  /// The upper bounds
  std::vector<double> m_pdUpperBounds;

private:

  finance::CalibrationProcess& m_calibrationProcess;

  NO_COPY_CLASS(CalibratorGeneral);

}; // class CalibratorGeneral


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_CALIBRATORGENERAL_H_
