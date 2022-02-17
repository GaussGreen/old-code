/////////////////////////////////////////////////////////////////////////////
// Name:        hg/calibrator_cdsrecovery.h
// Purpose:     General calibration class for HG
// Created:     2005/07/24
// RCS-ID:      $Id: calibrator_cdsrecovery.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/calibrator_cdsrecovery.h
   @brief General calibration class for HG
 */

#ifndef _HG_CALIBRATOR_CDSRECOVERY_H_
#define _HG_CALIBRATOR_CDSRECOVERY_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/calibratorgeneral.h"

namespace ito33
{


namespace hg
{

class ITO33_HG_DLLDECL UnderlyingProcess;
class Translator;
class TheoreticalModel;

/**
  General calibration class

  Fits a regime/jump structure to an arbitrary list of contracts. It is 
  assumed that all contracts in the list can be priced by the hg model.
*/
class CalibratorCDSRecovery : public pricing::CalibratorGeneral
{

public:

  CalibratorCDSRecovery(finance::CalibrationProcess& calibrationProcess)
      : pricing::CalibratorGeneral(calibrationProcess),
        m_bUseSpreads(true)
  { }

  /**
    Calibrate the specified derivatives using the specified process structure.

    @param pTranslator Translator from underlying process to parameter arrays
    @param derivatives The collection of derivatives for price matching

    @return the calibrated model
  */
  shared_ptr<UnderlyingProcess>
  Calibrate(Translator* pTranslator, const finance::Derivatives& derivatives);

  /**
     Objective function, for a minimization routine.

     @param pdX the current guess for the minimum point
     @param dObjectif (output) the value of the objective function
     @param pdGradients (output) the gradient of the objective function
     @param bComputeGradients whether or not to compute the gradient
   */
  void DoComputeObjectif(const double *pdX, double& dObjectif,
                         double* pdGradients, bool bComputeGradient);

  /**
    Return the last calibrated underlying process (or the last best guess).

    Useful if the calibration fails, but the best calibrated data
    is still desired.

    @return the last calibrated model
   */
  shared_ptr<UnderlyingProcess> GetLastCalibratedProcess() const;


  /**
    Return the last calibrated CDS recovery value

    @return the last calibrated cds recovery value
   */
  double GetCDSRecovery() const;

  /**
    Bracket operator needed by simulated annealing.

    @param pdParam parameter list
    @dF Computed objective function value (output)
   */
  void operator () (const std::vector<double> &pdParam, double &dF);

  /**
    For eds and cds, base objective function of spreads or price

    @param bUseSpreads if true, use spreads.  Otherwise use prices.
  */
  void CalibrateWithSpreads(bool bUseSpreads)
  {
    m_bUseSpreads = bUseSpreads;
  }

protected:

  /**
     Run the underlying calibration code (NAG, simulated annealing, etc).

     @param bUserGradient does the objective function compute gradients
   */
  void RunCalibrator(bool bUserGradients);


  /**
    Initialize the calibration arrays (initial guess, bounds, etc)

    @param pTranslator the translator to get model params
   */
  void InitializeArrays(const Translator* pTranslator);

  /** 
     Add relative tolerance to the weights.

     Call the base class version, but also use spreads instead of prices
     for cds and eds.
   */
  void AddRelTolToWeights();

  /** 
     Objective function calculations for a CDS.
   */
  void ObjectiveForCDS(const double *pdX, double& dObjectif,
                       double* pdGradients, bool bComputeGradient,
                       shared_ptr<finance::Derivative>& pDerivative,
                       double dWeight,
                       shared_ptr<TheoreticalModel>& pModel);

  /** 
     Objective function calculations for a CDS.
   */
  void ObjectiveForEDS(const double *pdX, double& dObjectif,
                       double* pdGradients, bool bComputeGradient,
                       shared_ptr<finance::Derivative>& pDerivative,
                       double dWeight,
                       shared_ptr<TheoreticalModel>& pModel);

  
  /// The translator from parameter array into underlying process 
  Translator* m_pTranslator;

  /// The (constant) underlying process used for mesh construction
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcessForMesh;

  /// For cds/eds, calibrate to spreads or price
  bool m_bUseSpreads;

}; // class CalibratorCDSRecovery


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CALIBRATOR_CDSRECOVERY_H_
