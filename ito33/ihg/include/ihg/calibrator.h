/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/calibrator.h
// Purpose:     calibration class for ihg
// Author:      David
// Created:     2004/01/14
// RCS-ID:      $Id: calibrator.h,v 1.11 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/calibrator.h
    @brief calibration class

    Implementation of the calibration class for ihg.
 */

#ifndef _IHG_CALIBRATOR_H_
#define _IHG_CALIBRATOR_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33 
{

namespace ihg
{

  class VolatilityTimeOnly;
  class HazardRateTimeOnly;

/**
  @brief Calibration class

  Calibration for the ihg model.  The primary contracts for calibration
  are call options and CDS contracts, since forward equations exist.

  See the articles
  - Tavella, Klopfer, "Implying Local Volatility", Wilmott Magazine, Aug. 2001. 
  - Andersen, Buffum, "Calibration and Implementation of Convertible Bond Models"
    (April 2002 working paper, but also appeared in Journal of Computational
    Finance) for a description of some of the methods used.
*/
class Calibrator
{

public:
  
  // Some useful typedefs
  typedef std::vector< shared_ptr<finance::CDS> > CDSElements;
  typedef std::vector< shared_ptr<finance::Option> > OptionElements;

  Calibrator(shared_ptr<ihg::TheoreticalModel>& pModel)
    : m_pTheoreticalModel(pModel)
  { 
    m_pTheoreticalModel->SetExternalFlagsToDefaults();
  }

  virtual ~Calibrator() { }

  /**
      The main calibration routine for cds and option contracts.

      @param pCDSs List (or vector) of cds contracts
      @param poptions List (or vector) of option contracts
      @param pVol The volatility to use during the calibration
      @param pHR The hazard rate to use during the calibration
  */
  void Calibrate(OptionElements pOptions, CDSElements pCDSs, 
                 shared_ptr<Volatility> pVol, shared_ptr<HazardRate> pHR);

  /**
      Get the calibrated volatility

      @return The calibrated time only volatility
  */
  shared_ptr<VolatilityTimeOnly> GetTimeOnlyVol()
  {
    return m_pVolTimeOnly;
  }

  /**
      Get the calibrated hazard rate

      @return The calibrated time only hazard rate
  */
  shared_ptr<HazardRateTimeOnly> GetTimeOnlyHR()
  {
    return m_pHRTimeOnly;
  }

protected:

  /// Check if the list/vector of CDS contracts is valid for calibration
  void CheckCDS(CDSElements& pCDS);

  /// Check if the list/vector of option contracts is valid for calibration
  void CheckOptions(OptionElements& pOptions);

  /**
     Timeonly hazard rate calibration

     @param pCDS List (vector) of CDS contracts to calibrate
     @param pVol The volatility to use during the calibration
  */
  void CalibrateHR(CDSElements pCDSs, shared_ptr<Volatility> pVol, 
                   shared_ptr<HazardRate> pHR);

  /**
     Timeonly volatility calibration

     @param pOptions List (vector) of option contracts to calibrate
     @param pHR The hazard rate to use during the calibration
  */
  void CalibrateVol(OptionElements pOptions, shared_ptr<HazardRate> pHR,
                    shared_ptr<Volatility> pVol);

  /**
     Timeonly volatility and hazard rate calibration

     @param pCDSs List (or vector) of cds contracts
     @param poptions List (or vector) of option contracts
     @param pVol The volatility to use during the calibration
     @param pHR The hazard rate to use during the calibration
  */
  void CalibrateVolAndHR(OptionElements pOptions, CDSElements pCDSs, 
                 shared_ptr<Volatility> pVol, shared_ptr<HazardRate> pHR);


  /// The theoretical model used for pricing/calibration
  shared_ptr<TheoreticalModel>& m_pTheoreticalModel;

  /// The calibrated time only volatlity
  shared_ptr<VolatilityTimeOnly> m_pVolTimeOnly;

  /// The calibrated time only hazard rate
  shared_ptr<HazardRateTimeOnly> m_pHRTimeOnly;

  NO_COPY_CLASS(Calibrator);

}; // class Calibrator


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CALIBRATOR_H_
