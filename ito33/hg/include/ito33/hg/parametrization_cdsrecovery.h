/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/parametrization_cdsrecovery.h
// Purpose:     HG parametrization to be used for calibration with cds recovery
// Created:     2005/07/24
// RCS-ID:      $Id: parametrization_cdsrecovery.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/hg/parametrization_cdsrecovery.h
   @brief HG parametrization to be used for calibration with cds recovery

   NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   This is not meant to be production code.  It is only to be used
   for internal testing.
 */

#ifndef _ITO33_HG_PARAMETRIZATIONCDSRECOVERY_H_
#define _ITO33_HG_PARAMETRIZATIONCDSRECOVERY_H_

#include "ito33/sharedptr.h"

#include "ito33/dlldecl.h"

#include "ito33/hg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
}

namespace hg
{

  class ITO33_HG_DLLDECL UnderlyingProcess;

/**
   HG parametrization to be used for calibration with cds recovery.
 */
  class ITO33_HG_DLLDECL ParametrizationCDSRecovery : public Parametrization
{
public:

  /**
     The initial guess(also the model, since it will give the structure
     of the jumps). 

     By default, all parameters will be calibrated.

     @param pUP The initial guess of the calibration
   */
  ParametrizationCDSRecovery(shared_ptr<UnderlyingProcess> pUP);

  /// virtual dtor
  virtual ~ParametrizationCDSRecovery() { }


  /**
     Calibrates an underlying process by trying to match the given derivatives.
     
     @param derivatives The collection of derivatives whose prices are to be
                        matched.

     @return The calibrated underlying process.
   */
  shared_ptr<UnderlyingProcess>
  Calibrate(const finance::Derivatives& derivatives);

  /** 
    Get the calibrated cds recovery value

    @return the calibrated cds recovery
  */
  double GetCalibratedCDSRecovery()
  {
    return m_dCDSRecovery;
  }

  /**
    For eds and cds, base objective function of spreads or price

    @param bUseSpreads if true, use spreads.  Otherwise use prices.
  */
  void CalibrateWithSpreads(bool bUseSpreads)
  {
    m_bUseSpreads = bUseSpreads;
  }

private:
  
  /// The calibrated CDS recovery value
  double m_dCDSRecovery;

  /// For cds/eds, calibrate to spreads or price
  bool m_bUseSpreads;

}; // class ParametrizationCDSRecovery


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_PARAMETRIZATION_CDSRECOVERY_H_
