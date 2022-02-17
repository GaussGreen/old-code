/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/parametrization.h
// Purpose:     HG parametrization to be used for calibration
// Created:     2005/05/20
// RCS-ID:      $Id: parametrization.h,v 1.13 2006/08/23 10:29:50 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/parametrization.h
    @brief HG parametrization to be used for calibration
 */

#ifndef _ITO33_HG_PARAMETRIZATION_H_
#define _ITO33_HG_PARAMETRIZATION_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/hg/common.h"
#include "ito33/finance/parametrization.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{
  class ITO33_DLLDECL Derivatives;
  class BasketGoodType;
}

namespace hg
{

  class ITO33_HG_DLLDECL UnderlyingProcess;

/**
    HG parametrization to be used for calibration.
 */
class ITO33_HG_DLLDECL Parametrization : public finance::Parametrization
{
public:

  /**
      The initial guess(also the model, since it will give the structure
      of the jumps). 

      By default, all parameters will be calibrated.

      @param pUP The initial guess of the calibration
   */
  Parametrization(const shared_ptr<UnderlyingProcess>& pUP);

  /// virtual dtor
  virtual ~Parametrization() { }

  /**
      The location of the debug output file.

      If EnableDebugOutput() had been called before, the value set by it is
      returned. Otherwise, the default location of the debug file for this
      platform is returned.

      This method may be called when debug output is enabled or not, i.e.
      independently of the value returne by IsDebugOutputEnabled().

      @return the full path to the which is or would be used for debug output
   */
  virtual std::string GetDebugOutputFile() const;

  /**
      Sets the flags indicating which parameter to calibrate.
     
      The ordering of the parameters: vols, default intensities, jumps 
      from regime1 to regime2 with intensity before amplitude

      @param calibrationFlags Flags incidate which parameter to calibrate
   */
  void SetCalibrationFlags(const std::vector<bool>& calibrationFlags);

  /**
      Calibrates an underlying process by trying to match the given derivatives.
     
      @param derivatives The collection of derivatives whose prices are to be
                         matched.

      @return The calibrated underlying process.
   */
  shared_ptr<UnderlyingProcess>
  Calibrate(const finance::Derivatives& derivatives);
  
  virtual void Calibrate(const finance::BasketGoodType& basket); 

  /**
      Gets the calibrated underlying process even if the calibration failed. 

       @return The calibrated underlying process.
   */
  shared_ptr<UnderlyingProcess> GetCalibratedUnderlyingProcess() const
  {
    return m_pCalibratedUP;
  }

  /**
      Dumps all data of this parametriztion in XML format.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  void Dump(ito33::XML::Tag& tagParent) const;
  

  virtual shared_ptr<finance::TheoreticalModel> GetTheoreticalModel();

protected:

  /// The initial guess and the chosen model
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcess;

  /// The calibration flags for each parameter
  std::vector<bool> m_calibrationFlags;

  /// The calibrated underlying process
  shared_ptr<UnderlyingProcess> m_pCalibratedUP;

}; // class Parametrization


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_PARAMETRIZATION_H_
