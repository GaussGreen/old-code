/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/parametrization.h
// Purpose:     Base parametrization class declaration
// Created:     2005/07/25
// RCS-ID:      $Id: parametrization.h,v 1.24 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/parametrization.h
    @brief Base parametrization class declaration
 */

#ifndef _ITO33_FINANCE_PARAMETRIZATION_H_
#define _ITO33_FINANCE_PARAMETRIZATION_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/string.h"

#include "ito33/finance/calibrationcontrol.h"

//--------------------------------------------------------------------------
// forward declarations
//--------------------------------------------------------------------------
namespace ito33
{
  class ITO33_DLLDECL Date;
  struct DebugParameters;

namespace finance
{
  class BasketGoodType;

  class ITO33_DLLDECL TheoreticalModel;
  
  class CalibrationProcess;
  class ModelParametersConsumer;
}

}
// end of forward declarations
//--------------------------------------------------------------------------

namespace ito33
{

namespace finance
{

/**
   Basic class for all parametrizations.

   @rename ParametrizationBase
   @nocreate
   @noexport COM
 */
class ITO33_DLLDECL Parametrization
{
public:

  /// Virtual dtor for base class 
  virtual ~Parametrization();

  /**
      @name Debugging output.

      The parametrization can produce debug output containing the value of all 
      input parameters if EnableDebugOutput(true) is called. The output is 
      dumped to a default location (system-dependent) but can be explicitly 
      changed by calling SetDebugOutputFile().
   */
  //@{

  /**
      Flag telling us whether or not to produce debugging output.

      @property DebugOutput
   */
  void EnableDebugOutput(bool bEnable = true);

  /**
      Returns the current value of "debug output" flag.

      @property DebugOutput
   */
  bool IsDebugOutputEnabled() const;

  /**
      The location of the debug output file.

      If the file name is empty, it reverts to default. If the file name is
      just the basename (i.e. without path components), it is created in the
      default directory, otherwise it is used as is.

      Also calls EnableDebugOutput(true), so it is not necessary to call it in
      addition to this method.
   */
  void SetDebugOutputFile(const std::string& filename);

  /**
      The location of the debug output file.

      If EnableDebugOutput() had been called before, the value set by it is
      returned. Otherwise, the default location of the debug file for this
      platform is returned.

      This method may be called when debug output is enabled or not, i.e.
      independently of the value returne by IsDebugOutputEnabled().

      @return the full path to the which is or would be used for debug output
   */
  virtual std::string GetDebugOutputFile() const = 0;
  
  //@}

#ifndef __CPP2ANY__

  /**
      Unified interface for different models.

      @param basket The calibration basket of the parametrization

      @noexport
   */
  virtual void 
  Calibrate(const BasketGoodType& basket) = 0;

#endif

  /**
      Sets the control pointer so that the calibration process could be 
      cancelled from outside by the control.      
      @see CalibrationControl.

      @param pControl control of calibration process.
    
      @noexport
   */
  void SetCalibrationControl(CalibrationControl* pControl);

  /**
      Gets the theoretical model obtained by the calibration.

      @noexport
   */
  virtual shared_ptr<finance::TheoreticalModel> GetTheoreticalModel() = 0;

protected:

  Parametrization();

  shared_ptr<CalibrationProcess> m_pProcess;

  shared_ptr<DebugParameters> m_pDebug;

}; // class Parametrization


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PARAMETRIZATION_H_
