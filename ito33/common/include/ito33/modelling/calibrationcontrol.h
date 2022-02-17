/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/CalibrationControl.h
// Purpose:     CalibrationControl class declaration
// Created:     2005/08/08
// Author:      ZHANG Yunzhi
// RCS-ID:      $Id: calibrationcontrol.h,v 1.3 2006/06/06 16:44:14 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/CalibrationControl.h
    @brief CalibrationControl class declaration
 */

#ifndef _ITO33_FINANCE_CALIBRATIONCONTROL_H_
#define _ITO33_FINANCE_CALIBRATIONCONTROL_H_

#include "ito33/common.h"

namespace ito33
{

namespace finance
{

/**
    Control of calibration process.
    
    In particular, this class permits to stop a calibration process.
   
    Note XXX a calibration class. In general, it will support two member
    functions Calibrate() and SetCalibrationControl(). Before running the
    calibration, the user should call XXX::SetCalibrationControl() to pass
    the pointer of CalibrationControl to the calibration class XXX. 
    Then he should run Parametrization::Calibrate() in a thread
    and should call CalibrationControl::StopProcess() in another thread when
    he wants to stop the calibration process.
 */
class ITO33_DLLDECL CalibrationControl
{
public:

  /// default constructor
  CalibrationControl() { Reset(); }
  
  /// Stops the calibration process
  void StopProcess() { m_IsCancelled = true; }

  /// Gets the process state: whether it has been requested to be cancelled.
  bool ProcessCancelled() { return m_IsCancelled; }

  /// Resets the control
  void Reset() { m_IsCancelled = false; }

private:

  /// flag to check if calibration process should be cancelled
  bool m_IsCancelled;

}; // class CalibrationControl


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_CALIBRATIONCONTROL_H_
