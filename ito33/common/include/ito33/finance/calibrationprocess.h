/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/calibrationprocess.h
// Purpose:     Definition of CalibrationProcess class
// Created:     2005/07/25
// RCS-ID:      $Id: calibrationprocess.h,v 1.2 2006/05/30 16:06:34 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/calibrationprocess.h
    @brief Definition of CalibrationProcess class
 */

#ifndef _ITO33_FINANCE_CALIBRATIONPROCESS_H_
#define _ITO33_FINANCE_CALIBRATIONPROCESS_H_

#include "ito33/finance/calibrationcontrol.h"

namespace ito33
{

namespace finance
{
/**
    Wrapper of the pointer to a CalibrationControl object.
   
    In particular, the wrapper hides CalibrationControl::StopProcess()
    that the user uses to stop the calibration process.

    The developer should use this wrapper as he is only allowed to use the
    getter function IsCancelled().
 */
class CalibrationProcess
{
public:
  CalibrationProcess() : m_pControl(0) {}

  void SetControl(CalibrationControl* pControl)
  {
    m_pControl = pControl;
  }

  bool IsCancelled()
  {
    if ( !m_pControl )
      return false;

    return m_pControl->ProcessCancelled();
  }

private:
  CalibrationControl* m_pControl;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_CALIBRATIONPROCESS_H_
