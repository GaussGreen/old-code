/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/crashrpt.h
// Purpose:     generate crash reports containing call stack traces
// Author:      Vadim Zeitlin
// Modified by:
// Created:     19.07.03
// RCS-ID:      $Id: crashrpt.h,v 1.5 2004/04/28 16:15:03 zeitlin Exp $
// Copyright:   (c) 2003 TT-Solutions
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/win32/crashrpt.h
    @brief Declares the class generating crash reports under Win32.

    This class is normally not used directly, the user code should only use
    the ITO33_HANDLE_EXCEPTION() and related macros from ito33/crashrpt.h.
 */

#ifndef _ITO33_WIN32_CRASHRPT_H_
#define _ITO33_WIN32_CRASHRPT_H_

#include "ito33/win32/winwrap.h"

namespace ito33
{

/// report generation flags
enum
{
  /// minimal crash report: just the call stack (this is always included)
  CRASH_REPORT_MINIMAL = 0,

  /// also include the global variable (data segment)
  CRASH_REPORT_WITH_DATA = 1,

  /// include everything: warning, the crash file may be HUGE!
  CRASH_REPORT_ALL = 2,

  /// default flags combination: generate the smallest possible dump file
  CRASH_REPORT_DEFAULT = CRASH_REPORT_MINIMAL
};

/**
    The wrapper class which sets up crash reporting.
 */
class CrashReport
{
public:
  CrashReport();
  ~CrashReport();

  /**
      Generate the crash report in a file.

      This can be called from anywhere.
   */
  bool Generate(int flags = CRASH_REPORT_DEFAULT);

  /**
      Generate the crash report in a file.

      This can only be called from exception handler which has the
      PEXCEPTION_POINTERS.

      @param pExceptionInfo the pointer to exception information
      @param flags report generation flags
      @return true if report has been generated successfully, false otherwise
   */
  bool Generate(PEXCEPTION_POINTERS pExceptionInfo,
         int flags = CRASH_REPORT_DEFAULT);

private:
  // function called when an unhandled exception occurs
  static LONG WINAPI
    UnhandledExceptionFilter(PEXCEPTION_POINTERS pExceptionInfo);

  // the name of the file we're goign to write crash report to
  std::string m_strFileName;

  // the exception filter we replaced
  LONG (WINAPI *m_pFilterOld)(PEXCEPTION_POINTERS);
};

/// Global crash reporter, shouldn't be used directly
extern CrashReport g_crashReporter;

} // namespace ito33

#endif // _ITO33_WIN32_CRASHRPT_H_

