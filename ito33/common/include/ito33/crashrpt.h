/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/crashrpt.h
// Purpose:     portable wrapper around Win32-only CrashReport class
// Author:      Vadim Zeitlin
// Modified by:
// Created:     19.07.03
// RCS-ID:      $Id: crashrpt.h,v 1.4 2003/11/03 10:32:19 zhang Exp $
// Copyright:   (c) 2003 TT-Solutions
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/crashrpt.h
    @brief Portable macros for generating crash reports.

    Although we can only create crash reports under Win32 right now we don't
    want the code using CrashReport to force to use #ifdefs everywhere so we
    provide the macros in this file instead.
 */
#ifndef _ITO33_CRASHRPT_H_
#define _ITO33_CRASHRPT_H_

#ifndef _WIN32
  // crash report generation is not supported unless under Win32
  #undef ITO33_GENERATE_CRASH_REPORT
#endif

#ifdef ITO33_GENERATE_CRASH_REPORT
  #include "ito33/win32/crashrpt.h"

  /**
      This macro should occur at the start of any function which wants to
      handle the exceptions.

      It may also occur at the start of a block inside a function and be
      nested. It can't occur at global scope however.
   */
  #define ITO33_PREPARE_HANDLE_EXCEPTION __try {

  /**
      For each use of ITO33_PREPARE_HANDLE_EXCEPTION there must be a matching
      occurence of this macro.

      @param rc the return code to be returned from the function if an
                exception occurs
   */
  #define ITO33_HANDLE_EXCEPTION(rc)                                        \
    }                                                                     \
    __except ( g_crashReporter.Generate(GetExceptionInformation()),       \
         EXCEPTION_EXECUTE_HANDLER )                                \
    {                                                                     \
      return rc;                                                        \
    }

  /**
      Produce the stack back trace for the current location.
   */
  #define ITO33_GENERATE_BACKTRACE() g_crashReporter.Generate()

#else // ITO33_GENERATE_CRASH_REPORT
  #define ITO33_PREPARE_HANDLE_EXCEPTION
  #define ITO33_HANDLE_EXCEPTION(x)
  #define ITO33_GENERATE_BACKTRACE()
#endif // ITO33_GENERATE_CRASH_REPORT/!ITO33_GENERATE_CRASH_REPORT

#endif // _ITO33_CRASHRPT_H_

