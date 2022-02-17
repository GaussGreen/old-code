/////////////////////////////////////////////////////////////////////////////
// Name:        utils/debug.cpp
// Purpose:     various debugging helpers
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: debug.cpp,v 1.15 2006/05/27 20:05:06 zhang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  // 'sprintf' was declared deprecated in VC8. But we can still use it in this
  // file
  #define _CRT_SECURE_NO_DEPRECATE 1
#endif

#include "ito33/common.h"
#include "ito33/debug.h"

// must do this to suppress stupid VC++ 6 warnings about unreferences inline
// functions being removed
#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

using namespace ito33;

#ifndef NDEBUG

#ifdef _WIN32
  #include <windows.h>
#else // Unix
  #include <signal.h>      // for SIGTRAP used by wxTrap()
#endif // platform

#ifdef _MSC_VER
  #include <crtdbg.h>

  static class DebugMemChecker
  {
  public:
    DebugMemChecker()
    {
      _CrtSetDbgFlag(
        _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG) | _CRTDBG_LEAK_CHECK_DF);
    }
  } s_debugMemChecker;

#endif // VC++

#include <stdio.h>


// ============================================================================
// implementation
// ============================================================================

OnAssertFunc ito33::AssertFunction = ito33::OnAssert;

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

void
ito33::OnAssert(const char *filename,
                int line,
                const char *cond,
                const char *msg)
{
  // having this variable allows to ignore the asserts: either by choosing
  // "Cancel" in the dialog under Windows or by changing the value of the
  // variable directly under debugger under the other systems
  static bool s_ignoreAsserts = false;

#ifdef _WIN32
  char buf[4096];

  // in any case show the assert string in the debug window
  //
  // NB: the format of the string is such that you can double click it
  //     in the debug window and immediately jump to the code in VC++ IDE
  sprintf(buf,
      "%.1024s(%d) : assert failure: %.1024s "
      "(condition \"%.1024s\" is false).\r\n",
      filename, line, msg, cond);

  ::OutputDebugString(buf);

  if ( !s_ignoreAsserts )
  {
    // also show the assert to the user and offer him the possibility to
    // break into the debugger
    sprintf(buf,
        "File '%.1024s', line %d\n\n"
        "Condition '%.1024s' is false:\n\n%.1024s\n\n"
        "Press \"Yes\" to debug the program, "
        "\"Cancel\" to ignore all asserts\n"
        "or \"No\" to ignore just this one.",
        filename, line, cond, msg);

    switch ( MessageBox(NULL, buf, "ASSERTION FAILURE",
              MB_YESNOCANCEL | MB_ICONSTOP) )
    {
      case IDNO:
        // skip Trap() call
        return;

      case IDCANCEL:
        s_ignoreAsserts = true;
        break;

      // case IDYES: -- nothing to do
    }
  }
#else // Unix
  fprintf(stderr,
      "%.1024s[%d]: assert failure: %.1024s "
      "(condition \"%.1024s\" is false).\n",
      filename, line, msg, cond);
  fflush(stderr);
#endif // Win32/Unix

  if ( !s_ignoreAsserts )
  {
    Trap();
  }
}

void ito33::Trap()
{
#ifdef _WIN32
  ::DebugBreak();
#else // Unix
  raise(SIGTRAP);
#endif // platform
}

#endif // NDEBUG

// ----------------------------------------------------------------------------
// MemoryUsageSnapshot
// ----------------------------------------------------------------------------

#if defined(_MSC_VER) && !defined(NDEBUG)

MemoryUsageSnapshot::MemoryUsageSnapshot()
{
  _CrtMemState ms;
  _CrtMemCheckpoint(&ms);

  m_blocksCur = ms.lCounts[_NORMAL_BLOCK];
  m_bytesCur = ms.lSizes[_NORMAL_BLOCK];
  m_bytexMax = ms.lHighWaterCount;
}

#else // !VC++ || !Debug

MemoryUsageSnapshot::MemoryUsageSnapshot()
{
  m_blocksCur =
  m_bytesCur =
  m_bytexMax = 0;
}

#endif // VC++ && Debug/everything else

