/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/dllmain_impl.cpp
// Purpose:     support for COM DLLs
// Author:      Vadim Zeitlin
// Created:     19.01.03
// RCS-ID:      $Id: dllmain_impl.cpp,v 1.17 2006/05/15 15:01:23 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

// logging should be explicitly enabled before compiling this file
#ifndef ITO33_USE_LOG
  #define ITO33_NO_LOG
#endif

#include "ito33/log.h"
#include "ito33/atexit.h"

#include "ito33/win32/module.h"
#include "ito33/win32/seh.h"

using namespace ito33;

ITO33_DEFINE_LOG_CATEGORY(LogDLL, "dll");

#define DLL_TRACE   ITO33_TRACE_CATEGORY(LogDLL)

ITO33_DEFINE_RUN_AT_EXIT();

// include the standard COM DLL functions implementations here (they used to be
// in this file and including them like this avoids the need to modify all the
// existing projects)
#include "comdll.cpp"

// ============================================================================
// DllMain implementation: the DLL entry point
// ============================================================================

BOOL APIENTRY
DllMain(HMODULE hModule, DWORD dwReason, void* /* lpReserved */)
{
#ifndef ITO33_NO_LOG
  static log::Sink *s_logSink = NULL;
#endif 

  switch ( dwReason )
  {
    case DLL_PROCESS_ATTACH:
      {
#ifndef ITO33_NO_LOG
        s_logSink = new log::DebugSink;
        s_logSink->SubscribeAll();
#endif // !ITO33_NO_LOG
        DLL_TRACE("process attaching");

        // store the DLL handle: we may need it for GetModuleFileName() &c
        Win32::Module::SetHandle(hModule);

#ifdef NDEBUG
        Win32::InitializeSEH();
#endif // !debug

        DLL_TRACE("process attached");
      }
      break;

    case DLL_PROCESS_DETACH:
      {
        DLL_TRACE("process detaching");

        ITO33_RUN_AT_EXIT();

        DLL_TRACE("process detached");

#ifndef ITO33_NO_LOG
        delete s_logSink;
        s_logSink = NULL;
#endif // !ITO33_NO_LOG
      }
      break;

    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
      // nothing to do
      DLL_TRACE("thread %s", dwReason == DLL_THREAD_ATTACH ? "attached"
                                                           : "detached");
  }

  return TRUE;
}

