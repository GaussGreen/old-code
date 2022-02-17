/////////////////////////////////////////////////////////////////////////////
// Name:        win32/crashrpt.cpp
// Purpose:     generate crash reports containing call stack traces
// Author:      Vadim Zeitlin
// Modified by:
// Created:     19.07.03
// RCS-ID:      $Id: crashrpt.cpp,v 1.10 2006/06/01 10:02:25 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
    The code in this file is heavily based on Matt Pietrek's column from
    the March 2002 issue of MSDN Magazine.
 */

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  // Allow use of strcpy and localtime with VC8
  #define _CRT_SECURE_NO_DEPRECATE 1
#endif

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/longlong.h"
#include "ito33/dynlib.h"

#include "ito33/string.h"

#include "ito33/win32/seh.h"

#include "ito33/win32/crashrpt.h"
#include <imagehlp.h>

#include <time.h>

#ifdef _MSC_VER
  #include <eh.h>
#endif // VC++

#include <imagehlp.h>
#include "ito33/win32/winwrap.h"

#ifndef _T
  #define _T(x) x
#endif

// we need to determine whether we have the declarations for the function in
// debughlp.dll version 5.81 (at least) and we check for DBHLPAPI to test this
//
// reasons:
//  - VC6 version of imagehlp.h doesn't define it
//  - VC7 one does
//  - testing for compiler version doesn't work as you can install and use
//    the new SDK headers with VC6
//
// in any case, the user may override by defining USE_DBGHELP himself
#ifndef USE_DBGHELP
  #ifdef DBHLPAPI
    #define USE_DBGHELP 1
  #else
    #define USE_DBGHELP 0
  #endif
#endif

using namespace ito33;

// ----------------------------------------------------------------------------
// types of imagehlp.h functions
// ----------------------------------------------------------------------------

#if USE_DBGHELP

typedef BOOL (WINAPI *MiniDumpWriteDump_t)(HANDLE, DWORD, HANDLE,
                                           MINIDUMP_TYPE,
                                           CONST PMINIDUMP_EXCEPTION_INFORMATION,
                                           CONST PMINIDUMP_USER_STREAM_INFORMATION,
                                           CONST PMINIDUMP_CALLBACK_INFORMATION);

#endif // USE_DBGHELP

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

// the real crash report generator
class CrashReportImpl
{
public:
  CrashReportImpl(const TCHAR *filename);

  bool Generate(PEXCEPTION_POINTERS pExceptionInfo, int flags);

  ~CrashReportImpl()
  {
    if ( m_hFile != INVALID_HANDLE_VALUE )
    {
      ::CloseHandle(m_hFile);
    }
  }

private:
  // formatted output to m_hFile
  void Output(const TCHAR *format, ...);

  // output end of line
  void OutputEndl() { Output(_T("\r\n")); }

#if USE_DBGHELP
  // load all the functions we need from dbghelp.dll, return true if all ok
  bool BindDbgHelpFunctions(const DynamicLibrary& dllDbgHelp);

  // dynamically loaded dbghelp.dll functions
  #define DECLARE_SYM_FUNCTION(func) static func ## _t func

  DECLARE_SYM_FUNCTION(MiniDumpWriteDump);

  #undef DECLARE_SYM_FUNCTION
#endif // USE_DBGHELP

  // the handle of the report file
  HANDLE m_hFile;
};


// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

CrashReport ito33::g_crashReporter;

// ============================================================================
// implementation
// ============================================================================

#if USE_DBGHELP

#define DEFINE_SYM_FUNCTION(func) func ## _t CrashReportImpl::func = 0

DEFINE_SYM_FUNCTION(MiniDumpWriteDump);

#undef DEFINE_SYM_FUNCTION

#endif // USE_DBGHELP

// ----------------------------------------------------------------------------
// CrashReportImpl
// ----------------------------------------------------------------------------

CrashReportImpl::CrashReportImpl(const TCHAR *filename)
{
  m_hFile = ::CreateFile
        (
          filename,
          GENERIC_WRITE,
          0,                          // no sharing
          NULL,                       // default security
          CREATE_ALWAYS,
          FILE_FLAG_WRITE_THROUGH,
          NULL                        // no template file
        );
}

void CrashReportImpl::Output(const TCHAR *format, ...)
{
  va_list argptr;
  va_start(argptr, format);

  DWORD cbWritten;

  std::string s = String::PrintfV(format, argptr);
  ::WriteFile
    (
      m_hFile,
      s.c_str(),
      static_cast<DWORD>(s.length()) * sizeof(TCHAR),
      &cbWritten,
      0
    );

  va_end(argptr);
}

#if USE_DBGHELP

bool CrashReportImpl::BindDbgHelpFunctions(const DynamicLibrary& dllDbgHelp)
{
  #define LOAD_SYM_FUNCTION(name)                                         \
    name = (name ## _t) dllDbgHelp.GetSymbol(#name);                      \
    if ( !name )                                                          \
    {                                                                     \
      Output(_T("\r\nFunction ") #name                                    \
         _T("() not found.\r\n"));                                        \
      return false;                                                       \
    }

  LOAD_SYM_FUNCTION(MiniDumpWriteDump);

  #undef LOAD_SYM_FUNCTION

  return true;
}

#endif // USE_DBGHELP

bool CrashReportImpl::Generate(PEXCEPTION_POINTERS pExceptionInfo, int flags)
{
  if ( m_hFile == INVALID_HANDLE_VALUE )
    return false;

#if USE_DBGHELP
  if ( !pExceptionInfo )
    return false;

  DynamicLibrary dllDbgHelp;
  if ( dllDbgHelp.Load(_T("dbghelp.dll"), DynamicLibrary::Load_Verbatim) )
  {
    if ( BindDbgHelpFunctions(dllDbgHelp) )
    {
      MINIDUMP_EXCEPTION_INFORMATION minidumpExcInfo;

      minidumpExcInfo.ThreadId = ::GetCurrentThreadId();
      minidumpExcInfo.ExceptionPointers = pExceptionInfo;
      minidumpExcInfo.ClientPointers = FALSE; // in our own address space

      // do generate the dump
      MINIDUMP_TYPE dumpFlags;
      if ( flags & CRASH_REPORT_ALL )
      {
          dumpFlags = MiniDumpWithFullMemory;
      }
      else if ( flags & CRASH_REPORT_WITH_DATA )
      {
          dumpFlags = MiniDumpWithDataSegs;
      }
      else // minimal dump
      {
          dumpFlags = (MINIDUMP_TYPE)(MiniDumpScanMemory |
                                      MiniDumpWithIndirectlyReferencedMemory);
      }

      if ( !MiniDumpWriteDump
            (
              ::GetCurrentProcess(),
              ::GetCurrentProcessId(),
              m_hFile,                    // file to write to
              dumpFlags,                  // kind of dump to craete
              &minidumpExcInfo,
              NULL,                       // no extra user-defined data
              NULL                        // no callbacks
            ) )
      {
          Output(_T("MiniDumpWriteDump() failed."));

          return false;
      }

      return true;
    }
    else // failed to resolve functions
    {
      Output(_T("\r\nPlease update your dbghelp.dll version, ")
             _T("at least version 5.1 is needed!\r\n")
             _T("(if you already have a new version, please ")
             _T("put it in the same directory where the program is.)\r\n"));
    }
  }
  else // failed to load DLL
  {
    Output(_T("Please install dbghelp.dll available free of charge ")
           _T("from Microsoft to get more detailed crash information!"));
  }

  Output(_T("\r\nLatest dbghelp.dll is available at "
            "http://www.microsoft.com/whdc/ddk/debugging/\r\n"));

#else // !USE_DBGHELP
  Output(_T("Support for crash report generation was not included ")
         _T("in this program version."));
#endif // USE_DBGHELP/!USE_DBGHELP

  return false;
}

bool CrashReport::Generate(int flags)
{
  __try
  {
    ::RaiseException(EXCEPTION_BREAKPOINT, 0, 0, NULL);
  }
  __except ( Generate(GetExceptionInformation(), flags),
       EXCEPTION_EXECUTE_HANDLER )
  {
    return true;
  }

  return false;
}

// ----------------------------------------------------------------------------
// CrashReport
// ----------------------------------------------------------------------------

CrashReport::CrashReport()
{
  m_pFilterOld =
    ::SetUnhandledExceptionFilter(CrashReport::UnhandledExceptionFilter);

  // try to find a place where we can put out report file later
  TCHAR szTempDir[MAX_PATH];
  if ( !::GetTempPath(SIZEOF(szTempDir), szTempDir) )
  {
    // when all else fails...
    strcpy(szTempDir, _T("c:\\"));
  }

  // use PID and date to make the report file name more unique
  std::string basename;
  TCHAR szFileName[MAX_PATH];
  if ( ::GetModuleFileName(NULL, szFileName, SIZEOF(szFileName)) )
  {
    // we need just the base name
    char *pDot = strrchr(szFileName, '.');
    if ( pDot )
      *pDot = 0;

    char *pBaseName = strrchr(szFileName, '\\');
    if ( pBaseName )
      pBaseName++;
    else
      pBaseName = szFileName;

    basename = pBaseName;
  }

  char szTimestamp[32];
  time_t t;
  time(&t);
  tm *pTm = localtime(&t);
  strftime(szTimestamp, SIZEOF(szTimestamp), "%Y%m%d_%H%M%S", pTm);

  m_strFileName = String::Printf
          (
            _T("%s%s_%s_%lu.dmp"),
            szTempDir,
            basename.c_str(),
            szTimestamp,
            ::GetCurrentProcessId()
          );
}

CrashReport::~CrashReport()
{
  (void)::SetUnhandledExceptionFilter(m_pFilterOld);
}

bool CrashReport::Generate(PEXCEPTION_POINTERS pExceptionInfo, int flags)
{
  MessageBox
  (
    NULL,
    String::Printf
    (
      _T("An unexpected error has occured, we will try to generate\n")
      _T("a crash report in the file\n")
      _T("\n")
      _T("            %s\n")
      _T("\n")
      _T("Please press \"Ok\" to continue."),
      m_strFileName.c_str()
    ).c_str(),
    _T("ITO33 Crash Report"),
    MB_OK | MB_ICONEXCLAMATION
  );

  CrashReportImpl impl(m_strFileName.c_str());

  if ( impl.Generate(pExceptionInfo, flags) )
  {
    MessageBox
    (
      NULL,
      String::Printf
      (
        _T("Crash report has been successfully generated in the file\n")
        _T("\n")
        _T("            %s\n")
        _T("\n")
        _T("Please send it to ITO 33 technical support.\n\n")
        _T("Thank you in advance and sorry for the inconvenience!"),
        m_strFileName.c_str()
      ).c_str(),
      _T("ITO33 Crash Report"),
      MB_OK | MB_ICONEXCLAMATION
    );

    return true;
  }

  return false;
}

/* static */
LONG CrashReport::UnhandledExceptionFilter(PEXCEPTION_POINTERS pExceptionInfo)
{
  g_crashReporter.Generate(pExceptionInfo);

  return g_crashReporter.m_pFilterOld
       ? g_crashReporter.m_pFilterOld(pExceptionInfo)
       : EXCEPTION_CONTINUE_SEARCH;
}

// ----------------------------------------------------------------------------
// SEH helpers
// ----------------------------------------------------------------------------

namespace ito33
{

namespace Win32
{

#ifdef _MSC_VER

extern void SEtranslator(unsigned int nCode, struct _EXCEPTION_POINTERS *pEP)
{
  g_crashReporter.Generate(pEP);

  throw Win32::StructuredException(nCode);
}

void InitializeSEH()
{
  // we only want to generate crash reports in release builds, in debug ones
  // we'd rather break into the debugger
#ifdef NDEBUG
    // transform Win32 structured exceptions (SE) into C++ exceptions: this
    // *must* be done with VC++ because "catch ( ... )" catches all
    // exceptions, including the SE with it, and it is impossible to catch the
    // SE separately as we want to do in order to generate a crash report for
    // them
    _set_se_translator(Win32::SEtranslator);
#endif // NDEBUG
}

#else // !VC++

void InitializeSEH()
{
  // nothing to do
}

#endif // VC++/!VC++

} // namespace Win32

} // namespace ito33

