/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/winwrap.h
// Purpose:     header to be included instead of <windows.h>
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: winwrap.h,v 1.12 2006/02/08 18:46:16 cosmin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/win32/winwrap.h
    @brief  Wrapper header around &lt;windows.h&gt;

    Instead of directly including &lt;windows.h&gt; in the code needing Win32
    stuff, you should always include this header. The reason for this is that
    &lt;windows.h&gt; #define's many common names which leads to difficult to
    find clashes later and this header does its best to undo the damage by
    #undef'ining these identifiers and declaring inline functions doing the
    same thing instead.

    There are also some lesser advantages in including this header insteadof
    &lt;windows.h&gt; directly -- it includes a few other useful headers as well.
 */

#ifndef _ITO33_WIN32_WINWRAP_H_
#define _ITO33_WIN32_WINWRAP_H_

#ifndef ITO33MFC
#define WIN32_LEAN_AND_MEAN     // exclude some rarely used stuff
#define NOMINMAX                // don't define min() and max() as macros

#include <windows.h>
#include <windowsx.h>
#else
#include "stdafx.h"
#endif

// ----------------------------------------------------------------------------
// #undef'ine macros which conflict with the identifiers used elsewhere
// ----------------------------------------------------------------------------

// NB: only ANSI mode is currently supported

#ifdef FindWindow
  #undef FindWindow

  inline HWND FindWindow(LPCSTR classname, LPCSTR windowname)
  {
    return FindWindowA(classname, windowname);
  }
#endif // FindWindow


#ifdef GetMessage
  #undef GetMessage

    inline BOOL
    GetMessage(LPMSG lpMsg, HWND hWnd, UINT wMsgFilterMin, UINT wMsgFilterMax)
    {
      return GetMessageA(lpMsg, hWnd, wMsgFilterMin, wMsgFilterMax);
    }
#endif // GetMessage


#ifdef AddMonitor
  #undef AddMonitor

    inline BOOL
    AddMonitor(LPTSTR pName, DWORD Level, LPBYTE pMonitors)
    {
      return AddMonitorA(pName, Level, pMonitors);
    }
#endif // AddMonitor

// add more later if needed...

#endif //  _ITO33_WIN32_WINWRAP_H_

