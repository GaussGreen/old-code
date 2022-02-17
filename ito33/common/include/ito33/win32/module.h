/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/module.h
// Purpose:     global (module-level) stuff used in Win32 programs
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: module.h,v 1.5 2004/10/05 09:13:40 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/win32/module.h
    @brief  Global WIN32 stuff (mostly only used by the library internally).

    Various global functions used by the WIN32 servers are gathered here.
 */

#ifndef _ITO33_WIN32_MODULE_H_
#define _ITO33_WIN32_MODULE_H_

#include "ito33/string.h"

#include "ito33/win32/winwrap.h"

namespace ito33
{

namespace Win32
{

/**
    Module namespace contains the global functions.
 */
namespace Module
{

/**
    Returns the full file name of our EXE or DLL.

    SetHandle() must have been called to initialize the module handle (used by
    this function) before calling it.

    @return the full file name (i.e. including path) of this module
 */
extern std::string GetFileName();

/**
    Returns the global handle (in Win32 sense) of the WIN32 module.

    This global variable is MT-safe because it is set only once at the program
    start (or when the DLL is loaded) using SetHandle(), before any threads
    are created, and remains constant later.

    @return handle to the module or NULL if it hadn't been set
 */
extern HMODULE GetHandle();

/**
    Stores the HMODULE in the internal variable.

    This function must be called exactly once as soon as possible after the
    program start -- this usually means from WinMain or DllMain.

    @param hModule the handle of this module
 */
extern void SetHandle(HMODULE hModule);

} // namespace Module

} // namespace Win32

} // namespace ito33

#endif // _ITO33_WIN32_MODULE_H_

