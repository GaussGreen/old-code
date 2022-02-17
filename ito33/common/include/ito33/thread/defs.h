/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/defs.h
// Purpose:     various thread-related declarations
// Author:      Vadim Zeitlin
// Created:     2005-10-07
// RCS-ID:      $Id: defs.h,v 1.1 2005/10/10 12:22:16 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/thread/defs.h
    @brief Common declarations for all thread-related classes.
 */

#ifndef _ITO33_THREAD_DEFS_H_
#define _ITO33_THREAD_DEFS_H_

#ifdef _WIN32
  #include "ito33/win32/winwrap.h"

  // most Windows compilers provide a beginthread() or similar function which
  // not only creates a new thread but also initializes C RTL so that it may
  // be used safely in the new thread -- for this reason, using beginthread()
  // is preferred, but as not compilers have it we can't do it
  // unconditionally hence the tests
  #ifdef _MSC_VER // add more compilers supporting _beginthreadex() here
    // _MT must be defined for multithreaded programs in VC++
    #ifndef _MT
      #error "Please add /D_MT to the project settings and rebuild."
    #endif

    /**
        This macro indicates that we have _beginthreadex() and not only
        (more standard but more buggy and dangerous as well)
        beginthread().

        If the compiler used supports it, you should add code here to
        define it for it as well -- currently this is only defined for
        VC++.
     */
    #define HAVE__BEGINTHREADEX
  #endif // VC++

  // _beginthreadex() and CreateThread() have different types for start
  // routine parameter, define symbolic names to hide these differences
  #ifdef HAVE__BEGINTHREADEX
    /// the return type of the thread entry function
    typedef unsigned THREAD_RETTYPE;

    /// the type of the thread id used with _beginthreadex()/CreateThread()
    typedef unsigned THREAD_IDTYPE;

    /// the calling convention for the thread entry function
    #define THREAD_CALLCONV __stdcall
  #else   // use CreateThread()
    typedef DWORD THREAD_RETTYPE;
    typedef DWORD THREAD_IDTYPE;
    #define THREAD_CALLCONV WINAPI
  #endif
#else // !_WIN32
  #include <unistd.h>         // for usleep()
  #include <pthread.h>

  /// the return type of the thread function
  typedef void *THREAD_RETTYPE;

  /// the type of the thread id
  typedef pthread_t THREAD_IDTYPE;

  /// calling convention of the thread function
  #define THREAD_CALLCONV
#endif // _WIN32/!_WIN32

#include "ito33/debug.h"

namespace ito33
{

/**
    Classes used in multi-thread programming.
 */
namespace thread
{

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_DEFS_H_

