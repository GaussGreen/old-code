// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 6/24/99 Robert Buff
//
// $Header$
//

#if ! defined(_KPLATFORM_)
#define _KPLATFORM_

// Platform.h is the only _header_ file that contains platform-
// dependent code. The following source files may also contain platform-
// dependent code: Platform.cpp, Debug.cpp, Memory.cpp, LowLevelError.cpp.

// The Platform.* files contain stuff like compiler settings, simple
// data types, and the like. The other files are about more complex
// and task-oriented.     

// ************************************************************

// Windows NT 4.0 and above, and Windows 98.

#if defined(_WIN32)

#if ! defined(_CPPRTTI)
    #error "CM: RTTI required (turn on in compiler settings)"
#endif

// It's possible to jump to the debugger directly. The trailing
// ";" needs to be appended.    

#pragma warning( disable : 4786 )


#endif

// ************************************************************

// SUNWspro compiler 5.0

#if defined(sun)


#endif



#endif
