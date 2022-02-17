// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 2/29/00 Bruce Broder
//
// $Header$
//

#if ! defined(_CM_LEGACY_LEGACY_)
#define _CM_LEGACY_LEGACY_

#include "cgeneral.h"  // GTO_EXPORT

#define GtoErrMsg SafeGtoErrMsg;

GTO_EXPORT(extern void) SafeGtoErrMsg(
    const char *format,    /* (I) printf style format string. */
    ...);                  /* (I) Variable arguments. */

#endif
