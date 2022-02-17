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

#include <stdarg.h>
#include <stdio.h>

#include "Legacy.h"
#undef GtoErrMsg
#include "cerror.h"  // GtoErrMsg


static const size_t bufferSize = 2048;

GTO_EXPORT(extern void )  SafeGtoErrMsg(
    const char *format,    /* (I) printf style format string. */
    ...)                   /* (I) Variable arguments. */
{
    char formattedString[bufferSize];
    int status;
    va_list args;
    va_start(args, format);


    status = vsprintf(formattedString, format, args);

    if (status>0)
    {
        GtoErrMsg((char*)"%s", formattedString);
    }
    else
    {
        GtoErrMsg((char*)"Message too long to print out.\n");
    }
}

