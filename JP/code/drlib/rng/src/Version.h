// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// SCException.h: Implementation of exception handling class for the
//              Supercube library.
// Author: M.Huq, Credit DR
// -----------------------------------------------------------------------
#include "edginc/coreConfig.hpp"
CORE_BEGIN_NAMESPACE
char *myVersion="$Name$";
char *getSCVersion()
{
    return &myVersion[7];
}// char *getSCVersion()
CORE_END_NAMESPACE
