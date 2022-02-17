// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// SCErrorHandlingInterface.h
//
// Author: M.Huq
// Interface class around error logging as well as exception handling.
// Provides a simple interface for initializing error logging and tying
// exception handling to the error log if so desired.
//------------------------------------------------------------------------
//
#ifndef _SC_SCERRORHANDLINGINTERFACE_H
#define _SC_SCERRORHANDLINGINTERFACE_H

#include "edginc/coreConfig.hpp"
#include "edginc/SCErrorLog.h"
#include "edginc/SCException.h"

CORE_BEGIN_NAMESPACE

extern int InitSCErrorLogging();

extern int InitSCErrorLogging(const char *fileName = "supercube.log");
extern int SCLogMessage(const char *str,const  char status);
CORE_END_NAMESPACE

#endif //_SC_SCERRORHANDLINGINTERFACE_H
