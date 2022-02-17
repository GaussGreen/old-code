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

#include "edginc/SCErrorLog.h"
#include "edginc/SCException.h"

#include <string>
using namespace std;

CORE_BEGIN_NAMESPACE

// An instance of the error logger
ErrorLog theLogger;


int InitSCErrorLogging(const char *fileName)
{
    // Initialize the logger
    try {
        theLogger.initialize(fileName);
        //exceptionStreamPtr = theLogger.getOutputStream();
    } catch (SCException &e) {
        e.printError();
        return 1;
    }
    return 0;
}

int SCLogMessage(const char *str, const char status)
{
    try {
        theLogger.LogMessage(str,status);
    } catch (SCException &e) {
        e.printError();
        return 1;
    }
    return 0;
}

CORE_END_NAMESPACE
