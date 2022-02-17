// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
//
// SCErrorLog.cpp
//
// Author: M.Huq
// Methods for SCErrorLog class.
// -----------------------------------------------------------------------

#include "edginc/coreConfig.hpp"
#include "edginc/SCException.h"
#include "edginc/SCErrorLog.h"

#include <time.h>

CORE_BEGIN_NAMESPACE

using namespace std;
//---------------------------------------------
// Default constructor
// Uses default log file name from header file.
ErrorLog::ErrorLog()
{
    dateString = (char *)0;
    sprintf(defaultLogFileName,"supercube.log");
    // Just create an instance. Do not initialize.
}//ErrorLog::ErrorLog

//---------------------------------------------
// Users passes in new filename for log file
ErrorLog::ErrorLog(const char *fileName)
{
    initialize(fileName);
} //ErrorLog::ErrorLog

//---------------------------------------------
// User passes in a stream to which the error
// messages should be sent.
ErrorLog::ErrorLog( ostream *outputStream)
{
    initialize(outputStream);
} //ErrorLog::ErrorLog

//---------------------------------------------
// Destructor
ErrorLog::~ErrorLog()
{
    if(errorLogFileStream) {
        errorLogFileStream.close();
    }
} // ErrorLog::~ErrorLog

//---------------------------------------------
// Various initialization methods
void ErrorLog::initialize(const char *fileName)
{
    sprintf(defaultLogFileName,"supercube.log");
    // First check that the filename is a defined string
    if(strlen(fileName)==0) {
        throw SCException(__FILE__, __LINE__,
                          "Invalid filename for log file!");
    }
    // First step: Try to create the log file as specified. Relative
    // to whatever path specified in the filename.
    // Mechanism set up so that if we cannot log to file the default
    // becomes standard error.
    errorLogFileStream.open(fileName,ios::app);
    if(!errorLogFileStream.is_open()) {
        cerr <<"ErrorLog:initialize(char)" << endl;
        cerr <<"       Could not open log file " << fileName << endl;
        cerr <<"       Using default log file "
        <<defaultLogFileName << endl;
        errorLogFileStream.close();
        errorLogFileStream.open(defaultLogFileName,ios::app);
        if(!errorLogFileStream) {
            errorLogStream = &cerr;
            cerr << "ErrorLog: Logging to standard error! Could not open log file"
            << endl;
        } else {
            errorLogStream = &errorLogFileStream;
        }
    } else {
        errorLogStream = &errorLogFileStream;
    }
    LogMessage("****Logging start. <N> = NORMAL, <F> = FAILURE.",NORMAL);
} // ErrorLog::initialize(char)
void ErrorLog::initialize( ostream *theStream)
{
    sprintf(defaultLogFileName,"supercube.log");
    // Check that theStream is defined
    if(!theStream) {
        throw SCException(__FILE__, __LINE__,
                          "Undefined ostream pointer passed in for error logging");
    }
    // Attach given pointer to errorLogStream
    errorLogStream = theStream;
    LogMessage("****Logging start. <N> = NORMAL, <F> = FAILURE",NORMAL);
}
//---------------------------------------------
// Method to obtain a date string or timestamp
char * ErrorLog::getTimeStamp()
{
    time_t t;
    t = time(NULL);
    char *theDate = asctime(localtime(&t));
    return theDate;
} // ErrorLog::getTimeStamp

//---------------------------------------------
// Method to log a message
void ErrorLog::LogMessage(const char *str, const char status)
{
    string dateReturned = getTimeStamp();
    if(*errorLogStream) {
        *errorLogStream << dateReturned.data()
        << " : "
        << "<" << status << "> "
        << str
        << endl;
    }
    // If a failure then log to standard error as well
    if(status == 'F' && errorLogStream != &cerr) {
        cerr <<  dateReturned.data()
        << " : "
        << "< FAILURE > "
        << str
        << endl;
    }
} // ErrorLog::LogMessage

//---------------------------------------------
// Method to flush the output buffers
void ErrorLog::flush()
{
    (*errorLogStream).flush();
}// ErrorLog::flush

CORE_END_NAMESPACE
