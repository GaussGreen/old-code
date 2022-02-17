// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// SCErrorLog.h
//
// Author: M.Huq
// Class to implement own form of error logging.
//
#ifndef _SC_SCErrorLog__H
#define _SC_SCErrorLog__H

#include "edginc/coreConfig.hpp"

#include <fstream>
#include <string>

CORE_BEGIN_NAMESPACE

using namespace std;

// Basic error logging class implementation
// Will attempt to open a file for logging.
// If no filename is specified it tries to open
// a file in the current directory.
// If it cannot it tries to open one in the users
// home directory in UNIX and C:\supercube.log on NT
// sends to standard error.
//
static const char FAILURE = 'F';
static const char NORMAL = 'N';

class  RNG_DLL ErrorLog
{
protected:
    char logFileName[256];   // Filename for log file
    char defaultLogFileName[256]; // default name
    char *dateString;     // Date string for time stamp
    ostream *errorLogStream; // Output stream for error log
    ofstream errorLogFileStream; // Output filestream for error log

    virtual char * getTimeStamp(); // method to get time stamp for logging
public:
    // Constructors
    ErrorLog();  // Default
    ErrorLog(const char *fileName); // Constructor with output file
    ErrorLog( ostream *outputStream); // constructor with given output stream
    // Destructor
    virtual ~ErrorLog();
    // Methods...
    virtual void initialize(const char *fileName);
    virtual ostream * getOutputStream()
    {
        return errorLogStream;
    }
    virtual void initialize( ostream *outputStream);
    virtual void flush(); // Method to flush the stream
    virtual void LogMessage(const char *str, const char status);
    // overload << operators for logging as well.
    //virtual ostream &operator<<(const ErrorLog it,const char *);
    //ostream &operator<<(ErrorLog,const double);
    //ostream &operator<<(ErrorLog,const int);
    //ostream &operator<<(ErrorLog,const long);

};

CORE_END_NAMESPACE

#endif //_SC_SCErrorLog__H
