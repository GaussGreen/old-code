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

#ifndef _SC_SCException__H
#define _SC_SCException__H

#include "edginc/coreConfig.hpp"
#include "edginc/ModelException.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace std;


#if 1

CORE_BEGIN_NAMESPACE
class  RNG_DLL SCException : public ModelException
{
    static string toString(int n)
    {
        ostringstream myStream; //creates an ostringstream object
        myStream << n << flush;
        return(myStream.str()); //returns the string form of the stringstream object
    }
public:
    SCException( const char *file, const int line, string description ) :
            ModelException(string(file) + ":" + toString(line) + ": " + description)
    {}
    virtual void printError()
    {
        cerr << what() << endl;
    }
};

CORE_END_NAMESPACE

#else // Old code


// Static global variable within SuperCube to tie in logging feature to
// exception handling. If printError() is used then this variable is used.
// Basically logs to the log file as well as to standard error.
//const ostream *exceptionStreamPtr = NULL;

class  RNG_DLL SCBasicException : public std::exception
{
protected:
    std::string m_description;
public:
    SCBasicException( const char *file, const int line, string description )
    {
        ostringstream os;
        os << file << ": " << line << ": " << description ;
        m_description = os.str();
    }
    virtual ~SCBasicException() throw()
    {}
    // Overload what
    const char * what () const throw ()
    {
        return m_description.c_str();
    }
    // Method to print error. Utilizes log file as well through SCErrorLog
    virtual void printError()
    {
        cerr << what() << endl;
    }
};

class SCException : public SCBasicException
{
public:
    //------------------------------------------------------------
    // pointer to ostream : useful for logging to file.
    // By default we choose cerr as the stream of choice.
    // User can set the stream to be using to point to a log file.
    //------------------------------------------------------------
    // Constructor
    SCException( const char *file, const int line, string description ) : SCBasicException( file, line, description )
    {}
}
;
#endif







CORE_BEGIN_NAMESPACE


// -----------------------------------------------------------------------
// FileIO exception handling

// const int file_name_size = 160;  // Hopefully no file names with path get this big

class  RNG_DLL SCFileIOException : public SCException
{
    /*protected:
        char fileName[ file_name_size ];*/
public:
    // Constructor
    SCFileIOException( const char *file, const int line, const char *description, const char *thisFile ) :
            SCException( file, line,
                         std::string( ": Error processing ") + thisFile + ": " + description)
    {}
}
;

CORE_END_NAMESPACE
#endif //_SC_SCException__H
