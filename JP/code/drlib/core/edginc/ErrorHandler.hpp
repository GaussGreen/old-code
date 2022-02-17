//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ErrorHandler.hpp
//
//   Description : Prototype error handler
//
//   Author      : Andrew J Swain
//
//   Date        : 8 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_ERRORHANDLER_HPP
#define QLIB_ERRORHANDLER_HPP

// #include "edginc/AtomicArray.hpp"
#include <string>
#include <cstdio>
#include <vector>

CORE_BEGIN_NAMESPACE

using namespace std;

class CORE_DLL ErrorHandler {
public:
    /** Sets the error handler to write messages to the given file. If
        truncate is true the initial contents of the file are lost.
        Does not throw an exception if the file cannot be truncated */
    static void logFile(const string& fileName, bool truncate = false);

    /** Writes the supplied string to the error file. NB Returns error
    status as bool, true = success. Does not throw an exception (avoid
    exception being thrown when in an exception already) */
    static bool writeMsg(const string &msg);

    /** invokes flush on stream where messages are being sent */
    static void flush();

    /** turn error handler on */
    static void on();

    /** turn error handler off */
    static void off();

    /** Creates a default error handle - does not effect the current
        error handling */
    ErrorHandler();

    /** Creates an error handle - does not effect the current error handling */
    ErrorHandler(const string& fileName, bool truncate = false);

    /** Creates an error handle where error messages are appended to
        the supplied string array - does not effect the current error
        handling */
    ErrorHandler(vector<string>&  errorStack);

    /** Sets the current error handler in use. Returns the current
        one, If closeExisting is true the current error handler is
        closed. Does not throw an exception if any IO problems are
        encountered */
    static ErrorHandler set(const ErrorHandler& handler, bool closeExisting);
private:
    /// static members
    static ErrorHandler* handler;

    // fields
    string     fileName;
    bool       truncate;
    bool       useFile;
    bool       active;
    FILE*      errfile;
    // non null => write messages to this stack
    vector<string>* stack;
};

CORE_END_NAMESPACE

#endif
