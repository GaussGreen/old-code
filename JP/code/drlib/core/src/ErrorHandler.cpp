//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ErrorHandler.cpp
//
//   Description : Prototype error handler
//
//   Author      : Andrew J Swain
//
//   Date        : 8 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/coreConfig.hpp"
#include "edginc/ErrorHandler.hpp"

CORE_BEGIN_NAMESPACE

using namespace std;


ErrorHandler* ErrorHandler::handler = 0;

static FILE* openFile(const string& fileName, bool truncate) {
    string smode;
    if (truncate){
        smode = "w";
    } else {
        smode = "a";
    }
    // use open to create the file
    return (fopen(fileName.c_str(), smode.c_str()));
}

void ErrorHandler::logFile(const string& fileName, bool truncate) {
    ErrorHandler handler(fileName, truncate);
    set(handler, true);
}

/** Sets the current error handler in use. Returns the current one */
ErrorHandler ErrorHandler::set(
    const ErrorHandler& newHandler,
    bool                closeExisting){

    if (handler){
        if (closeExisting){
            // delete handler->errfile;
            if (handler->errfile){
                fclose(handler->errfile);
            }
            handler->errfile = 0;
        }
    } else{
        handler = new ErrorHandler();
    }
    ErrorHandler   previous = *handler;

    *handler = newHandler;
    if (!handler->errfile && !handler->stack){
        // see if we can open (and if required, truncate) the file
        handler->errfile = openFile(handler->fileName, handler->truncate);

        if (handler->errfile) {
            // now shut it
            fclose(handler->errfile);
            handler->errfile = 0;
        }
    }
    return previous;
}

ErrorHandler::ErrorHandler(const string& fileName, bool truncate):
    fileName(fileName),truncate(truncate),useFile(true),
    active(true),errfile(0), stack(0) {}

/** Creates an error handle where error messages are appended to
    the supplied string array - does not effect the current error
    handling */
ErrorHandler::ErrorHandler(vector<string>&  errorStack):
    fileName(""),truncate(false),useFile(false),
    active(true),errfile(0),stack(&errorStack) {}

ErrorHandler::ErrorHandler():truncate(true),active(true),errfile(0),stack(0){
#ifdef WIN32
    fileName = "c:\\error.log";
#else
    char*  homePtr = getenv("HOME");
    if (!homePtr){
        homePtr = "."; // any other ideas?
    }
    string home(homePtr);
    fileName = home + "/error.log";
#endif
    useFile = true;
}

void ErrorHandler::flush(){
    if (handler->errfile){
        fflush(handler->errfile);
    }
}

void ErrorHandler::on() {
    handler->active = true;
}
void ErrorHandler::off() {
    handler->active = false;
}

/** Writes the supplied string to the error file */
bool ErrorHandler::writeMsg(const string &msg) {
    if (!handler){
        set(ErrorHandler(), true);
    }
    if (handler->active) {
        if (handler->stack){
            handler->stack->push_back(msg+"\n");
        } else if (handler->useFile){
            handler->errfile = openFile(handler->fileName, false);

            if (!handler->errfile) {
                return false; // quiet failure
            }
            // NB can't use fprintf unless we protect any '%' in the string
            fwrite(msg.c_str(), msg.size(), 1, handler->errfile);
            fprintf(handler->errfile,"\n"); // then add new line

            fclose(handler->errfile);
            handler->errfile = 0;
        }
    }
    return true;
}

CORE_END_NAMESPACE
