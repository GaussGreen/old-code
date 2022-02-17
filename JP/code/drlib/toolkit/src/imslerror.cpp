//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : imslerror.cpp
//
//   Description : imsl error handling class
//
//   Date        : 24 Jul 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

ModelException IMSLError::imslException = ModelException();

IMSLError::IMSLError(){
    imsl_error_options(/*IMSL_SET_PRINT, IMSL_WARNING, 0,             // printing off
                       IMSL_SET_PRINT, IMSL_WARNING_IMMEDIATE, 0,   // printing off */
                       IMSL_ERROR_PRINT_PROC, print_proc,
                       0);
    int code = imsl_error_code();
    if (code){
        throw ModelException("IMSLError::IMSLError", 
                             "Failed to set imsl's error print procedure.\n"
                             "Error code " + Format::toString(code));
    }
}

bool IMSLError::isError(){
    return !imslException.empty() ;
}

int IMSLError::getErrorCode(){
    return imsl_error_code();
}

void IMSLError::appendToStack(const exception& e){
    if(imslException.empty()) {
        imslException = ModelException(e);
    } else {
        imslException.addMsg(e.what());
    }
}

void IMSLError::throwExceptionIfError(){
    if (isError()){
        ModelException imslExceptionCopy(imslException);
        imslException = ModelException();
        throw imslExceptionCopy;
    }
}

void IMSLError::print_proc(Imsl_error type,
                           long       code,
                           char*      function_name,
                           char*      message){
    // do not print out warnings (setting the warnings off on top in the 
    // imsl_error_options function does not work for some reason)
    if (type == IMSL_WARNING
        || type == IMSL_WARNING_IMMEDIATE){
        return;
    }
    
    string errormsg = Format::toString(function_name)
                    + Format::toString(": Failed with error code %ld (error type %ld)\n", code, type)
                    + Format::toString(message);

    appendToStack(ModelException(errormsg));
}

IMSLError IMSLErrorInstance;


DRLIB_END_NAMESPACE
