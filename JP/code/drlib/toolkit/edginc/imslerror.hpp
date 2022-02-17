//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : imslerror.hpp
//
//   Description : imsl error handling class
//
//   Date        : 24 Jul 02
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMSL_ERROR_H
#define EDG_IMSL_ERROR_H

#include "edginc/Object.hpp"
#include "edginc/imsl.h"

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL IMSLError{
public:
    IMSLError();

    static bool isError();

    static int getErrorCode();

    // static void appendToStack(const string& errormsg);
    static void appendToStack(const exception& e);

    static void throwExceptionIfError();

private:
    static void print_proc(Imsl_error type,
                           long       code,
                           char*      function_name,
                           char*      message);

    static ModelException imslException;
};

DRLIB_END_NAMESPACE

#endif
