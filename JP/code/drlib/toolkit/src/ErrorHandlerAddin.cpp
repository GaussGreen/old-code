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

#include "edginc/config.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/Addin.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class ErrorHandlerAddin: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    string filename;

    /** the 'addin function' */
    static bool logFile(ErrorHandlerAddin* params){
        static const string routine = "ErrorHandlerAddin::logFile";
        try {
            ErrorHandler::logFile(params->filename, true);
            return true;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ErrorHandlerAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ErrorHandlerAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultErrorHandlerAddin);
        FIELD(filename, "filename");

        Addin::registerClassBoolMethod("ERROR_LOG_FILE",
                                       Addin::UTILITIES,
                                       "sets error log file",
                                       TYPE,
                                       (Addin::BoolMethod*)logFile);
    }

    static IObject* defaultErrorHandlerAddin(){
        return new ErrorHandlerAddin();
    }
};

bool ErrorHandlerAddinLoad() {
    return ErrorHandlerAddin::TYPE != 0;
}

CClassConstSP const ErrorHandlerAddin::TYPE = CClass::registerClassLoadMethod(
        "ErrorHandlerAddin", typeid(ErrorHandlerAddin), load);


class ErrorHandlerToggle: public CObject{
    static CClassConstSP const TYPE;

    bool toggle;

    /** the 'addin function' */
    static bool onAndOff(ErrorHandlerToggle* params){
        static const string routine = "ErrorHandlerToggle::onAndOff";
        try {
            if (params->toggle) {
                ErrorHandler::on();
            }
            else {
                ErrorHandler::off();
            }
            return true;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ErrorHandlerToggle():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ErrorHandlerToggle, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultErrorHandlerToggle);
        FIELD(toggle, "toggle");

        Addin::registerClassBoolMethod("ERROR_LOG",
                                       Addin::UTILITIES,
                                       "turns error log on/off",
                                       TYPE,
                                       (Addin::BoolMethod*)onAndOff);
    }

    static IObject* defaultErrorHandlerToggle(){
        return new ErrorHandlerToggle();
    }
};

CClassConstSP const ErrorHandlerToggle::TYPE = CClass::registerClassLoadMethod(
        "ErrorHandlerToggle", typeid(ErrorHandlerToggle), load);

DRLIB_END_NAMESPACE
