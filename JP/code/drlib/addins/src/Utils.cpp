//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Utils.cpp
//
//   Description : Utility Addin functions
//
//   Author      : Mark A Robson
//
//   Date        : 12 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/Handle.hpp"
DRLIB_BEGIN_NAMESPACE

/** class for taking an array of 'raw' strings */
class XLStringArray: public CObject{
public:
    static CClassConstSP const TYPE;

    Handle::RawStringArray   handles;

    XLStringArray(): CObject(TYPE){}

    static bool deleteHandle(XLStringArray* params){
        for (int i = 0; i < params->handles.size(); i++){
            Handle::destroy(params->handles[i]->getString());
        }
        return true;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLStringArray, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLStringArray);
        FIELD(handles, "Handles");
        Addin::registerClassBoolMethod("DELETE_HANDLE",
                                       Addin::UTILITIES,
                                       "Deletes list of handles",
                                       TYPE,
                                       (Addin::BoolMethod*)deleteHandle);

    }

    static IObject* defaultXLStringArray(){
        return new XLStringArray();
    }
};

CClassConstSP const XLStringArray::TYPE = CClass::registerClassLoadMethod(
    "XLStringArray", typeid(XLStringArray), load);

/** Addin to determine run-time type of a hGandle */
class XLIsSimpleTypeAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object;

    /** the 'addin function' - construct object from components */
    static bool isSimpleType(XLIsSimpleTypeAddin* params){
        static const string routine = "XLIsSimpleTypeAddin::isSimpleType";

        // look up how to convert this type of object
        const XLConvert* convertMethod = 
            XLConvertFactory::create(params->object->getClass());

        return(convertMethod->isSimpleType());
    }

    /** for reflection */
    XLIsSimpleTypeAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLIsSimpleTypeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSimpleTypeAddin);
        FIELD(object, "object handle");
        Addin::registerClassBoolMethod("TYPE_IS_SIMPLE",
                                       Addin::UTILITIES,
                                       "determines whether the handle "
                                       "represents a simple type",
                                       TYPE,
                                       (Addin::BoolMethod*)isSimpleType);
    }

    static IObject* defaultSimpleTypeAddin() {
        return new XLIsSimpleTypeAddin();
    }
};

CClassConstSP const XLIsSimpleTypeAddin::TYPE = CClass::registerClassLoadMethod(
    "XLIsSimpleTypeAddin", typeid(XLIsSimpleTypeAddin), load);





// symbol (referenced by Addin.cpp) to ensure file gets linked in
bool UtilsLink = true;


DRLIB_END_NAMESPACE
