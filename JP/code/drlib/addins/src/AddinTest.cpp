//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLAddin.hpp
//
//   Description : Class for creating and running regression tests for addins
//
//   Author      : Mark A Robson
//
//   Date        : 25 Feb 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/AddinTest.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/PrivateObject.hpp"

DRLIB_BEGIN_NAMESPACE

/** create an instance of an AddinTest - provides ability to write out
    regression file. Note does not take a copy of params */
AddinTest::AddinTest(const string&          addinName,
                     const IObjectSP&       params):
    CObject(TYPE),
    addinName(addinName), params(params){}


/** Run regression test for this addin */
IObjectSP AddinTest::runTest() const{
    static const string routine("AddinTest::run");
    /* need to look up addin name in list of all addins, and get the method
       as well as the return type (which tells us the signature of the
       function we need to call) */
    jvalue  output; // holds all sorts of return values

    try{
        const Addin* addin = Addin::lookUp(addinName);
        bool returnIsNative  = addin->isReturnNative();
        Addin::Method method = addin->getMethod();

        // validate that we've got the right type of object
        CClassConstSP clazz = addin->getDataClass();
        if ((!clazz->isInstance(params)) &&
            (method.objMethod || !IPrivateObject::TYPE->isInstance(params))){
            // Extra check for public->private addin functions
            // eg CASH_SWAP_CURVE. If type is different, then need to check that
            // method.objMethod is non null and
            // check that the current object is not a private type.
            throw ModelException(routine,"Addin has parameters represented by "+
                                 clazz->getName()+" but type "+
                                 params->getClass()->getName()+
                                 " supplied");
        }

        // we can reuse (the non XL specific part of) XLConvert to invoke
        // the method
        const XLConvert* convertMethod;
        if (returnIsNative){
            CClassConstSP  returnType = addin->getReturnType();
            convertMethod = XLConvertFactory::create(returnType);
        } else {
            convertMethod = XLConvertFactory::create(IObject::TYPE);
        }
        convertMethod->invoke(method, params, output);
        if (returnIsNative){
            // need to turn eg a double into an object
            convertMethod->nativeToObject(output);
        }
    } catch (exception& e){
        throw ModelException(e, routine, "Failed whilst executing "
                             "regression test for "+addinName);
    }
    return output.l;
}

// for reflection
AddinTest::AddinTest(): CObject(TYPE){}

class AddinTestHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AddinTest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddinTest);
        FIELD(addinName, "addin name");
        FIELD(params, "Class containing parameters for addin function");
    }

    static IObject* defaultAddinTest(){
        return new AddinTest();
    }
};

CClassConstSP const AddinTest::TYPE = CClass::registerClassLoadMethod(
    "AddinTest", typeid(AddinTest), AddinTestHelper::load);

DRLIB_END_NAMESPACE

