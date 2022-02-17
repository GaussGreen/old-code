/** Defines interface for the ability for one type to convert itself
    into another type.
    This interface is useful for the classes that implements the ITypeConvert.
    This interface provides the 'revert' method that convert a object into an object
    of a different type depending on the interface type (i.e. PYRAMID, KAPITAL,...)

    example: a DependenceMakerGauss object is reverted into a DependenceMaker object 
*/

#include "edginc/config.hpp"
#include "edginc/Null.hpp"
#include "edginc/PublicObject.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/Addin.hpp"
#include "edginc/RevertTypeConvert.hpp"

DRLIB_BEGIN_NAMESPACE

/** revert type for PYRAMID booking */
string const IRevertTypeConvert::PYRAMID = "PYRAMID";

/* static final value in IRevertTypeConvert class */
CClassConstSP const IRevertTypeConvert::TYPE = 
    CClass::registerInterfaceLoadMethod("IRevertTypeConvert", typeid(IRevertTypeConvert), 0);
IRevertTypeConvert::IRevertTypeConvert(){}
IRevertTypeConvert::~IRevertTypeConvert(){}

/** Addin for reverting an object */
class ObjectRevertAddin: public CObject{
public :
    static CClassConstSP const TYPE;

private:
    /**  parameters that our addin takes */
    string         interfaceName;
    IObjectSP      object;

    /** the 'addin function' - construct object from components */
    IObjectSP clientConvert(){
        static const string routine = "ObjectRevertAddin::createObject";

        if (IRevertTypeConvert::TYPE->isInstance(object)){
            IRevertTypeConvert* convertableObj = DYNAMIC_CAST(IRevertTypeConvert,object.get());
            return convertableObj->revert(interfaceName);
        }
        else {
            return object;
        }
    }

    /** for reflection */
    ObjectRevertAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectRevertAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectRevertAddin);
        FIELD(interfaceName, "Interface (PYRAMID for instance) by which the new type will be used");
        FIELD(object, "The object to revert");
        FIELD_MAKE_OPTIONAL(object);
        Addin::registerObjectMethod(
            "REVERT_INTERFACE",
            Addin::UTILITIES,
            "Convert an object into another object depending the interface required",
            true,
            Addin::returnHandle,
            &ObjectRevertAddin::clientConvert);
    }

    static IObject* defaultObjectRevertAddin(){
        return new ObjectRevertAddin();
    }
    
};

CClassConstSP const ObjectRevertAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectRevertAddin", typeid(ObjectRevertAddin), load);

/* for class loading (avoid having header file) */
bool RevertTypeConvertLoad() {
    return (ObjectRevertAddin::TYPE != 0);
}



DRLIB_END_NAMESPACE
