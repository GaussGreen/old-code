
/** Null exists so we can identify the non-existence of an object when 
converting to/from a DataDictionary
*/
#include "edginc/config.hpp"
#include "edginc/Null.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

/** Create an instance of CNull - forcing clients to go through this
    factory method allows us to reuse the same object for all clients. It must
    not be used before this class has been loaded */
IObjectSP CNull::create(){
    return commonCNullObject;
}

IObject* CNull::clone() const{
    return const_cast<CNull*>(this);
}

void CNull::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CNull, clazz);
    SUPERCLASS(CObject);
    // create our static copy
    commonCNullObject.reset(new CNull());
    // no fields
}

CClassConstSP const CNull::TYPE = CClass::registerClassLoadMethod(
    "Null", typeid(CNull), load);

CNull::CNull(): CObject(TYPE){}
// static field definition
IObjectSP CNull::commonCNullObject;

DRLIB_END_NAMESPACE
