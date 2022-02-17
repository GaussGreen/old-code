
#include "edginc/config.hpp"
#include "edginc/ScalarPerNameShift.hpp"

DRLIB_BEGIN_NAMESPACE

IScalarPerNameShift::IScalarPerNameShift(){}
IScalarPerNameShift::~IScalarPerNameShift(){}

void IScalarPerNameShift::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IScalarPerNameShift, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IScalarPerNameShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IScalarPerNameShift", typeid(IScalarPerNameShift), load);

DRLIB_END_NAMESPACE
